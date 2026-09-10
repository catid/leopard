#!/usr/bin/env python3
"""Run only the fixed untimed tower shape checks, each in its own memory scope."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def run(build, output):
    metadata = json.loads((build/'build.json').read_text())
    require(metadata['completed'] is True and metadata['timed'] is False, 'completed untimed build')
    output.mkdir()
    pins = {build/'source'/name: digest for name, digest in metadata['sources'].items()}
    for profile, record in metadata['profiles'].items():
        pins.update({build/profile/name: digest for name, digest in record['files'].items()})
        pins[Path(record['original_archive'])] = record['original_sha256']
    def verify():
        for path, digest in pins.items():
            require(sha(path) == digest, 'input drift: '+str(path))
    verify()
    results = []
    completed = False
    try:
        for profile in ('release', 'trace', 'sanitize', 'original-release', 'original-sanitize'):
            original = profile.startswith('original-')
            artifact_profile = profile.removeprefix('original-')
            for shape in range(21 if original else 27):
                verify()
                prefix = output/(profile+'-'+str(shape))
                parity = prefix.with_suffix('.parity')
                argv = ['systemd-run', '--user', '--scope', '--expand-environment=no',
                        '-p', 'MemoryMax=256M', '-p', 'MemorySwapMax=0',
                        'bash', '-c', '''set -uo pipefail
"$@"
check_status=$?
resource_group=$(awk -F: '$1 == "0" {print $3}' /proc/self/cgroup)
for name in memory.peak memory.max memory.events memory.swap.current memory.swap.max; do
    printf '%s\\n' "$name"
    cat "/sys/fs/cgroup$resource_group/$name"
done
exit "$check_status"
''', 'tower-check', 'flock', '-n', '/tmp/leopard-gf8-authoritative.lock',
                        'timeout', '--signal=TERM', '--kill-after=5', '60',
                        'prlimit', '--cpu=45:45', '--', 'env', 'OMP_NUM_THREADS=1', 'OMP_DYNAMIC=FALSE',
                        'ASAN_OPTIONS=detect_leaks=1:halt_on_error=1:quarantine_size_mb=8:thread_local_quarantine_size_kb=64',
                        'UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1',
                        str(build/artifact_profile/('test-original' if original else 'test'))]
                argv += (['--shape', str(shape), str(parity)] if shape < 21 else ['--special', str(shape-21)])
                with prefix.with_suffix('.stdout').open('xb') as out, prefix.with_suffix('.stderr').open('xb') as err:
                    status = subprocess.run(argv, stdout=out, stderr=err, timeout=75).returncode
                verify()
                require(status == 0, 'native failure: '+str(prefix))
                lines = prefix.with_suffix('.stdout').read_text().splitlines()
                if original:
                    guard = json.loads(lines.pop(0))
                    require(guard['cell'] == shape and guard['subset_masks'] == 6 and guard['timed'] is False,
                            'original guard record')
                record = json.loads(lines[0])
                require(record.get('shape', record.get('special', -100)+21) == shape and record['timed'] is False, 'native identity')
                require(record['trace'] is (not original and profile != 'release'), 'trace profile')
                require(record.get('original',False) is original, 'original/candidate identity')
                if shape < 21:
                    require(parity.stat().st_size == record['r']*record['bytes'], 'parity size')
                resources = lines[1:]
                expected = ['memory.max', '268435456', 'memory.events', 'low 0', 'high 0',
                            'max 0', 'oom 0', 'oom_kill 0', 'oom_group_kill 0',
                            'memory.swap.current', '0', 'memory.swap.max', '0']
                require(resources[0] == 'memory.peak' and resources[2:] == expected, 'resource counters')
                require(0 < int(resources[1]) < 268435456, 'memory peak')
                results.append(dict(profile=profile, shape=shape, native=record, exit_code=status,
                                    memory_peak=int(resources[1]), parity_sha256=sha(parity) if shape < 21 else None))
                print(profile, shape, record.get('parity_hash', 'special'), flush=True)
        for shape in range(21):
            rows = [r for r in results if r['shape'] == shape]
            require(len(rows) == 5, 'five artifact profiles')
            require(len({r['parity_sha256'] for r in rows}) == 1, 'cross-profile parity')
            require(len({r['native']['scratch_bytes'] for r in rows}) == 1, 'original scratch geometry')
        completed = True
    finally:
        with (output/'checks.json').open('x') as stream:
            json.dump(dict(tracker='leopard-79h.38.5.4.18.4', build=str(build),
                           build_sha256=sha(build/'build.json'), timed=False,
                           completed=completed, results=results), stream, indent=2)


if __name__ == '__main__':
    require(len(sys.argv) == 3, 'usage: run_tower_encoder_checks.py BUILD NEW_OUTPUT')
    run(Path(sys.argv[1]).resolve(strict=True), Path(sys.argv[2]).absolute())

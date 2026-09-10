#!/usr/bin/env python3
"""Fixed cache/selector checks and deliberate unsupported-callback aborts."""
import json
from pathlib import Path
import subprocess
import sys
from run_tower_encoder_checks import require, sha


def run(build, output):
    metadata = json.loads((build/'build.json').read_text())
    require(metadata['completed'] is True and metadata['timed'] is False, 'completed build')
    output.mkdir()
    results = []
    completed = False
    try:
        for profile in ('release', 'trace', 'sanitize'):
            program = build/profile/'test-kernels'
            digest = metadata['profiles'][profile]['files']['test-kernels']
            for case in range(13):
                args = [] if case == 0 else ['--abort',str(case-1)] if case < 10 else [('--measure','--abort','bogus')[case-10]]
                expected = 0 if case == 0 else 134 if case < 10 else 1
                require(sha(program) == digest, 'frozen executable')
                prefix = output/(profile+'-'+str(case))
                command = ['systemd-run','--user','--scope','--expand-environment=no',
                           '-p','MemoryMax=256M','-p','MemorySwapMax=0', 'bash', '-c', '''set -uo pipefail
"$@"
check_status=$?
resource_group=$(awk -F: '$1 == "0" {print $3}' /proc/self/cgroup)
for name in memory.peak memory.max memory.events memory.swap.current memory.swap.max; do
    printf '%s\\n' "$name"
    cat "/sys/fs/cgroup$resource_group/$name"
done
exit "$check_status"
''','tower-kernel-check','flock','-n','/tmp/leopard-gf8-authoritative.lock',
                           'timeout','--signal=TERM','--kill-after=5','60','prlimit','--cpu=45:45','--core=0:0','--',
                           'env','OMP_NUM_THREADS=1','OMP_DYNAMIC=FALSE',
                           'ASAN_OPTIONS=detect_leaks=1:halt_on_error=1:quarantine_size_mb=8:thread_local_quarantine_size_kb=64',
                           'UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1',str(program),*args]
                with prefix.with_suffix('.stdout').open('xb') as out, prefix.with_suffix('.stderr').open('xb') as err:
                    status = subprocess.run(command, stdout=out, stderr=err, timeout=75).returncode
                require(status == expected, 'unexpected exit '+str(prefix)+' '+str(status))
                require(sha(program) == digest, 'post-execution identity')
                lines = prefix.with_suffix('.stdout').read_text().splitlines()
                if case == 0:
                    record = json.loads(lines.pop(0))
                    require(record['selector_cases'] == 9000 and record['log_cases'] == 65536 and
                            record['basis_pairs'] == 2097152 and record['initializations'] == 1 and
                            record['timed'] is False, 'kernel record')
                expected_resources = ['memory.max','268435456','memory.events','low 0','high 0','max 0',
                                      'oom 0','oom_kill 0','oom_group_kill 0', 'memory.swap.current','0','memory.swap.max','0']
                require(lines[0] == 'memory.peak' and 0 < int(lines[1]) < 268435456 and
                        lines[2:] == expected_resources, 'kernel resource counters')
                results.append(dict(profile=profile,case=case,argv=command,exit_code=status,memory_peak=int(lines[1])))
                print(profile,case,status,flush=True)
        completed = True
    finally:
        with (output/'checks.json').open('x') as stream:
            json.dump(dict(tracker='leopard-79h.38.5.4.18.4',completed=completed,timed=False,
                           build=str(build),build_sha256=sha(build/'build.json'),results=results),stream,indent=2)


if __name__ == '__main__':
    require(len(sys.argv) == 3, 'usage: run_tower_kernel_checks.py BUILD NEW_OUTPUT')
    run(Path(sys.argv[1]).resolve(strict=True),Path(sys.argv[2]).absolute())

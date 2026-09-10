#!/usr/bin/env python3
"""Fixed 366-process clock-free qualification, local/serial/capped/locked."""
import json
import os
from pathlib import Path
import subprocess
import sys

from build_paired_epoch import copy_tools
from verify_paired_epoch import (BEAD, inventory, verify_build, verify_record, scope_command,
    sha, compare, split_scope, equal, require)


def check(root):
    out = root/'checks'; out.mkdir()
    source = Path(__file__).resolve().parent
    state = dict(bead=BEAD, root=str(root), completed=False, timed=False, records=[],
                 tools=copy_tools(source, root/'tools', [Path(__file__).name]))
    build, images, _ = verify_build(root)
    env = dict(os.environ, OMP_NUM_THREADS='1', OMP_DYNAMIC='FALSE',
        ASAN_OPTIONS='detect_leaks=1:halt_on_error=1:quarantine_size_mb=8:thread_local_quarantine_size_kb=64',
        UBSAN_OPTIONS='halt_on_error=1:print_stacktrace=1')
    for name in ('LEO_PAIRED_TEST_CLOCK_FAULT', 'LEO_PAIRED_TEST_CLOCK_FAULT_EPOCH'): env.pop(name, None)
    try:
        for row in inventory():
            label, p, kind = row['label'], row['profile'], row['binary']
            args = list(row['args']); parity = out/(label+'.parity') if row['parity'] else None
            if parity: args.append(str(parity))
            executable = root/'build'/p/kind
            equal(sha(executable), build['artifacts'][p+'/'+kind])
            if '--measure' in args: require(kind=='abort' and '-bad-' in label, 'clock-proof CLI refusal only')
            child_env = dict(env)
            if row['fault']:
                child_env['LEO_PAIRED_TEST_CLOCK_FAULT'] = row['fault']
                child_env['LEO_PAIRED_TEST_CLOCK_FAULT_EPOCH'] = str(row['fault_epoch'])
            command = scope_command(root, p, kind, args)
            stdout, stderr = out/(label+'.stdout'), out/(label+'.stderr')
            with stdout.open('xb') as a, stderr.open('xb') as b:
                code = subprocess.run(command, env=child_env, stdout=a, stderr=b, timeout=135).returncode
            record = dict(label=label, args=[p,kind,*args], fault=row['fault'], fault_epoch=row['fault_epoch'],
                returncode=code, stdout_sha256=sha(stdout), stderr_sha256=sha(stderr), scope_command=command)
            state['records'].append(record)
            equal(code, row['code'])
            output, error, peak = split_scope(stdout.read_text(), stderr.read_text(), code)
            record['memory_peak'] = peak
            verify_record(row, output, error, images[p,kind])
            if parity:
                record['parity_sha256'] = sha(parity)
                baseline = out/f'native-{int(row["args"][1])}-NNNN-1-check.parity'
                if parity != baseline: compare(parity, baseline)
                with parity.open('rb') as stream:
                    os.fsync(stream.fileno()); os.posix_fadvise(stream.fileno(), 0, 0, os.POSIX_FADV_DONTNEED)
            equal(sha(executable), build['artifacts'][p+'/'+kind])
            print(label, code, peak, flush=True)
        verify_build(root)
        for name, digest in state['tools'].items():
            equal(sha(root/'tools'/name), digest); equal(sha(source/name), digest)
        state['completed'] = True
    finally: (out/'checks.json').write_text(json.dumps(state, indent=2)+'\n')


if __name__ == '__main__': check(Path(sys.argv[1]).resolve(strict=True))

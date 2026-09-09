#!/usr/bin/env python3
"""Serial 256 MiB/no-swap native qualification. No --measure invocation."""
import json
import os
from pathlib import Path
import subprocess
import sys

from verify_avx2_adjacent_public import BEAD, equal, sha, inventory, verify_record, compare


def check(root):
    out = root/'checks'; out.mkdir()
    state = dict(bead=BEAD,root=str(root),completed=False,timed=False,records=[])
    build = json.loads((root/'build/build.json').read_text())
    equal(build['completed'],True)
    for name,digest in build['artifacts'].items(): equal(sha(root/'build'/name),digest)
    env = dict(os.environ,OMP_NUM_THREADS='1',OMP_DYNAMIC='FALSE',
               ASAN_OPTIONS='detect_leaks=1:quarantine_size_mb=8:thread_local_quarantine_size_kb=64',
               UBSAN_OPTIONS='halt_on_error=1:print_stacktrace=1')
    env.pop('LEO_PAIRED_TEST_CLOCK_FAULT',None)
    try:
        for row in inventory():
            p,binary,label = row['profile'],row['binary'],row['label']
            args = list(row['args'])
            if '--measure' in args: raise ValueError('benchmark invocation forbidden')
            parity = out/(label+'.parity') if row['parity'] else None
            if parity: args.append(str(parity))
            child_env = dict(env)
            if row['fault']: child_env['LEO_PAIRED_TEST_CLOCK_FAULT'] = row['fault']
            stdout,stderr = out/(label+'.stdout'),out/(label+'.stderr')
            with stdout.open('xb') as a, stderr.open('xb') as b:
                result = subprocess.run(['prlimit','--cpu=60:60','--',str(root/'build'/p/binary),*args],
                                        env=child_env,stdout=a,stderr=b,check=False,timeout=120)
            record = dict(label=label,args=[p,binary,*args],fault=row['fault'],returncode=result.returncode,
                          stdout_sha256=sha(stdout),stderr_sha256=sha(stderr))
            state['records'].append(record)
            equal(result.returncode,row['code'])
            verify_record(row,stdout.read_text(),stderr.read_text())
            if parity:
                record['parity_sha256'] = sha(parity)
                cell = int(row['args'][1]); baseline = out/f'native-{cell}-NNNN-1-plain.parity'
                if parity != baseline: compare(parity,baseline)
                with parity.open('rb') as stream:
                    os.fsync(stream.fileno())
                    os.posix_fadvise(stream.fileno(),0,0,os.POSIX_FADV_DONTNEED)
            print(label,result.returncode,flush=True)
        for name,digest in build['artifacts'].items(): equal(sha(root/'build'/name),digest)
        state['completed'] = True
    finally:
        (out/'checks.json').write_text(json.dumps(state,indent=2)+'\n')


if __name__ == '__main__': check(Path(sys.argv[1]).resolve(strict=True))

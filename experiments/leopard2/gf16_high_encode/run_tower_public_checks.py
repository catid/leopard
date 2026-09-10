#!/usr/bin/env python3
"""Run fixed clock-free/synthetic/abort cases serially, each in a 256M scope."""
import json
import os
from pathlib import Path
import subprocess
import sys

from build_tower_public import BEAD, require, sha
from verify_tower_public import inventory, split_scope, verify_record
from verify_avx2_adjacent_public import compare

def check(root):
    out = root/'checks'; out.mkdir()
    state = dict(bead=BEAD,root=str(root),completed=False,timed=False,records=[])
    build = json.loads((root/'build/build.json').read_text())
    require(build['completed'] is True and build['timed'] is False,'qualified build')
    def verify():
        for name,digest in build['artifacts'].items():
            require(sha(root/'build'/name) == digest,'artifact drift')
    verify()
    env = dict(os.environ,OMP_NUM_THREADS='1',OMP_DYNAMIC='FALSE',
        ASAN_OPTIONS='detect_leaks=1:halt_on_error=1:quarantine_size_mb=8:thread_local_quarantine_size_kb=64',
        UBSAN_OPTIONS='halt_on_error=1:print_stacktrace=1')
    env.pop('LEO_PAIRED_TEST_CLOCK_FAULT',None)
    try:
        for row in inventory():
            profile,binary,label = row['profile'],row['binary'],row['label']
            args = list(row['args'])
            require('--measure' not in args,'benchmark invocation forbidden')
            parity = out/(label+'.parity') if row['parity'] else None
            if parity: args.append(str(parity))
            executable = root/'build'/profile/binary
            require(sha(executable)==build['artifacts'][str(executable.relative_to(root/'build'))],
                    'executable before')
            child_env = dict(env)
            if row['fault']: child_env['LEO_PAIRED_TEST_CLOCK_FAULT'] = row['fault']
            argv = ['systemd-run','--user','--scope','--expand-environment=no',
                '-p','MemoryMax=256M','-p','MemorySwapMax=0','bash',str(root/'build/tower_public_scope.sh'),
                'flock','-n','/tmp/leopard-gf8-authoritative.lock',
                'timeout','--signal=TERM','--kill-after=5','120','prlimit','--cpu=60:60','--core=0:0',
                '--',str(executable),*args]
            stdout,stderr = out/(label+'.stdout'),out/(label+'.stderr')
            with stdout.open('xb') as a, stderr.open('xb') as b:
                code = subprocess.run(argv,env=child_env,stdout=a,stderr=b,timeout=135).returncode
            record = dict(label=label,args=[profile,binary,*args],fault=row['fault'],returncode=code,
                stdout_sha256=sha(stdout),stderr_sha256=sha(stderr),scope_command=argv)
            state['records'].append(record)
            require(code==row['code'],'native failure: '+label)
            native,error,peak = split_scope(stdout.read_text(),stderr.read_text(),code)
            record['memory_peak'] = peak
            verify_record(row,native,error)
            if parity:
                record['parity_sha256'] = sha(parity)
                cell = int(row['args'][1]); baseline = out/f'native-{cell}-NNNN-1-plain.parity'
                if parity != baseline: compare(parity,baseline)
                with parity.open('rb') as stream:
                    os.fsync(stream.fileno())
                    os.posix_fadvise(stream.fileno(),0,0,os.POSIX_FADV_DONTNEED)
            require(sha(executable)==build['artifacts'][str(executable.relative_to(root/'build'))],
                    'executable after')
            print(label,code,flush=True)
        verify(); state['completed'] = True
    finally:
        (out/'checks.json').write_text(json.dumps(state,indent=2)+'\n')


if __name__ == '__main__': check(Path(sys.argv[1]).resolve(strict=True))

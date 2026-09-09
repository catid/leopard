#!/usr/bin/env python3
"""Narrow untimed frontend qualification; run in a 256 MiB/no-swap lock scope."""
import json
import os
from pathlib import Path
import subprocess
import sys

from verify_paired_r19932 import ARCHIVES, BEAD, bad_arguments, expected, sha, witness, equal


def check(root):
    output = root/'checks'
    output.mkdir()
    state = dict(bead=BEAD,root=str(root),completed=False,timed=False,records=[])
    build = json.loads((root/'build/build.json').read_text())
    equal(build['completed'],True)
    for name,digest in build['artifacts'].items():
        equal(sha(root/'build'/name),digest)
    env = dict(os.environ, OMP_NUM_THREADS='1', OMP_DYNAMIC='FALSE',
               ASAN_OPTIONS='detect_leaks=1:quarantine_size_mb=8:thread_local_quarantine_size_kb=64',
               UBSAN_OPTIONS='halt_on_error=1:print_stacktrace=1')

    def run(profile,binary,args,label,code,wanted,parity=None):
        argv = ['prlimit','--cpu=60:60','--',str(root/'build'/profile/binary),*args]
        with (output/(label+'.stdout')).open('x') as out, (output/(label+'.stderr')).open('x') as err:
            result = subprocess.run(argv,env=env,stdout=out,stderr=err,check=False)
        record = dict(label=label,args=[profile,binary,*args],returncode=result.returncode,
                      stdout_sha256=sha(output/(label+'.stdout')),stderr_sha256=sha(output/(label+'.stderr')))
        state['records'].append(record)
        equal(result.returncode,code)
        lines = (output/(label+'.stdout')).read_text().splitlines()
        equal([json.loads(line) for line in lines],wanted)
        if parity is not None:
            record['parity_sha256'] = sha(parity)
        print(label, result.returncode, flush=True)

    try:
        for profile in ARCHIVES:
            schedules = ('NNNN',) if profile=='native' else ('0110','1001','0000','1111')
            for cell in range(9):
                for schedule in schedules:
                    for group in ((1,256) if cell==8 else (1,)):
                        for kind,binary in (('check','plain'),('exercise','witness')):
                            label = f'{profile}-{cell}-{schedule}-{group}-{kind}'
                            args = ['--'+kind,str(cell),schedule,str(group)]
                            parity = output/(label+'.parity') if kind=='exercise' else None
                            wanted = [expected(profile,cell,schedule,group,kind=='exercise')]
                            if parity is not None:
                                args.append(str(parity))
                                wanted.append(witness(profile,cell,schedule,group,kind))
                            run(profile,binary,args,label,0,wanted,parity)
            for cell,group in ((0,1),(8,256)):
                for schedule in schedules:
                    label = f'{profile}-{cell}-{schedule}-{group}-clock-guard'
                    run(profile,'witness',['--clock-guard',str(cell),schedule,str(group)],label,86,
                        [witness(profile,cell,schedule,group,'clock-guard')])
        for profile in ARCHIVES:
            for i,args in enumerate(bad_arguments(profile)):
                run(profile,'plain',args,f'{profile}-bad-{i}',1,[])
        for name,digest in build['artifacts'].items():
            equal(sha(root/'build'/name),digest)
        state['completed'] = True
    finally:
        (output/'checks.json').write_text(json.dumps(state,indent=2)+'\n')


if __name__=='__main__':
    check(Path(sys.argv[1]).resolve())

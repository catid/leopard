#!/usr/bin/env python3
"""Untimed native qualification. Never invokes --measure or any real clock."""
import json
import os
from pathlib import Path
import subprocess
import sys
from verify_paired_r19932 import ARCHIVES, BEAD, equal, sha


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

    def run(p,binary,args,label,code,fault=None,parity=None):
        if '--measure' in args: raise ValueError('benchmark invocation forbidden')
        child_env = dict(env)
        if fault: child_env['LEO_PAIRED_TEST_CLOCK_FAULT'] = fault
        with (out/(label+'.stdout')).open('x') as stdout, (out/(label+'.stderr')).open('x') as stderr:
            result = subprocess.run(['prlimit','--cpu=60:60','--',str(root/'build'/p/binary),*args],
                                    env=child_env,stdout=stdout,stderr=stderr,check=False)
        record = dict(label=label,args=[p,binary,*args],fault=fault,returncode=result.returncode,
                      stdout_sha256=sha(out/(label+'.stdout')),stderr_sha256=sha(out/(label+'.stderr')))
        state['records'].append(record)
        equal(result.returncode,code)
        if parity:
            record['parity_sha256'] = sha(parity)
            # Retained evidence need not stay in this job's page cache.
            with parity.open('rb') as stream:
                os.fsync(stream.fileno())
                os.posix_fadvise(stream.fileno(),0,0,os.POSIX_FADV_DONTNEED)
        print(label,result.returncode,flush=True)

    try:
        for p in ARCHIVES:
            schedules = ('NNNN',) if p=='native' else ('0110','1001','0000','1111')
            for c in range(9):
                for s in schedules:
                    for g in ((1,256) if c==8 else (1,)):
                        for variant,mode in (('plain','--exercise'),('synthetic','--clock-exercise')):
                            label = f'{p}-{c}-{s}-{g}-{variant}'
                            parity = out/(label+'.parity')
                            run(p,variant,[mode,str(c),s,str(g),str(parity)],label,0,parity=parity)
            for c,g in ((0,1),(8,256)):
                for s in schedules:
                    run(p,'abort',['--clock-guard',str(c),s,str(g)],f'{p}-{c}-{s}-{g}-abort',86)
            s = schedules[0]
            for fault in ('equal','reverse','negative','huge'):
                run(p,'synthetic',['--clock-exercise','8',s,'256'],f'{p}-fault-{fault}',1,fault=fault)
            bad = [[],['--check'],['--timing','0',s,'1'],['--exercise','0',s,'1','a','b'],
                   ['--check','9',s,'1'],['--check','00',s,'1'],['--check','0',s,'256'],
                   ['--check','8',s,'0256'],['--check','8',s,'0'],['--check','0','0101','1'],
                   ['--check','0',('0110' if p=='native' else 'NNNN'),'1'],
                   ['--clock-exercise','0',s,'1'],['--clock-guard','0',s,'1','forbidden'],
                   ['--clock-guard','8',s,'256']]
            for i,args in enumerate(bad): run(p,'plain',args,f'{p}-bad-{i}',1)
            if p!='native': run(p,'group-unit',[],p+'-unit',0)
        for name,digest in build['artifacts'].items(): equal(sha(root/'build'/name),digest)
        state['completed'] = True
    finally:
        (out/'checks.json').write_text(json.dumps(state,indent=2)+'\n')


if __name__=='__main__': check(Path(sys.argv[1]).resolve())

#!/usr/bin/env python3
"""Focused, untimed runtime matrix; no benchmark frontend is executed."""
import json
import os
from pathlib import Path
import subprocess
import sys
from verify_paired_r19932 import equal,parse,require,sha

BEAD = 'leopard-79h.38.5.4.18.3.1'
SELECTORS = ['--adjacent','--forward-ranges','--pairs','--split',*map(str,range(8)),
             '--roundtrip','--concurrent']


def check(root):
    out = root/'checks'; out.mkdir()
    build = json.loads((root/'build/build.json').read_text())
    equal([build['bead'],build['completed'],build['timed']],[BEAD,True,False])
    pins = {str(root/'build'/p/n):v for p,e in build['profiles'].items() for n,v in e['files'].items()}
    for p,e in build['profiles'].items():
        pins[str(root/'build'/p/'candidate.a')] = e['archive_sha256']
    def verify():
        for name,value in pins.items(): equal(sha(Path(name)),value)
    verify()
    env = dict(os.environ,OMP_NUM_THREADS='1',OMP_DYNAMIC='FALSE',OMP_THREAD_LIMIT='1',
        ASAN_OPTIONS='detect_leaks=1:quarantine_size_mb=8:thread_local_quarantine_size_kb=64:abort_on_error=1',
        UBSAN_OPTIONS='halt_on_error=1:print_stacktrace=1')
    state = dict(bead=BEAD,timed=False,completed=False,records=[],pins=pins)
    def run(profile,label,binary,args,code=0):
        directory = root/'build'/profile
        argv = ['/usr/bin/prlimit','--cpu=30:30','--fsize=1048576:1048576','--',str(directory/binary),*args]
        with (out/(label+'.stdout')).open('xb') as stdout, (out/(label+'.stderr')).open('xb') as stderr:
            child = subprocess.run(argv,stdout=stdout,stderr=stderr,env=env,timeout=60)
        row = dict(label=label,argv=argv,returncode=child.returncode,expected=code,
            stdout_sha256=sha(out/(label+'.stdout')),stderr_sha256=sha(out/(label+'.stderr')))
        state['records'].append(row); equal(child.returncode,code)
        if code==0:
            equal((out/(label+'.stderr')).read_text(),'')
            lines = (out/(label+'.stdout')).read_text().splitlines()
            record = parse(lines[-1])
            if binary=='control-unit':
                equal(record,dict(schema='adjacent-control-unit/v1',trace=profile!='release',timed=False))
            else:
                equal([record['schema'],record['mode'],record['trace'],record['timed']],
                      ['adjacent-runtime-focused/v1',int(args[0]=='on'),profile!='release',False])
                require(len(record['calls'])==len(record['blocks'])==2,'family dimensions')
                for family in range(2):
                    for key in ('calls','blocks'):
                        equal(record[key][family][int(args[0]!='on')],0)
                        if profile=='release': equal(record[key][family],[0,0])
        else:
            require((out/(label+'.stderr')).stat().st_size>0,'missing refusal')
        verify(); print(label+' '+str(child.returncode),flush=True)
    try:
        for profile in ('release','trace','sanitize'):
            run(profile,profile+'-unit','control-unit',[])
            for mode in ('off','on'):
                for selector in SELECTORS:
                    run(profile,profile+'-'+mode+'-'+selector.removeprefix('--'),'focused',[mode,selector])
            for i,args in enumerate(([],['off'],['bad','--pairs'],['on','--pairs','extra'],
                                      ['off','8'],['on','--measure'])):
                run(profile,profile+'-bad-'+str(i),'focused',args,1)
        state['completed'] = True
    finally:
        with (out/'checks.json').open('x') as stream: json.dump(state,stream,indent=2); stream.write('\n')


if __name__=='__main__':
    require(len(sys.argv)==2,'usage: check_avx2_adjacent_runtime.py ROOT')
    check(Path(sys.argv[1]).resolve())

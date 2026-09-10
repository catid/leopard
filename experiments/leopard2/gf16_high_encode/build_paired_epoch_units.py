#!/usr/bin/env python3
"""Build/check supplemental units; each native job is serialized and capped."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from verify_paired_epoch import BEAD, parse, equal, require, sha, scope_command, split_scope
from verify_paired_epoch_units import SOURCE_NAMES, inventory, verify_record


def build(root):
    out=root/'units'; out.mkdir(); source=Path(__file__).resolve().parent
    state=dict(bead=BEAD,completed=False,timed=False,source_sha256={},commands=[],artifacts={})
    frontend=parse((root/'build/build.json').read_text())
    try:
        for name in (*SOURCE_NAMES,Path(__file__).name,'verify_paired_epoch_units.py'):
            shutil.copyfile(source/name,out/name)
            if name in SOURCE_NAMES: state['source_sha256'][name]=sha(out/name)
        for p in ('native','release','sanitize'):
            directory=root/'build'/p
            driver=next(c for c in frontend['commands'] if c[-1]==str(directory/'driver.o'))
            flags=driver[driver.index('c++')+1:driver.index('-c')]
            for kind,filename in zip(('metadata','clock'),SOURCE_NAMES):
                obj=str(out/(p+'-'+kind+'.o')); binary=str(out/(p+'-'+kind))
                compile_command=['prlimit','--cpu=120:120','--','c++',*flags,'-I'+str(root/'build'),
                                 '-c',str(out/filename),'-o',obj]
                if kind=='metadata':
                    link=next(c for c in frontend['commands'] if c[-1]==str(directory/'abort'))
                    link=[obj if a==str(directory/'driver.o') else a for a in link]; link[-1]=binary
                else: link=['prlimit','--cpu=120:120','--','c++',*flags,obj,'-o',binary]
                for command in (compile_command,link):
                    state['commands'].append(command); print(json.dumps(command),flush=True)
                    subprocess.run(command,check=True)
        for name,digest in frontend['artifacts'].items(): equal(sha(root/'build'/name),digest)
        for path in sorted(out.iterdir()):
            state['artifacts'][path.name]=sha(path)
            path.chmod(0o555 if path.suffix=='' else 0o444)
        state['completed']=True
    finally: (out/'build.json').write_text(json.dumps(state,indent=2)+'\n')


def check(root):
    out=root/'units'; build_state=parse((out/'build.json').read_text())
    require(build_state['completed'] is True,'units built')
    state=dict(bead=BEAD,completed=False,timed=False,records=[])
    env=dict(os.environ,OMP_NUM_THREADS='1',OMP_DYNAMIC='FALSE',
        ASAN_OPTIONS='detect_leaks=1:halt_on_error=1:quarantine_size_mb=8:thread_local_quarantine_size_kb=64',
        UBSAN_OPTIONS='halt_on_error=1:print_stacktrace=1')
    for name in ('LEO_PAIRED_TEST_CLOCK_FAULT','LEO_PAIRED_TEST_CLOCK_FAULT_EPOCH'): env.pop(name,None)
    try:
        for row in inventory():
            p,kind,label=row['profile'],row['kind'],row['label']; binary=out/(p+'-'+kind)
            equal(sha(binary),build_state['artifacts'][binary.name])
            command=scope_command(root,p,'abort',row['args'])
            command[len(command)-len(row['args'])-1]=str(binary)
            stdout,stderr=out/(label+'.stdout'),out/(label+'.stderr')
            with stdout.open('xb') as a,stderr.open('xb') as b:
                code=subprocess.run(command,env=env,stdout=a,stderr=b,timeout=135).returncode
            record=dict(label=label,command=command,returncode=code,stdout_sha256=sha(stdout),stderr_sha256=sha(stderr))
            state['records'].append(record); equal(code,row['code'])
            output,error,peak=split_scope(stdout.read_text(),stderr.read_text(),code); record['memory_peak']=peak
            verify_record(row,output,error); equal(sha(binary),build_state['artifacts'][binary.name])
            print(label,code,peak,flush=True)
        state['completed']=True
    finally: (out/'checks.json').write_text(json.dumps(state,indent=2)+'\n')


if __name__=='__main__':
    root=Path(sys.argv[1]).resolve(strict=True)
    if len(sys.argv)==3 and sys.argv[2]=='--check': check(root)
    elif len(sys.argv)==2: build(root)
    else: raise SystemExit('usage: ROOT [--check]')

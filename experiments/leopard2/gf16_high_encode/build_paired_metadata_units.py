#!/usr/bin/env python3
"""Supplemental same-header boundary builds; requires 512MiB serial locked scope."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from paired_metadata_overlay import BEAD
from verify_paired_metadata import parse, split_scope
from verify_paired_r19932 import require, sha, equal


def build(root):
    out=root/'units'; out.mkdir()
    source=Path(__file__).parent/'test_paired_runtime_metadata.cpp'
    shutil.copyfile(source,out/source.name)
    shutil.copyfile(Path(__file__),out/Path(__file__).name)
    state=dict(bead=BEAD,completed=False,timed=False,source_sha256=sha(source),commands=[],artifacts={})
    original=parse((root/'build/build.json').read_text())
    try:
        for p in ('native','release','sanitize'):
            directory=root/'build'/p
            command=next(c for c in original['commands'] if c[-1]==str(directory/'driver.o'))
            flags=command[command.index('c++')+1:command.index('-c')]
            compile_command=['prlimit','--cpu=120:120','--','c++',*flags,'-I'+str(root/'build'),'-c',
                str(out/source.name),'-o',str(out/(p+'.o'))]
            link=next(c for c in original['commands'] if c[-1]==str(directory/'abort'))
            link=[str(out/(p+'.o')) if a==str(directory/'driver.o') else a for a in link]
            link[-1]=str(out/p)
            for argv in (compile_command,link):
                state['commands'].append(argv); print(json.dumps(argv),flush=True)
                subprocess.run(argv,check=True)
        for name,digest in original['artifacts'].items(): require(sha(root/'build'/name)==digest,'original artifact drift')
        for path in sorted(out.iterdir()):
            state['artifacts'][path.name]=sha(path); path.chmod(0o555 if path.name in ('native','release','sanitize') else 0o444)
        state['completed']=True
    finally: (out/'build.json').write_text(json.dumps(state,indent=2)+'\n')


def check(root):
    out=root/'units'
    state=dict(bead=BEAD,completed=False,timed=False,records=[])
    env=dict(os.environ,ASAN_OPTIONS='detect_leaks=1:halt_on_error=1:quarantine_size_mb=8:thread_local_quarantine_size_kb=64',
             UBSAN_OPTIONS='halt_on_error=1:print_stacktrace=1',OMP_NUM_THREADS='1',OMP_DYNAMIC='FALSE')
    env.pop('LEO_PAIRED_TEST_CLOCK_FAULT',None)
    build_state=parse((out/'build.json').read_text()); require(build_state['completed'] is True,'unit build')
    try:
        for profile in ('native','release','sanitize'):
            command=['systemd-run','--user','--scope','--expand-environment=no','-p','MemoryMax=256M','-p','MemorySwapMax=0',
                'bash',str(root/'build/tower_public_scope.sh'),'flock','-n','/tmp/leopard-gf8-authoritative.lock',
                'timeout','--signal=TERM','--kill-after=5','120','prlimit','--cpu=60:60','--core=0:0','--',str(out/profile)]
            stdout,stderr=out/(profile+'.stdout'),out/(profile+'.stderr')
            equal(sha(out/profile),build_state['artifacts'][profile])
            with stdout.open('xb') as a,stderr.open('xb') as b:
                code=subprocess.run(command,env=env,stdout=a,stderr=b,timeout=135).returncode
            record=dict(profile=profile,command=command,returncode=code,
                        stdout_sha256=sha(stdout),stderr_sha256=sha(stderr))
            state['records'].append(record); equal(code,0)
            native,error,peak=split_scope(stdout.read_text(),stderr.read_text(),0)
            record['memory_peak']=peak; equal(error,'')
            equal([parse(s) for s in native.splitlines()],
                [dict(schema='paired-metadata-unit/v1',cases=22,timed=False),
                 dict(schema='leopard-paired-witness/v1',calls=0,states=[0,0,0],apis=[0,0,0],order_hash='cbf29ce484222325')])
            equal(sha(out/profile),build_state['artifacts'][profile])
            print(profile,peak,flush=True)
        state['completed']=True
    finally: (out/'checks.json').write_text(json.dumps(state,indent=2)+'\n')


if __name__=='__main__':
    root=Path(sys.argv[1]).resolve(strict=True)
    if len(sys.argv)==3 and sys.argv[2]=='--check': check(root)
    else: build(root)

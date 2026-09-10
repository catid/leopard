#!/usr/bin/env python3
"""Run fixed untimed/abort/synthetic checks locally and serially. No real clocks."""
import ast
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

from paired_metadata_overlay import BEAD
from verify_paired_metadata import inventory, split_scope, verify_record, elf, parse
from verify_paired_r19932 import sha, require, compare


def check(root):
    out=root/'checks'; out.mkdir()
    source=Path(__file__).resolve().parent
    retained=root/'tools'; retained.mkdir()
    state=dict(bead=BEAD,root=str(root),completed=False,timed=False,records=[],tools={})
    def retain(name):
        if name in state['tools']: return
        path=source/name
        state['tools'][name]=sha(path); shutil.copyfile(path,retained/name)
        (retained/name).chmod(0o444)
        for node in ast.walk(ast.parse(path.read_text())):
            names=([node.module] if isinstance(node,ast.ImportFrom) else
                   [a.name for a in node.names] if isinstance(node,ast.Import) else [])
            for module in names:
                if module and (source/(module+'.py')).is_file(): retain(module+'.py')
    retain(Path(__file__).name)
    build=parse((root/'build/build.json').read_text())
    require(build['completed'] is True and build['timed'] is False,'build complete')
    def verify():
        for name,digest in build['artifacts'].items(): require(sha(root/'build'/name)==digest,'build artifact drift')
        for name,digest in state['tools'].items(): require(sha(retained/name)==digest,'retained verifier drift')
    verify()
    env=dict(os.environ,OMP_NUM_THREADS='1',OMP_DYNAMIC='FALSE',
        ASAN_OPTIONS='detect_leaks=1:halt_on_error=1:quarantine_size_mb=8:thread_local_quarantine_size_kb=64',
        UBSAN_OPTIONS='halt_on_error=1:print_stacktrace=1')
    env.pop('LEO_PAIRED_TEST_CLOCK_FAULT',None)
    images={(p,k):elf(root/'build'/p/k,p=='native') for p in ('native','release','sanitize') for k in ('abort','synthetic')}
    try:
        for row in inventory():
            label,p,kind=row['label'],row['profile'],row['binary']
            args=list(row['args']); parity=out/(label+'.parity') if row['parity'] else None
            if parity: args.append(str(parity))
            executable=root/'build'/p/kind
            require(sha(executable)==build['artifacts'][p+'/'+kind],'executable before')
            # --measure is an intentional pre-allocation refusal test, on an
            # abort-linked binary whose real-clock imports were independently rejected.
            if '--measure' in args: require(kind=='abort' and '-bad-' in label,'only clock-proof refusal')
            child_env=dict(env)
            if row['fault']: child_env['LEO_PAIRED_TEST_CLOCK_FAULT']=row['fault']
            command=['systemd-run','--user','--scope','--expand-environment=no','-p','MemoryMax=256M','-p','MemorySwapMax=0',
                'bash',str(root/'build/tower_public_scope.sh'),'flock','-n','/tmp/leopard-gf8-authoritative.lock',
                'timeout','--signal=TERM','--kill-after=5','120','prlimit','--cpu=60:60','--core=0:0','--',str(executable),*args]
            stdout,stderr=out/(label+'.stdout'),out/(label+'.stderr')
            with stdout.open('xb') as a,stderr.open('xb') as b:
                code=subprocess.run(command,env=child_env,stdout=a,stderr=b,timeout=135).returncode
            record=dict(label=label,args=[p,kind,*args],fault=row['fault'],returncode=code,
                stdout_sha256=sha(stdout),stderr_sha256=sha(stderr),scope_command=command)
            state['records'].append(record)
            require(code==row['code'],'native failure: '+label)
            output,error,peak=split_scope(stdout.read_text(),stderr.read_text(),code)
            record['memory_peak']=peak
            verify_record(row,output,error,images[p,kind])
            if parity:
                record['parity_sha256']=sha(parity)
                cell=int(row['args'][1]); compare(parity,out/f'native-{cell}-NNNN-1-check.parity')
                with parity.open('rb') as stream:
                    os.fsync(stream.fileno()); os.posix_fadvise(stream.fileno(),0,0,os.POSIX_FADV_DONTNEED)
            require(sha(executable)==build['artifacts'][p+'/'+kind],'executable after')
            print(label,code,flush=True)
        verify(); state['completed']=True
    finally:
        (out/'checks.json').write_text(json.dumps(state,indent=2)+'\n')


if __name__=='__main__': check(Path(sys.argv[1]).resolve(strict=True))

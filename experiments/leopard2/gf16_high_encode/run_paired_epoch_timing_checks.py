#!/usr/bin/env python3
"""One bounded clock-free frontend qualification; no real --measure dispatch."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

from build_paired_epoch_timing import copy_tools, REFERENCE
from verify_paired_epoch_timing import (BEAD,inventory,verify_build,verify_record,executable,
    scope_command,sha,parse,equal,require,split_scope,native_reference,NATIVE_CHECKS_PIN)
from verify_paired_metadata import compare


def check(root):
    out = root/'checks'; out.mkdir()
    source = Path(__file__).resolve().parent
    state = dict(bead=BEAD,root=str(root),completed=False,real_clocks_read=False,records=[],
                 tools=copy_tools(source,root/'check_tools',[Path(__file__).name]))
    build,images,_ = verify_build(root)
    reference = REFERENCE.parent/'checks'
    equal(sha(reference/'checks.json'),NATIVE_CHECKS_PIN)
    target = root/'native_reference'; target.mkdir()
    for name in ['checks.json']+[f'native-{c}-NNNN-1-check.parity' for c in range(9)]:
        shutil.copyfile(reference/name,target/name); (target/name).chmod(0o444)
    references = native_reference(root)
    env = dict(os.environ,OMP_NUM_THREADS='1',OMP_DYNAMIC='FALSE',OMP_THREAD_LIMIT='1',
        ASAN_OPTIONS='detect_leaks=1:halt_on_error=1:quarantine_size_mb=8:thread_local_quarantine_size_kb=64',
        UBSAN_OPTIONS='halt_on_error=1:print_stacktrace=1')
    for name in ('LEO_PAIRED_TEST_CLOCK_FAULT','LEO_PAIRED_TEST_CLOCK_FAULT_EPOCH'): env.pop(name,None)
    try:
        rows = inventory()
        for row in rows:
            p,v,label = row['profile'],row['variant'],row['label']
            name = executable(p,v); binary = root/'build'/name
            args = list(row['args'])
            parity = out/(label+'.parity') if row['code']==0 and args[0]!='--measure' else None
            if parity: args.append(str(parity))
            command = scope_command(root,p,v,args)  # fail closed BEFORE dispatch
            equal(sha(binary),build['artifacts'][name])
            # Variant source/recipe proof and actual ELF bindings were checked
            # above. A fake-steady stdout label never licenses a real binary.
            child_env = dict(env)
            if row['fault']:
                child_env['LEO_PAIRED_TEST_CLOCK_FAULT'] = row['fault']
                child_env['LEO_PAIRED_TEST_CLOCK_FAULT_EPOCH'] = str(row['fault_epoch'])
            stdout,stderr = out/(label+'.stdout'),out/(label+'.stderr')
            with stdout.open('xb') as a,stderr.open('xb') as b:
                code = subprocess.run(command,env=child_env,stdout=a,stderr=b,timeout=135).returncode
            record = dict(label=label,variant=v,command=command,executable_sha256=build['artifacts'][name],
                          fault=row['fault'],fault_epoch=row['fault_epoch'],returncode=code,
                          stdout_sha256=sha(stdout),stderr_sha256=sha(stderr))
            state['records'].append(record)
            equal(code,row['code'])
            output,error,peak = split_scope(stdout.read_text(),stderr.read_text(),code)
            record['memory_peak'] = peak
            verify_record(row,output,error,images[name])
            if parity:
                record['parity_sha256'] = sha(parity)
                compare(parity,references[int(row['args'][1])])
                with parity.open('rb') as stream:
                    os.fsync(stream.fileno()); os.posix_fadvise(stream.fileno(),0,0,os.POSIX_FADV_DONTNEED)
            equal(sha(binary),build['artifacts'][name])
            print(label,code,peak,flush=True)
        verify_build(root)
        native_reference(root)
        for name,digest in state['tools'].items():
            equal(sha(source/name),digest); equal(sha(root/'check_tools'/name),digest)
        state['completed'] = True
    finally:
        (out/'checks.json').write_text(json.dumps(state,indent=2)+'\n')


if __name__ == '__main__': check(Path(sys.argv[1]).resolve(strict=True))

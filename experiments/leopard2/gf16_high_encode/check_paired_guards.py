#!/usr/bin/env python3
"""Build/check negative controls for the paired frontend's actual guards."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

from verify_paired_r19932 import BEAD, equal, sha


def main(action, root):
    folder = root/'guards'
    source = Path(__file__).resolve().with_name('paired_guard_test.cpp')
    if action == 'build':
        folder.mkdir()
        shutil.copyfile(source,folder/source.name)
        state = dict(bead=BEAD,completed=False,timed=False,commands=[],artifacts={})
        try:
            original = json.loads((root/'build/build.json').read_text())
            for profile in ('release','sanitize'):
                # Reuse exactly the qualified frontend flags, headers, object and archive pins.
                command = next(c for c in original['commands'] if '-c' in c and
                               c[-1] == str(root/'build'/profile/'driver.o'))
                flags = command[command.index('c++')+1:command.index('-c')]
                args = ['prlimit','--cpu=120:120','--','c++',*flags,'-I'+str(root/'build'),
                        str(folder/source.name),str(root/'build'/profile/'codec.a'),
                        str(root/'build/clock_guard.cpp'),'-Wl,--wrap=_ZNSt6chrono3_V212steady_clock3nowEv',
                        '-o',str(folder/profile)]
                state['commands'].append(args)
                subprocess.run(args,check=True)
            for path in sorted(folder.iterdir()):
                state['artifacts'][path.name] = sha(path)
                path.chmod(0o555 if path.name in ('release','sanitize') else 0o444)
            equal(sha(source),sha(folder/source.name))
            state['completed'] = True
        finally:
            (folder/'build.json').write_text(json.dumps(state,indent=2)+'\n')
    elif action == 'check':
        records = []
        state = dict(bead=BEAD,completed=False,timed=False,records=records)
        env = dict(os.environ,ASAN_OPTIONS='detect_leaks=1:quarantine_size_mb=8:thread_local_quarantine_size_kb=64',
                   UBSAN_OPTIONS='halt_on_error=1:print_stacktrace=1')
        try:
            for profile,case,code in (('release','canary',0),('sanitize','canary',0),
                                      ('sanitize','underflow',1),('sanitize','overflow',1)):
                label = profile+'-'+case
                args = ['prlimit','--cpu=10:10','--',str(folder/profile),case]
                with (folder/(label+'.stdout')).open('x') as out, (folder/(label+'.stderr')).open('x') as err:
                    result = subprocess.run(args,env=env,stdout=out,stderr=err,check=False)
                records.append(dict(label=label,args=args,returncode=result.returncode,
                    stdout_sha256=sha(folder/(label+'.stdout')),stderr_sha256=sha(folder/(label+'.stderr'))))
                equal(result.returncode,code)
                stderr = (folder/(label+'.stderr')).read_text()
                if code:
                    if 'ERROR: AddressSanitizer: use-after-poison' not in stderr:
                        raise ValueError('poison boundary did not trigger ASan')
                else:
                    equal(stderr,'')
                    equal((folder/(label+'.stdout')).read_text(),
                          'both canaries rejected; zero-size and restored buffers pass\n')
                print(label,'expected',code,flush=True)
            state['completed'] = True
        finally:
            (folder/'checks.json').write_text(json.dumps(state,indent=2)+'\n')
    else:
        raise ValueError('action')


if __name__=='__main__':
    main(sys.argv[1],Path(sys.argv[2]).resolve())

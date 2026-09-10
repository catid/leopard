#!/usr/bin/env python3
"""Build new metadata frontends only, serial under 512MiB/no-swap canonical lock."""
import json
from pathlib import Path
import shutil
import subprocess
import sys

from paired_metadata_overlay import BEAD, adapt
from verify_paired_r19932 import ARCHIVES, require, sha

SOURCE = Path(__file__).resolve().parent
REPO = SOURCE.parents[2]
QUALIFIED = REPO / '.research/leopard-79h/auto-r19932-qualified.515m59j2/build'


def build(root):
    out = root / 'build'
    out.mkdir()
    state = dict(bead=BEAD, completed=False, timed=False, inputs={}, commands=[], artifacts={})

    def pin(path):
        digest = sha(path)
        state['inputs'][str(path)] = digest
        return digest

    def run(args):
        command = ['prlimit', '--cpu=120:120', '--nofile=65536:1048576', '--', *map(str,args)]
        state['commands'].append(command)
        print(json.dumps(command), flush=True)
        subprocess.run(command, check=True)

    try:
        for name in ('build_paired_metadata.py','paired_metadata_overlay.py','PairedRuntimeMetadata.h',
                     'paired_timer_r19932.cpp','paired_timer_witness.cpp','paired_public_witness.cpp',
                     'paired_timer_clock.cpp','PairedGroupTiming.h','tower_public_scope.sh'):
            pin(SOURCE/name); shutil.copyfile(SOURCE/name,out/name)
        (out/'paired_metadata.cpp').write_text(adapt((out/'paired_timer_r19932.cpp').read_text()))
        guard = QUALIFIED/'drivers/clock_guard.cpp'
        pin(guard); shutil.copyfile(guard,out/'clock_guard.cpp')
        for profile,digest in ARCHIVES.items():
            native,sanitize = profile == 'native', profile == 'sanitize'
            directory = out/profile; directory.mkdir()
            archive = QUALIFIED/profile/'candidate.a'
            require(pin(archive) == digest,'qualified archive identity')
            shutil.copyfile(archive,directory/'codec.a')
            include = (REPO/'.research/leopard-79h/avx2-isa-attribution.qqNHjs/pure-checks-v2/source'
                       if native else QUALIFIED/'source')
            headers = directory/'include'; headers.mkdir()
            for header in sorted(include.glob('*.h')):
                pin(header); shutil.copyfile(header,headers/header.name)
            flags = ['-std=gnu++11','-Wall','-Wextra','-Werror','-fopenmp','-pthread',
                '-march=x86-64','-mtune=generic','-mavx2','-mno-avx512f','-I'+str(headers),
                '--param=ggc-min-expand=10','--param=ggc-min-heapsize=4096']
            flags += (['-O1','-g1','-fsanitize=address,undefined','-fno-sanitize-recover=all',
                       '-fno-omit-frame-pointer','-fno-pie','-no-pie'] if sanitize else ['-O3','-DNDEBUG'])
            macros = ['-DLEO_PAIRED_CODEC="'+profile+':'+digest+'"']
            if native: macros += ['-DLEO_PAIRED_NATIVE=1']
            for src,obj in (('paired_metadata.cpp','driver.o'),('paired_timer_witness.cpp','witness.o')):
                run(['c++',*flags,*macros,'-c',out/src,'-o',directory/obj])
            wrappers = ['-Wl,--wrap=leo_encode'] if native else ['-Wl,--wrap=leo2_encode','-Wl,--wrap=leo2_encode_batch']
            for kind in ('abort','synthetic'):
                clock = directory/(kind+'-clock.o')
                run(['c++',*flags,'-DLEO_PAIRED_'+kind.upper()+'=1','-c',out/'paired_timer_clock.cpp','-o',clock])
                extras = [out/'clock_guard.cpp'] if kind == 'abort' else []
                run(['c++',*flags,directory/'driver.o',directory/'witness.o',clock,*extras,*wrappers,
                    '-Wl,--wrap=_ZNSt6chrono3_V212steady_clock3nowEv',directory/'codec.a','-ldl','-o',directory/kind])
                undefined = subprocess.check_output(['nm','-u',directory/kind],text=True)
                require('_ZNSt6chrono3_V212steady_clock3nowEv' not in undefined,'real benchmark clock linkage')
                (directory/(kind+'-undefined.txt')).write_text(undefined)
                symbols = subprocess.check_output(['nm','-n','--defined-only',directory/kind],text=True)
                (directory/(kind+'-symbols.txt')).write_text(symbols)
        for name,digest in state['inputs'].items(): require(sha(Path(name)) == digest,'build input drift')
        for path in sorted(out.rglob('*')):
            if path.is_file():
                state['artifacts'][str(path.relative_to(out))] = sha(path)
                path.chmod(0o555 if path.name in ('abort','synthetic') else 0o444)
        state['completed'] = True
    finally:
        (out/'build.json').write_text(json.dumps(state,indent=2)+'\n')


if __name__ == '__main__': build(Path(sys.argv[1]).resolve(strict=True))

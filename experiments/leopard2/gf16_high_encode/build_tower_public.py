#!/usr/bin/env python3
"""Link tower public frontends, not codecs. Caller owns serial 512M/no-swap lock."""
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

from tower_public_overlay import adapt

BEAD = 'leopard-79h.38.5.4.18.4.2'
SOURCE = Path(__file__).resolve().parent
REPO = SOURCE.parents[2]
TOWER = REPO/'.research/leopard-79h/tower-encoder-qualified.kvqrsu/build'
QUALIFIED = REPO/'.research/leopard-79h/auto-r19932-qualified.515m59j2/build'
ARCHIVES = {
    'native': '3f83c55599cbfbeea5d15adf3627405bcde1cbed2ddf4a672fbc75b186f5c4c1',
    'original': '89f33d3d43f792de078c38469f9406dded2b69126af2c6fc8db6e2acb77cd334',
    'original-sanitize': 'c2aee8119a36bc647a450649106656437c55ca83fe08337a559edea6c533b0a9',
    'release': '24a78109405a25a5164e656bcaabbcc9aff586849a7899492cbaf290b3018ba9',
    'trace': 'fb07bdae206cfdd909c294e1c8d5910db8614349ab93103c60f819fcd4e25b87',
    'sanitize': 'a580b3b95460ff6a5104c88cd5bd8a40b9ba6058e8a3e4c3043e69d28da1b1c4',
}
SOURCES = ('tower_public_overlay.py','tower_public_link.h','tower_public_link.cpp',
           'tower_encoder.h','avx2_adjacent_public.cpp','avx2_adjacent_public_witness.cpp',
           'paired_timer_clock.cpp','PairedGroupTiming.h','test_paired_group_timing.cpp',
           'test_avx2_adjacent_public_guards.cpp','tower_public_scope.sh')


def require(ok, message):
    if not ok: raise ValueError(message)


def sha(path):
    with path.open('rb') as stream:
        digest = hashlib.file_digest(stream, 'sha256').hexdigest()
        os.posix_fadvise(stream.fileno(),0,0,os.POSIX_FADV_DONTNEED)
        return digest


def build(root):
    out = root/'build'; out.mkdir()
    state = dict(bead=BEAD,root=str(root),completed=False,timed=False,inputs={},commands=[],artifacts={})
    def pin(path):
        state['inputs'][str(path)] = sha(path)
        return state['inputs'][str(path)]
    def run(argv):
        command = ['prlimit','--cpu=120:120','--nofile=65536:1048576','--',*map(str,argv)]
        index = len(state['commands']); state['commands'].append(command)
        with (out/f'command-{index}.stdout').open('xb') as a, (out/f'command-{index}.stderr').open('xb') as b:
            subprocess.run(command,stdout=a,stderr=b,check=True,timeout=180)
    try:
        pin(Path(__file__).resolve())
        for name in SOURCES:
            pin(SOURCE/name); shutil.copyfile(SOURCE/name,out/name)
        for base,kind,target in (('avx2_adjacent_public.cpp','driver','tower_public.cpp'),
                                 ('avx2_adjacent_public_witness.cpp','witness','tower_public_witness.cpp')):
            (out/target).write_text(adapt((out/base).read_text(),kind))
        # This tests the actual generated driver's guards, without running a codec.
        (out/'test_tower_public_guards.cpp').write_text(
            (out/'test_avx2_adjacent_public_guards.cpp').read_text().replace(
                '#include "avx2_adjacent_public.cpp"','#include "tower_public.cpp"'))
        guard = QUALIFIED/'drivers/clock_guard.cpp'
        pin(guard); shutil.copyfile(guard,out/'clock_guard.cpp')
        compiled = {}
        for profile,digest in ARCHIVES.items():
            native, san = profile == 'native', 'sanitize' in profile
            original = profile.startswith('original')
            directory = out/profile; directory.mkdir()
            archive = (QUALIFIED/'native/candidate.a' if native else
                       TOWER/('original-sanitize.a' if san else 'original-release.a') if original else
                       TOWER/profile/'candidate.a')
            require(pin(archive) == digest,'archive identity: '+profile)
            shutil.copyfile(archive,directory/'codec.a')
            include = (REPO/'.research/leopard-79h/avx2-isa-attribution.qqNHjs/pure-checks-v2/source'
                       if native else QUALIFIED/'source')
            headers = out/('native-headers' if native else 'l2-headers')
            if not headers.exists():
                headers.mkdir()
                for header in sorted(include.glob('*.h')):
                    pin(header); shutil.copyfile(header,headers/header.name)
            flags = ['-std=gnu++11','-Wall','-Wextra','-Werror','-fopenmp','-pthread',
                     '-march=x86-64','-mtune=generic','-mavx2','-mno-avx512f','-I'+str(headers),
                     '--param=ggc-min-expand=10','--param=ggc-min-heapsize=4096']
            flags += (['-O1','-g1','-fsanitize=address,undefined','-fno-sanitize-recover=all',
                       '-fno-omit-frame-pointer','-fno-pie','-no-pie'] if san else ['-O3','-DNDEBUG'])
            macros = ['-DLEO_PAIRED_NATIVE=1'] if native else []
            common = 'native' if native else 'sanitize' if san else 'l2-release'
            if common not in compiled:
                shared = out/(common+'-objects'); shared.mkdir(); compiled[common] = shared
                for source,obj in (('tower_public.cpp','driver.o'),('tower_public_witness.cpp','witness.o')):
                    run(['c++',*flags,*macros,'-c',out/source,'-o',shared/obj])
                for kind in ('steady','synthetic','abort'):
                    define = [] if kind == 'steady' else ['-DLEO_PAIRED_'+kind.upper()+'=1']
                    run(['c++',*flags,*define,'-c',out/'paired_timer_clock.cpp','-o',shared/(kind+'-clock.o')])
            shared = compiled[common]
            identity = ['-DLEO_TOWER_IDENTITY="'+profile+':'+digest+'"']
            if original: identity += ['-DLEO_TOWER_ORIGINAL=1']
            run(['c++',*flags,*macros,*identity,'-c',out/'tower_public_link.cpp','-o',directory/'link.o'])
            wrappers = (['-Wl,--wrap=leo_encode'] if native else
                        ['-Wl,--wrap=leo2_encode','-Wl,--wrap=leo2_encode_batch'])
            for variant in ('plain','abort','synthetic'):
                extras = [] if variant == 'plain' else [shared/'witness.o',*wrappers,
                          '-Wl,--wrap=_ZNSt6chrono3_V212steady_clock3nowEv']
                extras += [shared/('steady-clock.o' if variant == 'plain' else variant+'-clock.o')]
                if variant == 'abort': extras += [out/'clock_guard.cpp']
                run(['c++',*flags,shared/'driver.o',directory/'link.o',*extras,directory/'codec.a',
                     '-o',directory/variant])
                symbols = subprocess.check_output(['nm','-u',str(directory/variant)],text=True)
                require(('_ZNSt6chrono3_V212steady_clock3nowEv' in symbols) == (variant == 'plain'),
                        'clock symbol binding')
                (directory/(variant+'-undefined.txt')).write_text(symbols)
            if profile in ('release','sanitize'):
                run(['c++',*flags,out/'test_paired_group_timing.cpp','-o',directory/'group-unit'])
                run(['c++',*flags,out/'test_tower_public_guards.cpp',directory/'link.o',
                     shared/'abort-clock.o',directory/'codec.a',out/'clock_guard.cpp',
                     '-Wl,--wrap=_ZNSt6chrono3_V212steady_clock3nowEv','-o',directory/'guards'])
            print(profile+' public frontend linked',flush=True)
        for name,digest in state['inputs'].items(): require(sha(Path(name)) == digest,'input drift')
        for path in sorted(out.rglob('*')):
            if path.is_file():
                state['artifacts'][str(path.relative_to(out))] = sha(path)
                path.chmod(0o555 if path.name in ('plain','abort','synthetic','group-unit','guards') else 0o444)
        state['completed'] = True
    finally:
        (out/'build.json').write_text(json.dumps(state,indent=2)+'\n')


if __name__ == '__main__':
    require(len(sys.argv) == 2,'usage: build_tower_public.py NEW_WORKSPACE')
    build(Path(sys.argv[1]).resolve(strict=True))

#!/usr/bin/env python3
"""Build public frontends only. Caller holds serial 512 MiB/no-swap lock."""
import json
from pathlib import Path
import shutil
import subprocess
import sys

from verify_paired_r19932 import equal, require, sha

BEAD = 'leopard-79h.38.5.4.18.3.1'
SOURCE = Path(__file__).resolve().parent
REPO = SOURCE.parents[2]
QUALIFIED = REPO/'.research/leopard-79h/auto-r19932-qualified.515m59j2/build'
RUNTIME = REPO/'.research/leopard-79h/avx2-adjacent-runtime-focused.yMDquC/build'
ARCHIVES = {
    'native': '3f83c55599cbfbeea5d15adf3627405bcde1cbed2ddf4a672fbc75b186f5c4c1',
    'original': '89f33d3d43f792de078c38469f9406dded2b69126af2c6fc8db6e2acb77cd334',
    'release': '542c340e53c1c607feb44e4d8b048bdcef41a7a04c2829c9fb62d27acd624644',
    'trace': 'e15e2c86f07b4713157c0a71ab3f3d2c5babe951cc17de9f3f08a98f033f9094',
    'sanitize': '501d5b0375455b137c8c24befb2e152a0ea2d773f8b094aa1388d0b5e6cb4410',
}
SOURCES = ('avx2_adjacent_public.cpp', 'avx2_adjacent_public_link.h',
           'avx2_adjacent_public_link.cpp', 'avx2_adjacent_public_witness.cpp',
           'avx2_adjacent_control.h', 'paired_timer_clock.cpp', 'PairedGroupTiming.h',
           'test_paired_group_timing.cpp')


def build(root):
    out = root/'build'; out.mkdir()
    state = dict(bead=BEAD, completed=False, timed=False, inputs={}, commands=[], artifacts={})

    def pin(path):
        state['inputs'][str(path)] = sha(path)
        return state['inputs'][str(path)]

    def run(argv):
        command = ['prlimit', '--cpu=120:120', '--nofile=65536:1048576', '--', *map(str, argv)]
        index = len(state['commands']); state['commands'].append(command)
        with (out/f'command-{index}.stdout').open('xb') as stdout, (out/f'command-{index}.stderr').open('xb') as stderr:
            subprocess.run(command, stdout=stdout, stderr=stderr, check=True, timeout=180)

    try:
        pin(Path(__file__).resolve())
        for name in SOURCES:
            pin(SOURCE/name); shutil.copyfile(SOURCE/name, out/name)
        guard = QUALIFIED/'drivers/clock_guard.cpp'
        pin(guard); shutil.copyfile(guard, out/'clock_guard.cpp')
        compiled = {}
        for profile, digest in ARCHIVES.items():
            native, san = profile == 'native', profile == 'sanitize'
            directory = out/profile; directory.mkdir()
            archive = (QUALIFIED/('native' if native else 'release') if profile in ('native','original')
                       else RUNTIME/profile)/'candidate.a'
            equal(pin(archive), digest)
            shutil.copyfile(archive, directory/'codec.a')
            include = (REPO/'.research/leopard-79h/avx2-isa-attribution.qqNHjs/pure-checks-v2/source'
                       if native else QUALIFIED/'source')
            headers = out/('native-headers' if native else 'l2-headers')
            if not headers.exists():
                headers.mkdir()
                for header in sorted(include.glob('*.h')):
                    pin(header); shutil.copyfile(header, headers/header.name)
            flags = ['-std=gnu++11', '-Wall', '-Wextra', '-Werror', '-fopenmp', '-pthread',
                     '-march=x86-64', '-mtune=generic', '-mavx2', '-mno-avx512f', '-I'+str(headers),
                     '--param=ggc-min-expand=10', '--param=ggc-min-heapsize=4096']
            flags += (['-O1', '-g1', '-fsanitize=address,undefined', '-fno-sanitize-recover=all',
                       '-fno-omit-frame-pointer', '-fno-pie', '-no-pie'] if san else ['-O3', '-DNDEBUG'])
            macros = ['-DLEO_PAIRED_NATIVE=1'] if native else []
            common = 'native' if native else 'sanitize' if san else 'l2-release'
            if common not in compiled:
                shared = out/(common+'-objects'); shared.mkdir(); compiled[common] = shared
                for source, obj in (('avx2_adjacent_public.cpp','driver.o'),
                                    ('avx2_adjacent_public_witness.cpp','witness.o')):
                    run(['c++', *flags, *macros, '-c', out/source, '-o', shared/obj])
                for kind in ('steady','synthetic','abort'):
                    define = {'steady':[], 'synthetic':['-DLEO_PAIRED_SYNTHETIC=1'],
                              'abort':['-DLEO_PAIRED_ABORT=1']}[kind]
                    run(['c++', *flags, *define, '-c', out/'paired_timer_clock.cpp',
                         '-o', shared/(kind+'-clock.o')])
            shared = compiled[common]
            identity = ['-DLEO_ADJACENT_IDENTITY="'+profile+':'+digest+'"']
            if profile == 'original': identity += ['-DLEO_ADJACENT_ORIGINAL=1']
            run(['c++', *flags, *macros, *identity, '-c', out/'avx2_adjacent_public_link.cpp',
                 '-o', directory/'link.o'])
            wrappers = (['-Wl,--wrap=leo_encode'] if native else
                        ['-Wl,--wrap=leo2_encode','-Wl,--wrap=leo2_encode_batch'])
            for variant in ('plain','abort','synthetic'):
                extras = [] if variant == 'plain' else [shared/'witness.o',*wrappers,
                          '-Wl,--wrap=_ZNSt6chrono3_V212steady_clock3nowEv']
                extras += [shared/('steady-clock.o' if variant == 'plain' else variant+'-clock.o')]
                if variant == 'abort': extras += [out/'clock_guard.cpp']
                run(['c++', *flags, shared/'driver.o', directory/'link.o', *extras,
                     directory/'codec.a', '-o', directory/variant])
                symbols = subprocess.check_output(['nm','-u',str(directory/variant)], text=True)
                equal('_ZNSt6chrono3_V212steady_clock3nowEv' in symbols, variant == 'plain')
                (directory/(variant+'-undefined.txt')).write_text(symbols)
            if profile in ('release','sanitize'):
                run(['c++', *flags, out/'test_paired_group_timing.cpp', '-o', directory/'group-unit'])
            print(profile+' public frontend linked', flush=True)
        for name, digest in state['inputs'].items(): equal(sha(Path(name)), digest)
        for path in sorted(out.rglob('*')):
            if path.is_file():
                state['artifacts'][str(path.relative_to(out))] = sha(path)
                path.chmod(0o555 if path.name in ('plain','abort','synthetic','group-unit') else 0o444)
        state['completed'] = True
    finally:
        (out/'build.json').write_text(json.dumps(state, indent=2)+'\n')


if __name__ == '__main__':
    require(len(sys.argv) == 2, 'usage: build_avx2_adjacent_public.py NEW_WORKSPACE')
    build(Path(sys.argv[1]).resolve(strict=True))

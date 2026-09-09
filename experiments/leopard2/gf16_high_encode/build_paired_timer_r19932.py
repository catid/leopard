#!/usr/bin/env python3
"""Build timer adapters only; no execution. 512 MiB/no-swap/serial lock required."""
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

from verify_auto_r19932_checks import ARCHIVES

BEAD = 'leopard-79h.38.5.4.19.1.1'
QUALIFIED = Path('/tmp/leopard-auto-r19932.9EGR0p/build')
REPO = Path(__file__).resolve().parents[3]
SOURCE = Path(__file__).resolve().parent


def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def build(root):
    out = root / 'build'
    out.mkdir()  # Never replace an earlier build, including a failed one.
    state = dict(bead=BEAD, completed=False, timed=False, inputs={}, commands=[], artifacts={})

    def pin(path):
        digest = sha(path)
        state['inputs'][str(path)] = digest
        return digest

    def run(argv):
        command = ['prlimit', '--cpu=120:120', '--nofile=65536:1048576', '--', *map(str, argv)]
        state['commands'].append(command)
        print(json.dumps(command), flush=True)
        subprocess.run(command, check=True)

    try:
        pin(Path(__file__).resolve())
        for name in ('paired_timer_r19932.cpp', 'paired_public_witness.cpp', 'paired_timer_witness.cpp',
                     'paired_timer_clock.cpp', 'PairedGroupTiming.h', 'test_paired_group_timing.cpp'):
            pin(SOURCE / name)
            shutil.copyfile(SOURCE / name, out / name)
        guard = QUALIFIED / 'drivers/clock_guard.cpp'
        pin(guard)
        shutil.copyfile(guard, out / 'clock_guard.cpp')
        for profile, digest in ARCHIVES.items():
            native, san = profile == 'native', profile == 'sanitize'
            directory = out / profile
            directory.mkdir()
            archive = QUALIFIED / profile / 'candidate.a'
            if pin(archive) != digest:
                raise ValueError('qualified archive changed: ' + profile)
            shutil.copyfile(archive, directory / 'codec.a')
            include = (REPO / '.research/leopard-79h/avx2-isa-attribution.qqNHjs/pure-checks-v2/source'
                       if native else QUALIFIED / 'source')
            for header in sorted(include.glob('*.h')):
                pin(header)
            flags = ['-std=gnu++11', '-Wall', '-Wextra', '-Werror', '-fopenmp', '-pthread',
                     '-march=x86-64', '-mtune=generic', '-mavx2', '-mno-avx512f', '-I'+str(include),
                     '--param=ggc-min-expand=10', '--param=ggc-min-heapsize=4096']
            flags += (['-O1', '-g1', '-fsanitize=address,undefined', '-fno-sanitize-recover=all',
                       '-fno-omit-frame-pointer', '-fno-pie', '-no-pie'] if san else ['-O3', '-DNDEBUG'])
            macros = ['-DLEO_PAIRED_CODEC="'+profile+':'+digest+'"']
            if native:
                macros += ['-DLEO_PAIRED_NATIVE=1']
            for source, obj in (('paired_timer_r19932.cpp', 'driver.o'), ('paired_timer_witness.cpp', 'witness.o')):
                run(['c++', *flags, *macros, '-c', out/source, '-o', directory/obj])
            for kind in ('steady','synthetic','abort'):
                define = {'steady':[], 'synthetic':['-DLEO_PAIRED_SYNTHETIC=1'], 'abort':['-DLEO_PAIRED_ABORT=1']}[kind]
                run(['c++', *flags, *define,
                     '-c',out/'paired_timer_clock.cpp','-o',
                     directory/(kind+'-clock.o')])
            wrappers = (['-Wl,--wrap=leo_encode'] if native else
                        ['-Wl,--wrap=leo2_encode','-Wl,--wrap=leo2_encode_batch'])
            for variant in ('plain','abort','synthetic'):
                extras = [] if variant=='plain' else [directory/'witness.o',*wrappers,
                          '-Wl,--wrap=_ZNSt6chrono3_V212steady_clock3nowEv']
                extras += [directory/('steady-clock.o' if variant=='plain' else variant+'-clock.o')]
                if variant=='abort':
                    extras += [out/'clock_guard.cpp']
                run(['c++',*flags,directory/'driver.o',*extras,directory/'codec.a','-o',directory/variant])
                symbols = subprocess.check_output(['nm','-u',str(directory/variant)],text=True)
                if ('_ZNSt6chrono3_V212steady_clock3nowEv' in symbols) != (variant=='plain'):
                    raise ValueError('unexpected clock linkage: '+variant)
                (directory/(variant+'-undefined.txt')).write_text(symbols)
            if not native:
                run(['c++',*flags,out/'test_paired_group_timing.cpp','-o',directory/'group-unit'])
        for name, digest in state['inputs'].items():
            if sha(Path(name)) != digest:
                raise ValueError('build input changed: '+name)
        for path in sorted(out.rglob('*')):
            if path.is_file():
                state['artifacts'][str(path.relative_to(out))] = sha(path)
                path.chmod(0o555 if path.name in ('plain', 'abort', 'synthetic', 'group-unit') else 0o444)
        state['completed'] = True
    finally:
        (out/'build.json').write_text(json.dumps(state, indent=2)+'\n')


if __name__ == '__main__':
    build(Path(sys.argv[1]).resolve())

#!/usr/bin/env python3
"""Bounded serial build recipe; no codec execution, leopard-79h.38.5.4.18.4."""
import hashlib
import json
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
from tower_encoder_overlay import overlay, BASE

REPO = Path('/home/catid/leopard')
EXP = REPO / 'experiments/leopard2/gf16_high_encode'
PRIOR = REPO / '.research/leopard-79h/auto-r19932-qualified.515m59j2/build'
ARCHIVES = {'release': '89f33d3d43f792de078c38469f9406dded2b69126af2c6fc8db6e2acb77cd334',
            'sanitize': 'c2aee8119a36bc647a450649106656437c55ca83fe08337a559edea6c533b0a9'}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def build(root):
    source = root / 'source'
    source.mkdir()
    originals = {}
    for path in list(REPO.glob('*.h')) + [REPO/'LeopardFF16.cpp'] + [EXP/name for name in (
            'tower_encoder.h', 'tower_encoder.cpp', 'tower_encoder_overlay.py',
            'tower_butterfly_probe.h', 'tower_butterfly_probe.cpp', 'tower_avx2_probe.h',
            'test_tower_encoder.cpp', 'test_gfni_boundary.cpp', 'gfni_boundary_screen.cpp',
            'test_tower_encoder_kernels.cpp', 'build_tower_encoder.py')]:
        originals[str(path)] = sha(path)
        shutil.copy2(path, source/path.name)
    require(sha(source/'LeopardFF16.cpp') == BASE, 'field base')
    (source/'LeopardFF16.cpp').rename(source/'LeopardFF16.original.cpp')
    (source/'LeopardFF16.cpp').write_text(overlay((source/'LeopardFF16.original.cpp').read_text()))
    for path in source.iterdir():
        path.chmod(0o444)
    records = dict(tracker='leopard-79h.38.5.4.18.4', base_commit='9763aee', timed=False,
                   completed=False, original_sources=originals,
                   sources={p.name: sha(p) for p in source.iterdir()}, commands=[], profiles={})
    def run(args, cwd):
        prefix = root / ('command-' + str(len(records['commands'])))
        command = ['prlimit', '--cpu=120:120', '--', *args]
        records['commands'].append(dict(argv=command, cwd=str(cwd)))
        with prefix.with_suffix('.stdout').open('xb') as out, prefix.with_suffix('.stderr').open('xb') as err:
            subprocess.run(command, cwd=cwd, stdout=out, stderr=err, check=True, timeout=180)
    try:
        for profile in ('release', 'trace', 'sanitize'):
            dest = root/profile
            dest.mkdir()
            san = profile == 'sanitize'
            original_profile = 'sanitize' if san else 'release'
            original = PRIOR/original_profile/'candidate.a'
            require(sha(original) == ARCHIVES[original_profile], 'archive identity')
            recipe_path = REPO / ('.research/leopard-79h/gf16-current-route-failed.STc10h/build-metadata/sanitize-gc-build/compile_commands.json'
                                  if san else '.research/leopard-79h/auto-gfni-boundary-production.B5NV9R/release/compile_commands.json')
            recipes = json.loads(recipe_path.read_text())
            recipe, = [x for x in recipes if x['output'] == 'CMakeFiles/leopard.dir/LeopardFF16.cpp.o']
            field = dest/'LeopardFF16.cpp.o'
            if profile == 'trace':
                shutil.copy2(root/'release/LeopardFF16.cpp.o', field)
            else:
                args = shlex.split(recipe['command'])
                args[args.index('-I'+str(REPO))] = '-I'+str(source)
                args[args.index('-c')+1] = str(source/'LeopardFF16.cpp')
                args[args.index('-o')+1] = str(field)
                run(args, dest)
            flags = ['-std=c++11', '-Wall', '-Wextra', '-Werror', '-fopenmp', '-pthread',
                     '-march=x86-64', '-mtune=generic', '-I'+str(source),
                     '--param=ggc-min-expand=10', '--param=ggc-min-heapsize=4096']
            flags += (['-O1', '-g1', '-fsanitize=address,undefined', '-fno-sanitize-recover=all',
                       '-fno-omit-frame-pointer', '-fno-pie', '-no-pie'] if san else ['-O3', '-DNDEBUG'])
            control = dest/'tower_encoder.cpp.o'
            probe = dest/'tower_butterfly_probe.cpp.o'
            trace = profile != 'release'
            run(['c++', *flags, '-DLEO_TOWER_TRACE='+str(int(trace)), '-c', str(source/'tower_encoder.cpp'),
                 '-o', str(control)], dest)
            if profile == 'trace':
                shutil.copy2(root/'release/tower_butterfly_probe.cpp.o', probe)
            else:
                run(['c++', *flags, '-mavx2', '-mno-avx512f', '-mno-gfni', '-c',
                     str(source/'tower_butterfly_probe.cpp'), '-o', str(probe)], dest)
            archive = dest/'candidate.a'
            shutil.copy2(original, archive)
            archive.chmod(0o600)
            run(['ar', 'r', str(archive), str(field), str(control), str(probe)], dest)
            run(['ar', 's', str(archive)], dest)
            before = subprocess.check_output(['ar','t',str(original)], text=True).splitlines()
            after = subprocess.check_output(['ar','t',str(archive)], text=True).splitlines()
            require(len(before) == len(set(before)) == 24, 'original member inventory')
            require(after == before + [control.name, probe.name], 'candidate member inventory')
            unchanged = {}
            for name in before:
                if name == field.name:
                    continue
                old = subprocess.check_output(['ar','p',str(original),name])
                require(old == subprocess.check_output(['ar','p',str(archive),name]), 'unrelated member drift')
                unchanged[name] = hashlib.sha256(old).hexdigest()
            run(['c++', *flags, '-DLEO_BOUNDARY_CODEC_COMMIT="tower-experiment"',
                 str(source/'test_tower_encoder.cpp'), str(archive), '-o', str(dest/'test')], dest)
            run(['c++', *flags, str(source/'test_tower_encoder_kernels.cpp'), str(archive),
                 '-o', str(dest/'test-kernels')], dest)
            if profile != 'trace':
                run(['c++', *flags, '-Wno-unused-function', '-DLEO_TOWER_ORIGINAL_CONTROL=1',
                     '-DLEO_BOUNDARY_CODEC_COMMIT="original-45e2eff"', str(source/'test_tower_encoder.cpp'),
                     str(original), '-o', str(dest/'test-original')], dest)
            for obj in (field, control, probe):
                with (dest/(obj.name+'.disassembly')).open('xb') as out:
                    subprocess.run(['objdump','-drwC','-Mintel',str(obj)], stdout=out, check=True, timeout=10)
            records['profiles'][profile] = dict(original_archive=str(original), original_sha256=sha(original),
                recipe=recipe, recipe_file=str(recipe_path), recipe_sha256=sha(recipe_path),
                unchanged_members=unchanged, files={p.name:sha(p) for p in dest.iterdir()})
            for path in dest.iterdir():
                path.chmod(0o555 if path.name in ('test','test-kernels','test-original') else 0o444)
            print(profile+' '+sha(archive), flush=True)
        for name, digest in originals.items():
            require(sha(Path(name)) == digest, 'original source drift')
        records['completed'] = True
    finally:
        with (root/'build.json').open('x') as out:
            json.dump(records, out, indent=2)


if __name__ == '__main__':
    require(len(sys.argv) == 2, 'usage: build_tower_encoder.py NEW_ROOT')
    build(Path(sys.argv[1]).resolve(strict=True))

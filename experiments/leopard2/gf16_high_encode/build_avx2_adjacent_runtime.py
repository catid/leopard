#!/usr/bin/env python3
"""Serial, 512 MiB capped experiment build; no codec or benchmark execution."""
import hashlib
import json
from pathlib import Path
import shlex
import shutil
import subprocess
import sys

from avx2_adjacent_runtime import overlay
from verify_paired_r19932 import equal,require,sha
from verify_auto_r19932_checks import ARCHIVES

BEAD = 'leopard-79h.38.5.4.18.3.1'
REPO = Path('/home/catid/leopard')
EXPERIMENT = REPO/'experiments/leopard2/gf16_high_encode'
QUALIFIED = REPO/'.research/leopard-79h/auto-r19932-qualified.515m59j2'
STATIC = REPO/'.research/leopard-79h/avx2-adjacent-qualified.wxkx9v2b'
MEMBER = 'Leopard2BackendAVX2.cpp.o'
DRIVERS = ('avx2_adjacent_control.cpp','avx2_adjacent_control.h','avx2_adjacent_runtime.py',
    'avx2_adjacent_schedule.h','test_avx2_adjacent_control.cpp','test_avx2_adjacent_runtime.cpp',
    'test_avx2_adjacent_schedule.cpp','test_avx2_pair_schedule.cpp','test_gfni_boundary.cpp',
    'gfni_boundary_screen.cpp','build_avx2_adjacent_runtime.py')


def build(root):
    root = root.resolve(strict=True)
    out = root/'build'; out.mkdir()
    source,drivers = out/'source',out/'drivers'
    shutil.copytree(QUALIFIED/'build/source',source)
    source.chmod(0o700)  # Only the new owned copy; retained input stays read-only.
    drivers.mkdir()
    for name in DRIVERS: shutil.copyfile(EXPERIMENT/name,drivers/name)
    focused = (drivers/'test_avx2_adjacent_schedule.cpp').read_text()
    require(focused.count('int main(int argc, char** argv)')==1,'focused entry anchor')
    (drivers/'adjacent_focused.inc').write_text(focused.replace(
        'int main(int argc, char** argv)','int LeoAdjacentFocusedMain(int argc, char** argv)',1))
    # The arithmetic overlay is the exact earlier-qualified mode3 source.
    base = STATIC/'source/Leopard2BackendAVX2.cpp'
    (source/'Leopard2BackendAVX2.cpp').chmod(0o600)
    (source/'Leopard2BackendAVX2.cpp').write_text(overlay(base.read_text()))
    for name in ('avx2_adjacent_schedule.h','avx2_adjacent_control.h'):
        if (source/name).exists():
            equal(sha(source/name),sha(EXPERIMENT/name))
        else:
            shutil.copyfile(EXPERIMENT/name,source/name)
    pins = {str(p):sha(p) for d in (source,drivers) for p in d.iterdir()}
    state = dict(bead=BEAD,timed=False,completed=False,source_base='45e2eff',inputs=pins,
                 static_source_sha256=sha(base),commands=[],profiles={},external_inputs={str(base):sha(base)})
    def run(argv,directory):
        argv = ['/usr/bin/prlimit','--cpu=120:120','--nofile=65536:1048576','--',*argv]
        index = len(state['commands']); state['commands'].append(argv)
        with (out/f'command-{index}.stdout').open('xb') as stdout, (out/f'command-{index}.stderr').open('xb') as stderr:
            subprocess.run(argv,cwd=directory,stdout=stdout,stderr=stderr,check=True,timeout=180)
    try:
        for profile in ('release','trace','sanitize'):
            directory = out/profile; directory.mkdir(); san = profile=='sanitize'
            original_profile = 'sanitize' if san else 'release'
            original = QUALIFIED/'build'/original_profile/'candidate.a'
            equal(sha(original),ARCHIVES[original_profile])
            state['external_inputs'][str(original)] = ARCHIVES[original_profile]
            recipes_path = (REPO/'.research/leopard-79h/gf16-current-route-failed.STc10h/build-metadata/sanitize-gc-build/compile_commands.json'
                            if san else REPO/'.research/leopard-79h/auto-gfni-boundary-production.B5NV9R/release/compile_commands.json')
            recipes = json.loads(recipes_path.read_text())
            state['external_inputs'][str(recipes_path)] = sha(recipes_path)
            recipe, = [r for r in recipes if r['output']=='CMakeFiles/leopard2_backend_avx2.dir/'+MEMBER]
            obj = directory/MEMBER
            argv = shlex.split(recipe['command'])
            argv[argv.index('-I'+str(REPO))] = '-I'+str(source)
            argv[argv.index('-c')+1] = str(source/'Leopard2BackendAVX2.cpp')
            argv[argv.index('-o')+1] = str(obj)
            trace = int(profile!='release')
            argv += ['-DLEO_AVX2_ADJACENT_SCHEDULE=3','-DLEO_ADJACENT_TRACE='+str(trace)]
            run(argv,directory)
            flags = ['-std=gnu++11','-Wall','-Wextra','-Werror','-fopenmp','-pthread',
                '-march=x86-64','-mtune=generic','-mavx2','-mno-avx512f','-I'+str(source),
                '--param=ggc-min-expand=10','--param=ggc-min-heapsize=4096','-DLEO_ADJACENT_TRACE='+str(trace)]
            flags += (['-O1','-g1','-fsanitize=address,undefined','-fno-sanitize-recover=all',
                       '-fno-omit-frame-pointer','-fno-pie','-no-pie'] if san else ['-O3','-DNDEBUG'])
            control = directory/'avx2_adjacent_control.cpp.o'
            run(['c++',*flags,'-c',str(drivers/'avx2_adjacent_control.cpp'),'-o',str(control)],directory)
            archive = directory/'candidate.a'; shutil.copyfile(original,archive)
            run(['ar','r',str(archive),str(obj),str(control)],directory)
            run(['ar','s',str(archive)],directory)
            before = subprocess.check_output(['ar','t',str(original)],text=True).splitlines()
            after = subprocess.check_output(['ar','t',str(archive)],text=True).splitlines()
            require(len(before)==24 and len(set(before))==24,'original inventory')
            equal(after,before+[control.name])
            unchanged = {}
            for name in before:
                if name==MEMBER: continue
                old = subprocess.check_output(['ar','p',str(original),name])
                require(old==subprocess.check_output(['ar','p',str(archive),name]),'unrelated member drift')
                unchanged[name] = hashlib.sha256(old).hexdigest()
            run(['c++',*flags,str(drivers/'test_avx2_adjacent_control.cpp'),str(archive),
                 '-o',str(directory/'control-unit')],directory)
            run(['c++',*flags,'-DLEO_BOUNDARY_CODEC_COMMIT="adjacent-runtime:'+sha(archive)+'"',
                 str(drivers/'test_avx2_adjacent_runtime.cpp'),str(archive),'-o',str(directory/'focused')],directory)
            with (directory/'disassembly.txt').open('xb') as stream:
                subprocess.run(['objdump','-drwC',str(obj)],stdout=stream,check=True)
            state['profiles'][profile] = dict(original=str(original),original_sha256=sha(original),
                archive_sha256=sha(archive),object_sha256=sha(obj),recipe=recipe,recipe_sha256=sha(recipes_path),
                unchanged_members=unchanged,files={n:sha(directory/n) for n in
                    ('control-unit','focused','avx2_adjacent_control.cpp.o','disassembly.txt')})
            print(profile+' runtime archive '+sha(archive),flush=True)
        for path,digest in pins.items(): equal(sha(Path(path)),digest)
        for path,digest in state['external_inputs'].items(): equal(sha(Path(path)),digest)
        state['completed'] = True
    finally:
        with (out/'build.json').open('x') as stream: json.dump(state,stream,indent=2); stream.write('\n')


if __name__=='__main__':
    require(len(sys.argv)==2,'usage: build_avx2_adjacent_runtime.py NEW_WORKSPACE')
    build(Path(sys.argv[1]))

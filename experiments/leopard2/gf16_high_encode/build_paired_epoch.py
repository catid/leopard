#!/usr/bin/env python3
"""Reuse pinned metadata compile recipes for a separate clock-free frontend."""
import ast
import json
from pathlib import Path
import shutil
import subprocess
import sys

from paired_epoch_overlay import BEAD, driver, header, clock, witness
from verify_paired_metadata import sha, parse, equal, require, elf

SOURCE = Path(__file__).resolve().parent
REFERENCE = SOURCE.parents[2] / '.research/leopard-79h/paired-runtime-metadata.oQioRv/build'
REFERENCE_PIN = 'f8f2529034e462ccbcfcc94ec3951182a3497a34e1258e01b32b49e3354ae543'
BASE_NAMES = ('paired_timer_r19932.cpp','PairedRuntimeMetadata.h','paired_timer_clock.cpp',
              'paired_timer_witness.cpp','paired_public_witness.cpp','PairedGroupTiming.h',
              'clock_guard.cpp','tower_public_scope.sh')

def copy_tools(source, output, names):
    output.mkdir()
    pins = {}
    def copy(name):
        if name in pins: return
        path = source/name
        pins[name] = sha(path); shutil.copyfile(path, output/name)
        (output/name).chmod(0o444)
        for node in ast.walk(ast.parse(path.read_text())):
            modules = ([node.module] if isinstance(node,ast.ImportFrom) else
                       [a.name for a in node.names] if isinstance(node,ast.Import) else [])
            for module in modules:
                if module and (source/(module+'.py')).is_file(): copy(module+'.py')
    for name in names: copy(name)
    for name, digest in pins.items(): equal(sha(source/name), digest)
    return pins

def recipes(reference, output):
    old = str(Path(reference['commands'][0][-1]).parents[1])
    require(old.endswith('/build'), 'reference build root')
    return [[arg.replace(old+'/', str(output)+'/') for arg in command]
            for command in reference['commands']]

def build(root):
    out = root/'build'; out.mkdir(); (out/'baseline').mkdir()
    equal(sha(REFERENCE/'build.json'), REFERENCE_PIN)
    old = parse((REFERENCE/'build.json').read_text())
    shutil.copyfile(REFERENCE/'build.json', out/'reference_build.json')
    state = dict(bead=BEAD, root=str(root), completed=False, timed=False,
                 reference_build_sha256=REFERENCE_PIN, baseline={}, commands=[], artifacts={})
    state['tools'] = copy_tools(SOURCE, out/'source_tools', ['build_paired_epoch.py'])
    try:
        for name in BASE_NAMES:
            equal(sha(REFERENCE/name), old['artifacts'][name])
            shutil.copyfile(REFERENCE/name, out/'baseline'/name)
            shutil.copyfile(REFERENCE/name, out/name)
            state['baseline'][name] = old['artifacts'][name]
        (out/'paired_metadata.cpp').write_text(driver((out/'baseline/paired_timer_r19932.cpp').read_text()))
        for name, adapt in (('PairedRuntimeMetadata.h',header),('paired_timer_clock.cpp',clock),
                            ('paired_timer_witness.cpp',witness)):
            (out/name).write_text(adapt((out/'baseline'/name).read_text()))
        for profile in ('native','release','sanitize'):
            folder = out/profile; folder.mkdir(); (folder/'include').mkdir()
            for path in [REFERENCE/profile/'codec.a', *sorted((REFERENCE/profile/'include').glob('*.h'))]:
                relative = str(path.relative_to(REFERENCE))
                equal(sha(path), old['artifacts'][relative])
                shutil.copyfile(path, out/relative)
                state['baseline'][relative] = old['artifacts'][relative]
        for command in recipes(old,out):
            state['commands'].append(command); print(json.dumps(command),flush=True)
            subprocess.run(command,check=True)
        for profile in ('native','release','sanitize'):
            for kind in ('abort','synthetic'):
                elf(out/profile/kind,profile=='native')  # also forbids a real clock import
        for path in sorted(out.rglob('*')):
            if path.is_file():
                state['artifacts'][str(path.relative_to(out))] = sha(path)
                path.chmod(0o555 if path.name in ('abort','synthetic') else 0o444)
        for name,digest in state['tools'].items(): equal(sha(SOURCE/name),digest)
        state['completed'] = True
    finally:
        (out/'build.json').write_text(json.dumps(state,indent=2)+'\n')

if __name__ == '__main__': build(Path(sys.argv[1]).resolve(strict=True))

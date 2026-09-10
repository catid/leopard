#!/usr/bin/env python3
"""Serial frontend-only build; never executes a codec or benchmark clock."""
import json
from pathlib import Path
import shutil
import subprocess
import sys

from build_paired_epoch import copy_tools
from paired_epoch_timing_overlay import BEAD, driver, fake_steady_clock, PLAIN_MARK
from paired_epoch_timing_elf import elf, VARIANTS, CLOCK
from verify_paired_metadata import sha, parse, equal, require
from verify_paired_r19932 import ARCHIVES

SOURCE = Path(__file__).resolve().parent
REFERENCE = SOURCE.parents[2]/'.research/leopard-79h/paired-epoch-qualified.unzse9/build'
REFERENCE_PIN = '7de095968e41289605e9164d29f585de934b822f570cca72d1bb4c3c44d5ca72'
ASSETS = ('paired_metadata.cpp','PairedRuntimeMetadata.h','PairedGroupTiming.h',
          'paired_timer_clock.cpp','paired_timer_witness.cpp','paired_public_witness.cpp',
          'clock_guard.cpp','tower_public_scope.sh')
REUSED = ('codec.a','witness.o','abort-clock.o','synthetic-clock.o')


def executable(profile,variant):
    require(profile in ARCHIVES and variant in VARIANTS, 'executable selection')
    return ('measurement' if variant == 'steady' else 'qualification')+'/'+profile+'/'+variant


def commands(reference,out):
    old = str(Path(reference['root'])/'build')
    result = []
    for profile in ARCHIVES:
        folder = out/profile
        original = [c for c in reference['commands'] if c[-1] == old+'/'+profile+'/driver.o']
        require(len(original) == 1, 'one qualified driver recipe')
        base = [a.replace(old+'/',str(out)+'/') for a in original[0][:-4]]
        for source,obj,extra in (
            ('paired_epoch_diagnostic.cpp','driver.o',()),
            ('plain_marks.cpp','plain-marks.o',()),
            ('paired_timer_clock.cpp','steady-clock.o',()),
            ('fake_steady_clock.cpp','fake-steady-clock.o',('-DLEO_PAIRED_SYNTHETIC=1',))):
            result.append([*base,*extra,'-c',str(out/source),'-o',str(folder/obj)])
        for variant in VARIANTS:
            witnessed = variant in ('synthetic','fake-steady')
            wrapper = (['-Wl,--wrap=leo_encode'] if profile == 'native' else
                       ['-Wl,--wrap=leo2_encode','-Wl,--wrap=leo2_encode_batch']) if witnessed else []
            # Metadata captures __real_* addresses in the common driver.o.
            # Plain links alias those names to the actual public functions;
            # there is no public-call wrapper or forwarding function.
            aliases = [] if witnessed else ['-Wl,--defsym=__real_'+name+'='+name for name in
                (('leo_encode',) if profile == 'native' else ('leo2_encode','leo2_encode_batch'))]
            result.append([*base,str(folder/'driver.o'),str(folder/('witness.o' if witnessed else 'plain-marks.o')),
                str(folder/(variant+'-clock.o')), *([str(out/'clock_guard.cpp')] if variant == 'abort' else []),
                *wrapper, *aliases, *(['-Wl,--wrap='+CLOCK] if variant != 'steady' else []),
                str(folder/'codec.a'),'-ldl','-o',str(out/executable(profile,variant))])
    return result


def build(root):
    out = root/'build'; out.mkdir(); baseline = out/'baseline'; baseline.mkdir()
    equal(sha(REFERENCE/'build.json'),REFERENCE_PIN)
    reference = parse((REFERENCE/'build.json').read_text())
    equal([reference['completed'],reference['timed']],[True,False])
    shutil.copyfile(REFERENCE/'build.json',out/'reference_build.json')
    state = dict(bead=BEAD,root=str(root),completed=False,real_clocks_read=False,
                 reference_build_sha256=REFERENCE_PIN,inputs={},commands=[],images={},artifacts={})
    state['tools'] = copy_tools(SOURCE,out/'source_tools',[Path(__file__).name])
    try:
        def copy(name,target):
            path = REFERENCE/name
            equal(sha(path),reference['artifacts'][name])
            require(not path.is_symlink() and not path.stat().st_mode&0o222, 'immutable reference')
            target.parent.mkdir(parents=True,exist_ok=True)
            shutil.copyfile(path,target); state['inputs'][name] = sha(path)
        for name in ASSETS:
            copy(name,baseline/name)
            if name != 'paired_metadata.cpp': shutil.copyfile(baseline/name,out/name)
        (out/'paired_epoch_diagnostic.cpp').write_text(driver((baseline/'paired_metadata.cpp').read_text()))
        (out/'fake_steady_clock.cpp').write_text(fake_steady_clock((baseline/'paired_timer_clock.cpp').read_text()))
        (out/'plain_marks.cpp').write_text(PLAIN_MARK)
        for profile,digest in ARCHIVES.items():
            for name in REUSED: copy(profile+'/'+name,out/profile/name)
            equal(sha(out/profile/'codec.a'),digest)
            for name in reference['artifacts']:
                if name.startswith(profile+'/include/'): copy(name,out/name)
            for variant in VARIANTS: (out/executable(profile,variant)).parent.mkdir(parents=True,exist_ok=True)
        for command in commands(reference,out):
            state['commands'].append(command); print(json.dumps(command),flush=True)
            subprocess.run(command,check=True)
        for profile in ARCHIVES:
            for variant in VARIANTS:
                name = executable(profile,variant)
                state['images'][name] = elf(out/name,profile=='native',variant)
        for name,digest in state['inputs'].items(): equal(sha(REFERENCE/name),digest)
        for name,digest in state['tools'].items(): equal(sha(SOURCE/name),digest)
        executables = set(state['images'])
        for path in sorted(out.rglob('*')):
            if path.is_file():
                name = str(path.relative_to(out)); state['artifacts'][name] = sha(path)
                path.chmod(0o555 if name in executables else 0o444)
        state['completed'] = True
    finally:
        (out/'build.json').write_text(json.dumps(state,indent=2)+'\n')


if __name__ == '__main__': build(Path(sys.argv[1]).resolve(strict=True))

#!/usr/bin/env python3
"""Freeze final replay/test sources only after native collection is terminal.

The captured build and collection closures remain unchanged. This does not
freeze a timing campaign, read clocks, or enable its collection interface.
"""
import json
from pathlib import Path
import shutil
import sys

from build_paired_epoch import copy_tools
import verify_paired_epoch_timing as verify


def freeze(root):
    source=Path(__file__).resolve().parent
    output=root/'final_tools'; manifest=root/'final-tools.json'
    verify.require(not output.exists() and not output.is_symlink() and
                   not manifest.exists() and not manifest.is_symlink(),'fresh final tool inventory')
    state=verify.parse((root/'checks/checks.json').read_text())
    verify.equal([state['bead'],state['completed'],state['real_clocks_read']],[verify.BEAD,True,False])
    reference=root/'build/reference_build.json'
    verify.equal(verify.sha(reference),verify.REFERENCE_PIN)
    baseline=verify.parse(reference.read_text())['baseline']
    # Validate all source assets before creating a partially populated freeze.
    for name in verify.FINAL_ASSETS: verify.equal(verify.sha(source/name),baseline[name])
    pins=copy_tools(source,output,verify.FINAL_ROOTS)
    verify.equal(sorted(pins),sorted(verify.FINAL_PYTHON))
    for name in verify.FINAL_ASSETS:
        shutil.copyfile(source/name,output/name); (output/name).chmod(0o444)
        pins[name]=verify.sha(output/name); verify.equal(pins[name],baseline[name])
    with manifest.open('x') as stream:
        json.dump(dict(bead=verify.BEAD,real_clocks_read=False,files=pins),stream,sort_keys=True,indent=2)
        stream.write('\n')
    manifest.chmod(0o444); output.chmod(0o555)
    verify.verify_final_tools(root)


if __name__=='__main__':
    verify.require(len(sys.argv)==2,'usage: freeze_paired_epoch_timing.py TERMINAL_ROOT')
    freeze(Path(sys.argv[1]).resolve(strict=True))

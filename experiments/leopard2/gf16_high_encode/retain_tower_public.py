#!/usr/bin/env python3
"""Freeze this lane's evidence; deduplicate only parity copies inside the bundle."""
import json
import os
from pathlib import Path
import shutil
import sys

from build_tower_public import require, sha


def copy(source, target):
    target.parent.mkdir(parents=True,exist_ok=True)
    with source.open('rb') as incoming, target.open('xb') as outgoing:
        shutil.copyfileobj(incoming,outgoing,65536)
        outgoing.flush(); os.fsync(outgoing.fileno())
        for stream in (incoming,outgoing):
            os.posix_fadvise(stream.fileno(),0,0,os.POSIX_FADV_DONTNEED)
    require(sha(source)==sha(target),'copy identity')


def retain(raw, output):
    require(output.is_dir() and not any(output.iterdir()),'empty owned retention directory required')
    state = json.loads((raw/'checks/checks.json').read_text())
    require(state['completed'] is True and state['timed'] is False,'completed untimed checks')
    seen = {}; shared = 0
    for source in sorted(raw.rglob('*')):
        require(not source.is_symlink(),'no source links')
        if not source.is_file(): continue
        target = output/source.relative_to(raw)
        target.parent.mkdir(parents=True,exist_ok=True)
        digest = sha(source)
        if source.suffix == '.parity' and digest in seen:
            # The first copy is privately owned by this bundle. Never link to
            # raw evidence or other lanes, and never deduplicate executables.
            os.link(seen[digest],target); shared += source.stat().st_size
        else:
            copy(source,target)
            if source.suffix == '.parity': seen[digest] = target
        require(sha(target)==digest,'retained identity')
    # Final verifier is separate from the original collector/import snapshot.
    for source in sorted(Path(__file__).parent.glob('*.py')):
        copy(source,output/'replay-tools'/source.name)
    prior = Path('/home/catid/leopard/.research/leopard-79h/avx2-adjacent-public.NFkTeG')
    copy(prior/'SHA256SUMS',output/'prior-native/SHA256SUMS')
    for c in range(9):
        name = f'checks/native-{c}-NNNN-1-plain.parity'
        copy(prior/name,output/'prior-native'/name)
    lines,logical,allocated = [],0,0
    inodes = set()
    for path in sorted(output.rglob('*')):
        if not path.is_file(): continue
        lines.append(sha(path)+'  '+str(path.relative_to(output)))
        info = path.stat(); logical += info.st_size
        if info.st_ino not in inodes: allocated += info.st_blocks*512; inodes.add(info.st_ino)
        path.chmod(0o444)
    manifest = output/'SHA256SUMS'
    manifest.write_text('\n'.join(lines)+'\n'); manifest.chmod(0o444)
    for path in sorted(output.rglob('*'),reverse=True):
        if path.is_dir(): path.chmod(0o555)
    output.chmod(0o555)
    print(json.dumps(dict(files=len(lines)+1,logical_bytes=logical+manifest.stat().st_size,
        allocated_bytes=allocated+manifest.stat().st_blocks*512,parity_duplicate_bytes_saved=shared,
        manifest_sha256=sha(manifest)),sort_keys=True))


if __name__ == '__main__':
    require(len(sys.argv)==3,'usage: retain_tower_public.py RAW EMPTY_OUTPUT')
    retain(Path(sys.argv[1]).resolve(strict=True),Path(sys.argv[2]).resolve(strict=True))

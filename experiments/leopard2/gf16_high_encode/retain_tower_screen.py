#!/usr/bin/env python3
"""Copy a terminal screen into a new readonly bundle; never runs native code.

Caller owns the canonical lock and a 256 MiB/no-swap scope after the timing
job has exited. Both complete and interrupted attempts remain retainable.
Tracker: leopard-79h.38.5.4.18.4.4.
"""
import hashlib
import json
import os
from pathlib import Path
import shutil
import sys


def require(ok, message):
    if not ok: raise ValueError(message)


def sha(path):
    with path.open('rb') as stream:
        value = hashlib.file_digest(stream,'sha256').hexdigest()
        os.posix_fadvise(stream.fileno(),0,0,os.POSIX_FADV_DONTNEED)
        return value


def retain(source, output):
    require(source.is_dir() and output.is_dir(),'existing directories required')
    require(not source.is_symlink() and not output.is_symlink(),'no root links')
    require(source.resolve()!=output.resolve() and source.resolve() not in output.resolve().parents
            and output.resolve() not in source.resolve().parents,'disjoint trees')
    require(not any(output.iterdir()),'new empty owned destination required')
    require(output.stat().st_uid==os.getuid() and output.stat().st_mode & 0o777==0o700,
            'private destination ownership')
    require((source/'attempt1/attempt.json').is_file(),'terminal attempt journal required')
    require(not (source/'SHA256SUMS').exists() and not (source/'SHA256SUMS').is_symlink(),
            'reserved root manifest name')
    # Actual job termination is checked by the caller; a file alone is not a
    # live-process check. Preserve interrupted attempts without inventing ratios.
    files = []
    for path in sorted(source.rglob('*')):
        require(not path.is_symlink(),'no source links')
        if path.is_dir(): continue
        require(path.is_file(),'regular files only')
        require(path.stat().st_size<64*1024**2,'bounded member size')
        files.append(path)
    require(0<len(files)<4096,'bounded file inventory')
    require(sum(p.stat().st_size for p in files)<256*1024**2,'bounded total input')
    manifest=[]; total=0
    for path in files:
        relative=path.relative_to(source)
        target=output/relative; target.parent.mkdir(parents=True,exist_ok=True)
        before=sha(path)
        with path.open('rb') as incoming,target.open('xb') as outgoing:
            shutil.copyfileobj(incoming,outgoing,65536)
            outgoing.flush(); os.fsync(outgoing.fileno())
            for stream in (incoming,outgoing):
                os.posix_fadvise(stream.fileno(),0,0,os.POSIX_FADV_DONTNEED)
        require(sha(path)==before==sha(target),'copy identity')
        require(path.stat().st_ino!=target.stat().st_ino or path.stat().st_dev!=target.stat().st_dev,
                'private copy, not hard link')
        target.chmod(0o444); total+=target.stat().st_size
        manifest.append(before+'  '+str(relative))
    require([p for p in sorted(source.rglob('*')) if p.is_file()]==files,'source inventory drift')
    for path,line in zip(files,manifest):
        require(not path.is_symlink() and sha(path)==line.split('  ',1)[0],
                'source changed after copy')
    checksum=output/'SHA256SUMS'
    with checksum.open('x') as stream: stream.write('\n'.join(manifest)+'\n')
    checksum.chmod(0o444)
    for path in sorted(output.rglob('*'),reverse=True):
        if path.is_dir(): path.chmod(0o555)
    output.chmod(0o555)
    print(json.dumps(dict(files=len(files)+1,bytes=total+checksum.stat().st_size,
                         manifest_sha256=sha(checksum)),sort_keys=True))


if __name__=='__main__':
    require(len(sys.argv)==3,'usage: retain_tower_screen.py TERMINAL_RAW EMPTY_DESTINATION')
    retain(Path(sys.argv[1]),Path(sys.argv[2]))

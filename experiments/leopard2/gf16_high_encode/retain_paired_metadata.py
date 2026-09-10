#!/usr/bin/env python3
"""Seal untimed metadata evidence only; no codec execution. Local capped lock required."""
import ast
import json
import os
from pathlib import Path
import shutil
import sys

from paired_metadata_overlay import BEAD
from verify_paired_metadata import parse
from verify_paired_r19932 import equal, require
from retain_tower_screen import sha


def freeze(root):
    source=Path(__file__).resolve().parent
    output=root/'final_tools'; output.mkdir()
    pins={}
    def copy(name):
        if name in pins: return
        path=source/name; pins[name]=sha(path)
        shutil.copyfile(path,output/name); (output/name).chmod(0o444)
        if path.suffix!='.py': return
        for node in ast.walk(ast.parse(path.read_text())):
            names=([node.module] if isinstance(node,ast.ImportFrom) else
                   [a.name for a in node.names] if isinstance(node,ast.Import) else [])
            for name in names:
                if name and (source/(name+'.py')).is_file(): copy(name+'.py')
    for name in ('verify_paired_metadata.py','test_paired_metadata.py','test_paired_metadata_replay.py','build_paired_metadata_units.py',
                 'retain_paired_metadata.py','paired_timer_r19932.cpp','test_paired_runtime_metadata.cpp'):
        copy(name)
    for name in ('paired_metadata_overlay.py','PairedRuntimeMetadata.h','build_paired_metadata.py'):
        equal(sha(source/name),sha(root/'build'/name))
    for name,digest in pins.items(): equal(sha(source/name),digest)
    (root/'final-tools.json').write_text(json.dumps(dict(bead=BEAD,timed=False,files=pins),sort_keys=True,indent=2)+'\n')


def retain(source,output):
    require(source.is_dir() and output.is_dir() and not source.is_symlink() and not output.is_symlink(),'real directories')
    require(source.resolve()!=output.resolve() and source.resolve() not in output.resolve().parents
            and output.resolve() not in source.resolve().parents,'disjoint trees')
    require(not any(output.iterdir()) and output.stat().st_uid==os.getuid() and output.stat().st_mode&0o777==0o700,
            'new private owned destination')
    for path in ('checks/checks.json','units/checks.json'):
        state=parse((source/path).read_text()); equal([state['completed'],state['timed']],[True,False])
    equal((source/'replay-normal.json').read_text(),(source/'replay-optimized.json').read_text())
    require(not (source/'SHA256SUMS').exists() and not (source/'SHA256SUMS').is_symlink(),'reserved manifest name')
    files=[]
    for path in sorted(source.rglob('*')):
        require(not path.is_symlink(),'no links')
        if path.is_dir(): continue
        require(path.is_file() and path.stat().st_size<64*1024**2,'bounded regular file')
        files.append(path)
    require(0<len(files)<4096 and sum(p.stat().st_size for p in files)<4*1024**3,'bounded evidence inventory')
    manifest=[]; total=0
    for path in files:
        relative=path.relative_to(source); target=output/relative; target.parent.mkdir(parents=True,exist_ok=True)
        before=sha(path)
        with path.open('rb') as a,target.open('xb') as b:
            shutil.copyfileobj(a,b,65536); b.flush(); os.fsync(b.fileno())
            for stream in (a,b): os.posix_fadvise(stream.fileno(),0,0,os.POSIX_FADV_DONTNEED)
        equal(sha(path),before); equal(sha(target),before)
        require((path.stat().st_dev,path.stat().st_ino)!=(target.stat().st_dev,target.stat().st_ino),'private copy')
        target.chmod(0o444); total+=target.stat().st_size
        manifest.append(before+'  '+str(relative))
    require([p for p in sorted(source.rglob('*')) if p.is_file()]==files,'source inventory drift')
    for path,line in zip(files,manifest): equal(sha(path),line.split('  ',1)[0])
    checksum=output/'SHA256SUMS'; checksum.write_text('\n'.join(manifest)+'\n'); checksum.chmod(0o444)
    for path in sorted(output.rglob('*'),reverse=True):
        if path.is_dir(): path.chmod(0o555)
    output.chmod(0o555)
    print(json.dumps(dict(bead=BEAD,files=len(files)+1,bytes=total+checksum.stat().st_size,
                         manifest_sha256=sha(checksum)),sort_keys=True))


if __name__=='__main__':
    if len(sys.argv)==3 and sys.argv[1]=='freeze': freeze(Path(sys.argv[2]).resolve(strict=True))
    elif len(sys.argv)==4 and sys.argv[1]=='retain': retain(Path(sys.argv[2]),Path(sys.argv[3]))
    else: raise SystemExit('usage: freeze ROOT | retain ROOT NEW_PRIVATE_DESTINATION')

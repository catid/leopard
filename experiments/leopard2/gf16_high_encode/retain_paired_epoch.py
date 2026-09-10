#!/usr/bin/env python3
"""Private readonly epoch evidence copies, including explicitly stopped history.

No codec execution. Run under the campaign lock and 256MiB/no-swap cap after
collection is terminal. Never modify raw evidence or an occupied destination.
"""
import json
import os
from pathlib import Path
import shutil
import sys
from build_paired_epoch import copy_tools
from verify_paired_epoch import BEAD, BASE_NAMES, parse, equal, require, sha


def freeze(root):
    source=Path(__file__).resolve().parent
    output=root/'final_tools'
    manifest=root/'final-tools.json'
    require(not output.exists() and not output.is_symlink() and not manifest.exists() and not manifest.is_symlink(),
            'fresh final tool inventory')
    pins=copy_tools(source,output,['verify_paired_epoch.py','verify_paired_epoch_units.py','test_paired_epoch.py',
                                 'retain_paired_epoch.py','test_retain_paired_epoch.py','build_paired_epoch_units.py'])
    build=parse((root/'build/build.json').read_text())
    for name in BASE_NAMES:
        path=root/'build/baseline'/name
        equal(sha(path),build['baseline'][name])
        require(not (output/name).exists(),'unique frozen source')
        shutil.copyfile(path,output/name); (output/name).chmod(0o444); pins[name]=sha(path)
    for name in ('test_paired_epoch_native.cpp','test_paired_epoch_clock.cpp'):
        path=root/'units'/name
        equal(sha(path),sha(source/name))
        shutil.copyfile(path,output/name); (output/name).chmod(0o444); pins[name]=sha(path)
    for name in ('paired_epoch_overlay.py','build_paired_epoch.py'):
        equal(sha(source/name),sha(root/'build/source_tools'/name))
    for name,digest in pins.items(): equal(sha(output/name),digest)
    with manifest.open('x') as stream:
        stream.write(json.dumps(dict(bead=BEAD,timed=False,files=pins),sort_keys=True,indent=2)+'\n')


def retain(source,output,stopped=None):
    inputs=[('',source)] + ([('stopped_attempt',stopped)] if stopped is not None else [])
    if stopped is not None:
        require(source.resolve()!=stopped.resolve() and source.resolve() not in stopped.resolve().parents and
                stopped.resolve() not in source.resolve().parents,'disjoint complete and stopped histories')
    require(output.is_dir() and not output.is_symlink(),'real destination')
    require(not any(output.iterdir()) and output.stat().st_uid==os.getuid() and output.stat().st_mode&0o777==0o700,
            'new private owned destination')
    for prefix,folder in inputs:
        require(folder.is_dir() and not folder.is_symlink(),'real source')
        require(folder.resolve()!=output.resolve() and folder.resolve() not in output.resolve().parents and
                output.resolve() not in folder.resolve().parents,'disjoint source and destination')
        if prefix:
            state=parse((folder/'checks/checks.json').read_text())
            equal([state['bead'],state['completed'],state['timed']],[BEAD,False,False])
            require(state['records'] and state['records'][-1]['returncode']==143,'explicitly stopped native child')
        else:
            for name in ('checks/checks.json','units/checks.json'):
                state=parse((folder/name).read_text())
                equal([state['bead'],state['completed'],state['timed']],[BEAD,True,False])
            equal((folder/'replay-normal.json').read_text(),(folder/'replay-optimized.json').read_text())
            require(not (folder/'SHA256SUMS').exists() and not (folder/'SHA256SUMS').is_symlink(),'reserved manifest')
            require(not (folder/'stopped_attempt').exists(),'reserved history prefix')
    files=[]
    def scan():
        found=[]
        for prefix,folder in inputs:
            for path in sorted(folder.rglob('*')):
                require(not path.is_symlink(),'no source links')
                if path.is_dir(): continue
                require(path.is_file() and path.stat().st_size<64*1024**2,'bounded regular source')
                relative=str(Path(prefix)/path.relative_to(folder))
                require('\n' not in relative and '\r' not in relative,'manifest-safe name')
                found.append((relative,path))
        require(0<len(found)<4096 and sum(p.stat().st_size for _,p in found)<4*1024**3,'bounded evidence tree')
        require(len({r for r,_ in found})==len(found),'unique evidence paths')
        return found
    files=scan(); manifest=[]; total=0
    for relative,path in files:
        target=output/relative; target.parent.mkdir(parents=True,exist_ok=True)
        before=sha(path)
        with path.open('rb') as a,target.open('xb') as b:
            shutil.copyfileobj(a,b,65536); b.flush(); os.fsync(b.fileno())
            for stream in (a,b): os.posix_fadvise(stream.fileno(),0,0,os.POSIX_FADV_DONTNEED)
        equal(sha(path),before); equal(sha(target),before)
        require((path.stat().st_dev,path.stat().st_ino)!=(target.stat().st_dev,target.stat().st_ino),'private copy')
        target.chmod(0o444); total+=target.stat().st_size; manifest.append(before+'  '+relative)
    equal([(r,str(p)) for r,p in scan()],[(r,str(p)) for r,p in files])
    for (_,path),line in zip(files,manifest): equal(sha(path),line.split('  ',1)[0])
    checksum=output/'SHA256SUMS'; checksum.write_text('\n'.join(manifest)+'\n'); checksum.chmod(0o444)
    for path in sorted(output.rglob('*'),reverse=True):
        if path.is_dir(): path.chmod(0o555)
    output.chmod(0o555)
    return dict(bead=BEAD,files=len(files)+1,bytes=total+checksum.stat().st_size,
                manifest_sha256=sha(checksum),stopped_attempt_included=stopped is not None)


if __name__=='__main__':
    if len(sys.argv)==3 and sys.argv[1]=='freeze': freeze(Path(sys.argv[2]).resolve(strict=True))
    elif len(sys.argv) in (4,5) and sys.argv[1]=='retain':
        print(json.dumps(retain(Path(sys.argv[2]),Path(sys.argv[3]),Path(sys.argv[4]) if len(sys.argv)==5 else None),sort_keys=True))
    else: raise SystemExit('usage: freeze ROOT | retain ROOT NEW_PRIVATE_DESTINATION [STOPPED_ROOT]')

#!/usr/bin/env python3
"""Private readonly frontend evidence, preserving the initial failed build.

Caller confirms all jobs terminal and owns the canonical capped lock. This
copies evidence; it neither runs a codec nor qualifies a timing collector.
"""
import json
import os
from pathlib import Path
import shutil
import sys
from verify_paired_epoch_timing import BEAD,sha,parse,require,equal


def retain(source,output,failed):
    inputs=[('',source),('failed_build',failed)]
    require(output.is_dir() and not output.is_symlink() and not any(output.iterdir()),'fresh destination')
    require(output.stat().st_uid==os.getuid() and output.stat().st_mode&0o777==0o700,'private destination')
    for a,b in ((source,failed),(source,output),(failed,output)):
        require(a.resolve()!=b.resolve() and a.resolve() not in b.resolve().parents and
                b.resolve() not in a.resolve().parents,'disjoint evidence trees')
    for _,folder in inputs:
        require(folder.is_dir() and not folder.is_symlink(),'real source directory')
    checks=parse((source/'checks/checks.json').read_text())
    equal([checks['bead'],checks['completed'],checks['real_clocks_read']],[BEAD,True,False])
    old=parse((failed/'build/build.json').read_text())
    equal([old['bead'],old['completed'],old['real_clocks_read']],[BEAD,False,False])
    require(not (failed/'checks').exists(),'failed build precedes native execution')
    equal((source/'replay-normal.json').read_text(),(source/'replay-optimized.json').read_text())
    for name in ('SHA256SUMS','failed_build'):
        require(not (source/name).exists() and not (source/name).is_symlink(),'reserved evidence name')
    def scan():
        files=[]
        for prefix,folder in inputs:
            for path in sorted(folder.rglob('*')):
                require(not path.is_symlink(),'no source links')
                if path.is_dir(): continue
                require(path.is_file() and path.stat().st_size<64*1024**2,'bounded regular source')
                name=str(Path(prefix)/path.relative_to(folder))
                require('\n' not in name and '\r' not in name,'manifest-safe path')
                files.append((name,path))
        require(0<len(files)<4096 and sum(p.stat().st_size for _,p in files)<4*1024**3,'bounded evidence')
        equal(len({n for n,_ in files}),len(files))
        return files
    files=scan(); manifest=[]; total=0
    for name,path in files:
        target=output/name; target.parent.mkdir(parents=True,exist_ok=True)
        digest=sha(path)
        with path.open('rb') as a,target.open('xb') as b:
            shutil.copyfileobj(a,b,65536); b.flush(); os.fsync(b.fileno())
            for stream in (a,b): os.posix_fadvise(stream.fileno(),0,0,os.POSIX_FADV_DONTNEED)
        equal(sha(path),digest); equal(sha(target),digest)
        require((path.stat().st_dev,path.stat().st_ino)!=(target.stat().st_dev,target.stat().st_ino),'private copy')
        target.chmod(0o444); total+=target.stat().st_size; manifest.append(digest+'  '+name)
    equal([(n,str(p)) for n,p in scan()],[(n,str(p)) for n,p in files])
    for (_,path),line in zip(files,manifest): equal(sha(path),line.split('  ',1)[0])
    checksum=output/'SHA256SUMS'
    with checksum.open('x') as stream: stream.write('\n'.join(manifest)+'\n')
    checksum.chmod(0o444)
    for path in sorted(output.rglob('*'),reverse=True):
        if path.is_dir(): path.chmod(0o555)
    output.chmod(0o555)
    return dict(bead=BEAD,files=len(files)+1,bytes=total+checksum.stat().st_size,
                manifest_sha256=sha(checksum),failed_build_included=True,collector_qualified=False)


if __name__=='__main__':
    require(len(sys.argv)==4,'usage: retain_paired_epoch_timing.py ROOT NEW_PRIVATE_DESTINATION FAILED_BUILD')
    print(json.dumps(retain(*(Path(p) for p in sys.argv[1:])),sort_keys=True))

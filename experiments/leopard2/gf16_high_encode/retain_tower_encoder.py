#!/usr/bin/env python3
"""Retain owned qualification evidence without deleting earlier failures."""
import hashlib
import json
import os
from pathlib import Path
import shutil
import sys


def copy(source, target):
    target.parent.mkdir(parents=True,exist_ok=True)
    with source.open('rb') as incoming, target.open('xb') as outgoing:
        shutil.copyfileobj(incoming,outgoing,65536)
        outgoing.flush(); os.fsync(outgoing.fileno())
        os.posix_fadvise(outgoing.fileno(),0,0,os.POSIX_FADV_DONTNEED)


def retain(raw, native, output):
    if output.exists():
        if not output.is_dir() or any(output.iterdir()): raise ValueError('retention output must be empty')
    else: output.mkdir()
    for source in sorted(raw.rglob('*')):
        if source.is_file(): copy(source,output/'build'/source.relative_to(raw))
    metadata = json.loads((raw/'build.json').read_text())
    for profile in ('release','sanitize'):
        copy(Path(metadata['profiles'][profile]['original_archive']), output/'build'/('original-'+profile+'.a'))
    # Preserve native oracle source metadata and the exact manifest entries
    # already checked by replay; do not copy or run any remote artifact.
    copy(native/'SHA256SUMS',output/'native/SHA256SUMS')
    for cell in (0,1,2,3,4):
        for extension in ('parity','stdout','stderr'):
            name = 'checks/native-'+str(cell)+'-NNNN-1-plain.'+extension
            copy(native/name,output/'native'/name)
    for name in ('build_tower_encoder.py','run_tower_encoder_checks.py','run_tower_kernel_checks.py',
                 'replay_tower_encoder.py','test_tower_encoder_replay.py','tower_encoder_overlay.py',
                 'retain_tower_encoder.py'):
        copy(Path(__file__).parent/name,output/'tools'/name)
    # Earlier source/build/failure records are retained as such, never relabeled
    # successful; large obsolete binaries/parity remain at the original roots.
    for label, root in (('initial','/tmp/leopard-tower-encoder.JfYtrY'),
                        ('overstrict-alias','/tmp/leopard-tower-encoder-final.Q23AuN')):
        source = Path(root)
        for name in ('build.json','build.log'):
            copy(source/name,output/'failures'/label/name)
        for path in (source/'source').iterdir(): copy(path,output/'failures'/label/'source'/path.name)
        checks = source/('initial-checks' if label == 'initial' else 'checks')
        for path in checks.iterdir():
            if path.suffix != '.parity': copy(path,output/'failures'/label/'checks'/path.name)
    lines = []
    for path in sorted(output.rglob('*')):
        if path.is_file():
            with path.open('rb') as stream:
                digest = hashlib.file_digest(stream,'sha256').hexdigest()
                os.posix_fadvise(stream.fileno(),0,0,os.POSIX_FADV_DONTNEED)
            lines.append(digest+'  '+str(path.relative_to(output)))
            path.chmod(0o444)
    (output/'SHA256SUMS').write_text('\n'.join(lines)+'\n')
    (output/'SHA256SUMS').chmod(0o444)
    for path in sorted(output.rglob('*'),reverse=True):
        if path.is_dir(): path.chmod(0o555)
    output.chmod(0o555)
    print(json.dumps(dict(files=len(lines)+1,bytes=sum(p.stat().st_size for p in output.rglob('*') if p.is_file()),
                         manifest_sha256=hashlib.sha256((output/'SHA256SUMS').read_bytes()).hexdigest())))


if __name__ == '__main__':
    if len(sys.argv) != 4: raise SystemExit('usage: retain_tower_encoder.py RAW NATIVE NEW_OUTPUT')
    retain(Path(sys.argv[1]).resolve(strict=True),Path(sys.argv[2]).resolve(strict=True),Path(sys.argv[3]).absolute())

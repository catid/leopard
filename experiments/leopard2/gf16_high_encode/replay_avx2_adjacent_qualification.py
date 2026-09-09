#!/usr/bin/env python3
"""Collector-free raw/native/parity replay; never execute an input codec.

leopard-79h.38.5.4.18.3. Only compile-time mode3 is semantically qualified.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
from avx2_adjacent_counts import validate
from audit_avx2_adjacent_schedule import audit as codegen_audit
from audit_avx2_pair_schedule import require, resource

REFERENCE = Path('/home/catid/leopard/.research/leopard-79h/avx2-pair-screen.m_vukpir')
ORIGINAL_ARCHIVES = {
    'release':'d8eb0985e7347e491ea069e2e39d64ec256a35808422d02c1df36641fb315c8a',
    'sanitize':'e2d1bacb1142f71f53ab8aaf927dc0254c37fa9abfdea89727633faf634724a6',
}


def same(a, b):
    require(json.dumps(a,sort_keys=True,allow_nan=False) ==
            json.dumps(b,sort_keys=True,allow_nan=False), 'value/type mismatch')


def read(path):
    return json.loads(path.read_text())


def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream,'sha256').hexdigest()


def parity_labels():
    return [p+'-'+n+'-'+str(c) for p in ('release','sanitize')
            for c in range(8) for n in ('no-clock','observer')]


def public_identity(label, record, native, codec_identity):
    same(record['cell'],int(label.rsplit('-',1)[1]))
    same([record['mode'],record['traced'],record['encode_calls'],record['samples_ns'],record['route_counts']],
         ['off',False,1,[],[0,0,0,0]])
    for key in ('cell','k','r','bytes','field','requested_backend','input_hash','output_hash','outer_guards'):
        same(record[key],native[key])
    same(record['codec_commit'],codec_identity)


def replay(root):
    state = read(root/'native/report.json')
    build = read(root/'checks-build/build.json')
    same([state['bead'],state['completed'],state['timed'],state['mode']],
         ['leopard-79h.38.5.4.18.3',True,False,3])
    same([build['bead'],build['timed'],build['mode']],[state['bead'],False,3])
    original_root = Path(next(iter(state['pins']))).parent.parent
    def locate(value):
        path = Path(value)
        return root/path.relative_to(original_root) if path.is_relative_to(original_root) else path
    for path, digest in state['pins'].items():
        same(sha(locate(path)),digest)
    for name, digest in build['driver_sources'].items():
        same(sha(root/'checks-build/drivers'/name),digest)
    codegen = codegen_audit(root)
    same(sha(REFERENCE/'SHA256SUMS'),'da76c223dacb2977699f39ab045150323b30ee020ab9b4cddcb19f5935495c46')
    reference_pins = dict((name,digest) for digest,name in
                         (line.split('  ',1) for line in (REFERENCE/'SHA256SUMS').read_text().splitlines()))
    same(sorted(build['profiles']),['release','sanitize'])
    labels = []
    expected_codes = {}
    for profile in ('release','sanitize'):
        selectors = ['adjacent','forward-ranges','pairs','split',*map(str,range(8)),
                     'roundtrip','concurrent','observer-unit']
        selectors += [name+'-'+str(cell) for cell in range(8) for name in ('no-clock','observer')]
        for name in selectors:
            label = profile+'-'+name
            labels.append(label); expected_codes[label] = 0
        label = profile+'-clock-guard'
        labels.append(label); expected_codes[label] = 86
        for name in [*('bad-'+str(i) for i in range(4)),*('observer-bad-'+str(i) for i in range(6))]:
            label = profile+'-'+name
            labels.append(label); expected_codes[label] = 1
    same([r['label'] for r in state['records']],labels)
    same(len(labels),84)
    for row in state['records']:
        label = row['label']
        same([row['returncode'],row['expected']],[expected_codes[label]]*2)
        for suffix in ('stdout','stderr'):
            same(sha(root/'native'/(label+'.'+suffix)),row[suffix+'_sha256'])
        stderr = (root/'native'/(label+'.stderr')).read_text()
        if label in state['observations']:
            same(json.loads(stderr),state['observations'][label])
        elif expected_codes[label] == 0:
            require(not stderr,'unexpected native stderr')
        elif expected_codes[label] == 86:
            same(stderr,'unexpected driver benchmark clock\n')
    same(sorted(state['observations']),sorted(p+'-observer-'+str(i) for p in ('release','sanitize') for i in range(8)))
    counts = []
    for cell in range(8):
        release = state['observations']['release-observer-'+str(cell)]
        same(release,state['observations']['sanitize-observer-'+str(cell)])
        counts.append(dict(cell=cell,**validate(release,cell)))
    for profile in ('release','sanitize'):
        prefix = root/'native'/profile
        def output(name):
            return Path(str(prefix)+'-'+name+'.stdout').read_text()
        same(output('adjacent'),'adjacent pairs: forward=65535 accumulating=65535 boundary_accumulations=918\n')
        same(output('forward-ranges'),'forward range cases: 384\n')
        same(output('pairs'),'pair kernel cases: 66147\n')
        same(output('split'),'split range cases: 256\n')
        same(output('roundtrip'),''.join('roundtrip cell '+str(i)+' passed\n' for i in range(3)))
        lines = output('concurrent').splitlines()
        same([len(lines),lines.count('roundtrip cell 1 passed'),lines.count('roundtrip cell 2 passed'),
              lines.count('four-thread both-field roundtrips passed')],[18,9,8,1])
        require('16 exact delegations' in output('observer-unit'),'observer delegation unit')
        for cell in range(8):
            guarded = json.loads(output(str(cell)))
            same([guarded['cell'],guarded['subset_masks'],guarded['timed']],[cell,6,False])
        entry = build['profiles'][profile]
        same(entry['original_sha256'],ORIGINAL_ARCHIVES[profile])
        if profile == 'release':
            same(entry['object_sha256'],codegen['modes']['3']['object_sha256'])
        original = Path(entry['original'])
        archive = root/'checks-build'/profile/'candidate.a'
        same(sha(original),entry['original_sha256'])
        same(sha(archive),entry['archive_sha256'])
        members = subprocess.check_output(['ar','t',str(original)],text=True).splitlines()
        same(subprocess.check_output(['ar','t',str(archive)],text=True).splitlines(),members)
        same([len(members),len(set(members)),len(entry['unchanged_members'])],[24,24,23])
        same(sha(root/'checks-build'/profile/'Leopard2BackendAVX2.cpp.o'),entry['object_sha256'])
        for member in members:
            actual = subprocess.check_output(['ar','p',str(archive),member])
            if member == 'Leopard2BackendAVX2.cpp.o':
                same(hashlib.sha256(actual).hexdigest(),entry['object_sha256'])
            else:
                require(actual == subprocess.check_output(['ar','p',str(original),member]),'unrelated archive member')
                same(hashlib.sha256(actual).hexdigest(),entry['unchanged_members'][member])
        for name, digest in entry['files'].items():
            same(sha(root/'checks-build'/profile/name),digest)
    total = 0
    same(len(state['parity_comparisons']),32)
    for row,label in zip(state['parity_comparisons'],parity_labels()):
        a,b = locate(row['path']),locate(row['original'])
        cell = int(label.rsplit('-',1)[1])
        same(str(a),str(root/'native'/(label+'.parity')))
        same(str(b),str(REFERENCE/'preparation'/f'native-native-{cell}.parity'))
        same(row['original_sha256'],reference_pins[f'preparation/native-native-{cell}.parity'])
        same(sha(a),row['sha256']); same(sha(b),row['original_sha256'])
        same([a.stat().st_size,b.stat().st_size],[row['bytes']]*2)
        with a.open('rb') as left,b.open('rb') as right:
            while True:
                chunk = left.read(65536)
                require(chunk == right.read(65536),'full native Leopard1 parity mismatch')
                if not chunk: break
                total += len(chunk)
            if a.is_relative_to(root):
                os.posix_fadvise(left.fileno(),0,0,os.POSIX_FADV_DONTNEED)
    same(total,217712384); same(total,state['parity_bytes'])
    same(sorted(state['public_records']),sorted(parity_labels()))
    for label, record in state['public_records'].items():
        same(read(root/'native'/(label+'.stdout')),record)
        cell = int(label.rsplit('-',1)[1])
        native_path = REFERENCE/'preparation'/f'native-native-{cell}-check.stdout'
        same(sha(native_path),reference_pins[f'preparation/native-native-{cell}-check.stdout'])
        native = read(native_path)
        profile = label.split('-')[0]
        digest = build['profiles'][profile]['original_sha256' if '-observer-' in label else 'archive_sha256']
        public_identity(label,record,native,('production:' if '-observer-' in label else 'adjacent-mode3:')+digest)
    peaks = {name:resource(root/(name+'.log'),536870912 if name=='checks-build' else 268435456)
             for name in ('checks-build','native','qualification-delivery-unit','qualification-delivery-unit-opt')}
    return dict(bead=state['bead'],timed=False,mode=3,positive_native_records=62,records=84,
        full_native_leopard1_parity_bytes=total,production_changed=False,observed_production_counts=counts,
        memory_peaks=peaks,archive_sha256={p:e['archive_sha256'] for p,e in build['profiles'].items()})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('workspace',type=Path)
    print(json.dumps(replay(parser.parse_args().workspace.resolve(strict=True)),indent=2))

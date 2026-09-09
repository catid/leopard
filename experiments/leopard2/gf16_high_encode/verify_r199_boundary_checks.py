#!/usr/bin/env python3
"""Collector-free qualification replay; no codec execution. leopard-79h.38.5.4.19."""
import hashlib
import json
import os
from pathlib import Path
import sys
from audit_avx2_pair_schedule import require, resource
from verify_gf16_callback_probe import validate_shape_counts

RAW = Path('/tmp/leopard-r199-boundary.v13Wsv')
REFERENCE = Path('/home/catid/leopard/.research/leopard-79h/avx2-pair-screen.m_vukpir')
ARCHIVES = {
    'native':'3f83c55599cbfbeea5d15adf3627405bcde1cbed2ddf4a672fbc75b186f5c4c1',
    'release':'d8eb0985e7347e491ea069e2e39d64ec256a35808422d02c1df36641fb315c8a',
    'sanitize':'e2d1bacb1142f71f53ab8aaf927dc0254c37fa9abfdea89727633faf634724a6',
}


def same(a,b):
    require(json.dumps(a,sort_keys=True,allow_nan=False) ==
            json.dumps(b,sort_keys=True,allow_nan=False),'value/type mismatch')


def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream,'sha256').hexdigest()


def read(path):
    return json.loads(path.read_text())


def codec_identity(profile):
    return ('native:6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198:' if profile == 'native'
            else 'production:') + ARCHIVES[profile]


def public_record(record,profile,variant,exercise=False):
    require(profile in ARCHIVES,'profile')
    native = profile == 'native'
    require(variant in (('native',) if native else ('auto','gfni')),'variant')
    same(record,dict(schema='leopard-r199-boundary/v1',profile='native' if native else 'current',
        codec_commit=codec_identity(profile),cell=0,k=1000,r=199,bytes=32768,requested=variant,
        execution_route='native_l1' if native else 'avx2' if variant=='auto' else 'gfni',
        scratch_bytes=16777216 if native else 16808512,
        output_semantics='first_r_work_buffers' if native else 'separate_output_buffers',
        input_hash='8b78decd04e67d27',output_hash='00648d9dbf0f2b20',outer_guards=True,
        encode_calls=26 if exercise else 1,samples_ns=[]))


def inventory():
    rows = []
    for profile in ('release','sanitize'):
        for legacy in (False,True):
            rows.extend((profile+('-original-' if legacy else '-focused-')+str(cell),0)
                        for cell in range(8))
    for profile in ARCHIVES:
        for variant in (('native',) if profile=='native' else ('auto','gfni')):
            prefix = profile+'-'+variant
            rows.extend((prefix+'-'+kind,0) for kind in ('check','exercise'))
            if profile != 'sanitize': rows.append((prefix+'-plain',0))
            if profile != 'native': rows.append((prefix+'-observer',0))
            rows.append((prefix+'-clock-guard',86))
        rows.extend((profile+'-bad-'+str(i),1) for i in range(6))
        if profile != 'native':
            rows.extend((profile+'-observer-bad-'+str(i),1) for i in range(5))
    return rows


def replay(root):
    state,build = read(root/'preparation/preparation.json'),read(root/'build/build.json')
    for value in (state,build):
        same([value['bead'],value['completed'],value['timed'],value['production_changed']],
             ['leopard-79h.38.5.4.19',True,False,False])
    def locate(value):
        path = Path(value)
        return root/path.relative_to(RAW) if path.is_relative_to(RAW) else path
    for name,digest in state['pins'].items(): same(sha(locate(name)),digest)
    same(sorted(build['profiles']),sorted(ARCHIVES))
    for name,digest in build['sources'].items(): same(sha(root/'build/drivers'/name),digest)
    for name,digest in build['headers'].items(): same(sha(root/'build/production-headers'/name),digest)
    for profile,entry in build['profiles'].items():
        same(entry['archive_sha256'],ARCHIVES[profile])
        same(entry['codec_identity'],codec_identity(profile))
        same(sha(Path(entry['original'])),ARCHIVES[profile])
        same(sha(root/'build'/profile/'codec.a'),ARCHIVES[profile])
        files = ['codec.a','no-clock','screen.o']
        if profile != 'sanitize': files.append('screen')
        if profile != 'native': files.extend(('observer','focused'))
        same(sorted(entry['files']),sorted(files))
        for name,digest in entry['files'].items(): same(sha(root/'build'/profile/name),digest)
    expected = inventory()
    same(len(expected),82)
    same([row['label'] for row in state['records']],[label for label,_ in expected])
    for row,(label,code) in zip(state['records'],expected):
        same([row['returncode'],row['expected']],[code,code])
        for suffix in ('stdout','stderr'):
            same(sha(root/'preparation'/(label+'.'+suffix)),row[suffix+'_sha256'])
        stderr = (root/'preparation'/(label+'.stderr')).read_text()
        if label in state['observations']:
            same(json.loads(stderr),state['observations'][label])
        elif code == 0: same(stderr,'')
        elif code == 86: same(stderr,'unexpected driver benchmark clock\n')
    public_labels = []
    parity_labels = []
    for profile in ARCHIVES:
        for variant in (('native',) if profile=='native' else ('auto','gfni')):
            prefix = profile+'-'+variant
            kinds = ['check','exercise']
            if profile != 'sanitize': kinds.append('plain')
            if profile != 'native': kinds.append('observer')
            for kind in kinds:
                label = prefix+'-'+kind
                record = read(root/'preparation'/(label+'.stdout'))
                public_record(record,profile,variant,kind=='exercise')
                same(record,state['public_records'][label])
                public_labels.append(label)
            parity_labels.append(prefix)
            if profile != 'native': parity_labels.append(prefix+'-observer')
    same(sorted(state['public_records']),sorted(public_labels))
    same(len(public_labels),17)
    for legacy in (False,True):
        for cell in range(8):
            name = ('original-' if legacy else 'focused-')+str(cell)
            actual = read(root/'preparation'/('release-'+name+'.stdout'))
            same(actual,read(root/'preparation'/('sanitize-'+name+'.stdout')))
            k = 1000 if cell < 6 else 17
            if legacy:
                r = (199 if cell%2 else 200) if cell < 6 else 7
                size = ((65536 if cell%2 else 32768)+(2 if cell>=4 else 0)) if cell<6 else (65 if cell==6 else 66)
                offset = 0 if cell<2 else 1 if cell<4 or cell>=6 else 2
            else:
                r = 199 if cell < 6 else 7
                size = (32768,32768,32770,32766,64,66,65,66)[cell]
                offset = (0,1,2,1,1,2,1,1)[cell]
            require(type(actual['scratch_bytes']) is int and actual['scratch_bytes']>0,'scratch size')
            same(actual,dict(schema='leopard-gfni-boundary-guards/v1',cell=cell,k=k,r=r,bytes=size,
                misalignment=offset,field=1 if cell==6 else 2,subset_masks=6,
                scratch_bytes=actual['scratch_bytes'],timed=False))
    same(sorted(state['observations']),sorted(p+'-'+v+'-observer'
         for p in ('release','sanitize') for v in ('auto','gfni')))
    counts = {}
    for variant,kind in (('auto',3),('gfni',6)):
        record = state['observations']['release-'+variant+'-observer']
        same(record,state['observations']['sanitize-'+variant+'-observer'])
        counts[variant] = validate_shape_counts(record,(1000,199,32768,kind,32768,1))
    same(sha(REFERENCE/'SHA256SUMS'),'da76c223dacb2977699f39ab045150323b30ee020ab9b4cddcb19f5935495c46')
    refs = dict((name,digest) for digest,name in
                (line.split('  ',1) for line in (REFERENCE/'SHA256SUMS').read_text().splitlines()))
    old = REFERENCE/'preparation/native-native-4.parity'
    same(sha(old),refs['preparation/native-native-4.parity'])
    same(len(state['comparisons']),9)
    total = 0
    for row,label in zip(state['comparisons'],parity_labels):
        path = root/'preparation'/(label+'.parity')
        reference = old if label=='native-native' else root/'preparation/native-native.parity'
        same(str(locate(row['path'])),str(path))
        same(str(locate(row['original'])),str(reference))
        same([path.stat().st_size,reference.stat().st_size,row['bytes']],[199*32768]*3)
        same(sha(path),row['sha256']); same(sha(reference),row['original_sha256'])
        with path.open('rb') as a,reference.open('rb') as b:
            while True:
                data = a.read(65536)
                require(data == b.read(65536),'full native parity mismatch')
                if not data: break
                total += len(data)
            os.posix_fadvise(a.fileno(),0,0,os.POSIX_FADV_DONTNEED)
    same(total,58687488); same(total,state['parity_bytes'])
    return dict(bead=state['bead'],positive_native_records=49,total_records=82,
        full_parity_comparisons=9,full_native_leopard1_parity_bytes=total,
        public_records=17,clock_guards=5,malformed_cli_refusals=28,
        guarded_shapes_per_profile=8,original_guarded_shapes_per_profile=8,
        production_changed=False,timed=False,observed_operations=counts,
        memory_peaks={name:resource(root/(name+'.log'),536870912 if name=='build' else 268435456)
                      for name in ('build','preparation')})


if __name__ == '__main__':
    require(len(sys.argv)==2,'usage: verify_r199_boundary_checks.py ROOT')
    print(json.dumps(replay(Path(sys.argv[1]).resolve()),indent=2))

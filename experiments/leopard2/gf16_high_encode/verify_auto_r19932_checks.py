#!/usr/bin/env python3
"""Read-only, collector-free qualification replay; leopard-79h.38.5.4.19.1.

Does not execute a codec or make a performance/promotion decision. The first
native scope hit memory.max without OOM; preserve that fact, not an all-zero
resource claim. Build records are pinned independently of the collector.
"""
import json
import os
from pathlib import Path
import re
import sys
import xml.etree.ElementTree as ET
import hashlib
from verify_auto_gfni_boundary_checks import require, same, sha, archive_members

RAW = Path('/tmp/leopard-auto-r19932.9EGR0p')
REFERENCE = Path('/home/catid/leopard/.research/leopard-79h/avx2-pair-screen.m_vukpir')
ARCHIVES = {
    'native': '3f83c55599cbfbeea5d15adf3627405bcde1cbed2ddf4a672fbc75b186f5c4c1',
    'release': '89f33d3d43f792de078c38469f9406dded2b69126af2c6fc8db6e2acb77cd334',
    'sanitize': 'c2aee8119a36bc647a450649106656437c55ca83fe08337a559edea6c533b0a9',
}
CELLS = [(1000,199,32768),(1000,199,32768),(1000,200,32768),
         (1000,199,65536),(1000,200,65536),(1000,198,32768),
         (1000,199,32768),(4096,512,4096),(17,7,64)]
INPUTS = ['8b78decd04e67d27']*3 + ['ba910426bb91dec4']*2 + \
         ['8b78decd04e67d27']*2 + ['8ab693ec532fd3fa','bef540029ca6f65b']
OUTPUTS = ['00648d9dbf0f2b20','00648d9dbf0f2b20','50d5c45deb48367a',
           '0b985bca4df709d5','b1d946fa500bbf12','9c2ebf307806ce91',
           '00648d9dbf0f2b20','12c692cff473f0ee','83958df25bd2c455']


def public_expected(profile, mode, cell, exercise=False):
    require(profile in ARCHIVES and type(mode) is int and mode in (0,1), 'profile/mode')
    require(type(cell) is int and 0 <= cell < 9, 'cell')
    native = profile == 'native'
    require(not native or mode == 0, 'native mode')
    gfni = 2 <= cell <= 4 or (cell < 2 and mode == 1)
    k,r,size = CELLS[cell]
    scratch = ([16777216]*3 + [33554432]*2 + [16777216]*2 + [4194304,1024]
               if native else [16808512]*7 + [4308992,1728])[cell]
    return dict(schema='leopard-auto-r19932-screen/v1',
        implementation='leopard1' if native else 'leopard2',
        codec_commit=('native:6e5725eb:' if native else 'auto-r19932:')+ARCHIVES[profile],
        cell=cell,k=k,r=r,bytes=size,boundary_mode=-1 if native else mode,
        api='leo_encode' if native else 'leo2_encode_batch_one_item' if cell==1 else 'leo2_encode',
        execution_route='native-leopard1' if native else 'gf8-auto' if cell==8 else 'gfni' if gfni else 'avx2',
        field=1 if cell==8 else 2,encode_calls=26 if exercise else 1,
        untimed_route_calls=0 if native else int(gfni),scratch_bytes=scratch,
        output_semantics='first_r_work_buffers' if native else 'separate_output_buffers',
        input_hash=INPUTS[cell],output_hash=OUTPUTS[cell],samples_ns=[])


def inventory():
    rows = []
    for p in ('release','sanitize'):
        rows += [(p+'-new-'+s,0) for s in ('routes','api','concurrent')]
        rows += [(p+'-new-fault-'+s,0) for s in ('host','unavailable','oom','kat')]
        rows += [(f'{p}-new-guard-{m}-{c}',0) for m in range(2) for c in range(8)]
        rows.append((p+'-old-routes',0))
        for c in range(2):
            rows += [(f'{p}-old-{s}-{c}',0) for s in ('api','concurrent')]
            rows += [(f'{p}-old-fault-{c}-{s}',0) for s in ('host','unavailable','oom','kat')]
        rows += [(f'{p}-old-guard-{m}-{c}',0) for m in range(2) for c in range(8)]
        rows.append((p+'-production',0))
        rows += [(f'{p}-focused-bad-{i}',1) for i in range(6)]
    for p in ARCHIVES:
        for m in ((0,) if p=='native' else (0,1)):
            for c in range(9):
                rows += [(f'{p}-{m}-{c}-{s}',0) for s in
                         (('check','exercise') if p=='sanitize' else ('check','exercise','plain'))]
            rows.append((f'{p}-{m}-clock-guard',86))
        rows += [(f'{p}-screen-bad-{i}',1) for i in range(6)]
        if p=='native': rows.append(('native-no-candidate-mode',1))
    return rows


def guard_expected(old, cell):
    sizes = [(32768,65536,32768,65536,32770,65538,65,66),
             (32768,32768,32770,32766,64,66,65,66)][not old]
    offsets = [(0,0,1,1,2,2,1,1),(0,1,2,1,1,2,1,1)][not old]
    scratch = [(16808512,)*4+(16872512,)*2+(2816,)*2,
               (16808512,16808512,16872512,16839744,64064,128064,2816,2816)][not old]
    return dict(schema='leopard-gfni-boundary-guards/v1',cell=cell,
        k=1000 if cell<6 else 17,r=(199 if not old or cell%2 else 200) if cell<6 else 7,
        bytes=sizes[cell],misalignment=offsets[cell],field=1 if cell==6 else 2,
        subset_masks=6,scratch_bytes=scratch[cell],timed=False)


def arguments(label):
    p,*parts = label.split('-')
    directory = RAW/'build'/p
    if parts[0] in ('new','old'):
        old = parts[0]=='old'
        kind = parts[1]
        argv = [str(directory/('original-check' if old else 'focused'))]
        if kind=='guard': return argv+['--guards',parts[3],parts[2]]
        if kind=='fault': return argv+['--fault',*parts[2:]]
        return argv+['--'+kind,*parts[2:]]
    if parts[0]=='production': return [str(directory/'production-check')]
    if parts[0]=='focused':
        bad = [[],['--api','extra'],['--fault','bad'],['--guards'],['--guards','8','0'],['--guards','0','2']]
        return [str(directory/'focused'),*bad[int(parts[-1])]]
    if parts[0]=='screen':
        bad = [[],['--check'],['--measure','0','0','forbidden-output'],['--check','9','0'],
               ['--check','0','2'],['--check','00','0']]
        return [str(directory/'no-clock'),*bad[int(parts[-1])]]
    if label=='native-no-candidate-mode': return [str(directory/'no-clock'),'--check','0','1']
    if parts[1]=='clock': return [str(directory/'no-clock'),'--measure','0',parts[0]]
    mode,cell,kind = parts
    argv = [str(directory/('screen' if kind=='plain' else 'no-clock')),
            '--check' if kind=='plain' else '--'+kind,cell,mode]
    if kind=='check': argv.append(str(RAW/'preparation'/f'{p}-{mode}-{cell}.parity'))
    return argv


def scope(text, maximum, max_events=0):
    def value(name):
        matches = re.findall(r'^'+re.escape(name)+r'\n([0-9]+)$',text,re.M)
        require(len(matches)==1, 'resource field: '+name)
        return int(matches[0])
    peak = value('memory.peak')
    require(0 < peak <= maximum and value('memory.max') == maximum,'memory bound')
    same(text.split('memory.events\n')[1],
         f'low 0\nhigh 0\nmax {max_events}\noom 0\noom_kill 0\noom_group_kill 0\n'
         'memory.swap.current\n0\nmemory.swap.max\n0\n')
    require('\tExit status: 0\n' in text,'scope exit')
    return dict(peak=peak,maximum=maximum,max_events=max_events,oom_events=0,swap=0)


def read(path):
    return json.loads(path.read_text())


def hook_members(path):
    """Independent streamed ar inventory for the actual 23-member hook build."""
    result, names = {}, b''
    with path.open('rb') as f:
        require(f.read(8)==b'!<arch>\n','hook ar format')
        while True:
            header = f.read(60)
            if not header: break
            require(len(header)==60 and header[58:]==b'`\n','hook ar header')
            size = int(header[48:58]); name = header[:16].decode('ascii').strip()
            require(0 <= size < 64*1024*1024,'hook member size')
            digest = hashlib.sha256(); remaining = size; chunks = []
            while remaining:
                chunk = f.read(min(65536,remaining))
                require(bool(chunk),'short hook archive')
                digest.update(chunk); remaining -= len(chunk)
                if name=='//': chunks.append(chunk)
            if size%2: require(f.read(1)==b'\n','hook ar padding')
            if name=='//': names = b''.join(chunks); continue
            if name in ('/','/SYM64/'): continue
            if name.startswith('/'):
                offset = int(name[1:]); require(0 <= offset < len(names),'hook long name')
                name = names[offset:].split(b'/\n',1)[0].decode('ascii')
            else: name = name.rstrip('/')
            require(name and name==Path(name).name and name not in result,'hook member name')
            result[name] = digest.hexdigest()
    require(len(result)==23,'hook archive has 23 members')
    return result


def replay(root):
    same(sha(root/'build/build.json'),'52af5c02bc1040bfe822d270159cfd9bfc704d72ed662a6678bc00b7d70dffeb')
    same(sha(root/'hooks/finished.json'),'e538981c16ecddd0070f326ab89c12a6f87413a0a62d12a43a95e1f6a9533076')
    build,state = read(root/'build/build.json'),read(root/'preparation/preparation.json')
    def locate(name):
        path = Path(name)
        return root/path.relative_to(RAW) if path.is_relative_to(RAW) else path
    for value in (build,state):
        same([value['bead'],value['completed'],value['timed'],value['default_enabled']],
             ['leopard-79h.38.5.4.19.1',True,False,False])
    expected_pins = {}
    for group,folder in (('source_pins','source'),('driver_pins','drivers')):
        for name,digest in build[group].items():
            same(sha(root/'build'/folder/name),digest)
            expected_pins[str(RAW/'build'/folder/name)] = digest
    for p,entry in build['profiles'].items():
        expected_pins[entry['original']] = entry['original_sha256']
        same(sha(Path(entry['original'])),entry['original_sha256'])
        for name,digest in entry['files'].items():
            expected_pins[str(RAW/'build'/p/name)] = digest
            same(sha(root/'build'/p/name),digest)
        same(sha(root/'build'/p/'candidate.a'),ARCHIVES[p])
        if p!='native':
            old,new = archive_members(Path(entry['original'])),archive_members(root/'build'/p/'candidate.a')
            same([n for n in old if old[n]!=new[n]],['leopard2.cpp.o'])
            same({n:d for n,d in new.items() if n!='leopard2.cpp.o'},entry['unchanged_members'])
            same(new['leopard2.cpp.o'],sha(root/'build'/p/'leopard2.cpp.o'))
    same(state['pins'],expected_pins)
    rows = inventory()
    same([r['label'] for r in state['records']],[n for n,_ in rows])
    same(sorted(p.stem for p in (root/'preparation').glob('*.stdout')),sorted(n for n,_ in rows))
    for row,(label,code) in zip(state['records'],rows):
        same([row['returncode'],row['expected']],[code,code])
        same(row['argv'],arguments(label))
        for suffix in ('stdout','stderr'):
            same(sha(root/'preparation'/(label+'.'+suffix)),row[suffix+'_sha256'])
        stderr = (root/'preparation'/(label+'.stderr')).read_text()
        if code==0: same(stderr,'')
        elif code==86: same(stderr,'unexpected driver benchmark clock\n')
        else: require(bool(stderr),'CLI rejection diagnostic')
        if code or label in state['public_records']: continue
        lines = (root/'preparation'/(label+'.stdout')).read_text().splitlines()
        if label.endswith('-production'):
            same(lines,[f'Production AUTO GF16 GFNI route passed: R={r} bytes={b}'
                        for r,b in ((200,65536),(200,32768),(199,65536))])
        else:
            parts = label.split('-'); old = parts[1]=='old'; case = parts[2]
            case = '--guards' if case=='guard' else '--'+case
            same(len(lines),2 if case=='--guards' else 1)
            same(json.loads(lines[-1]),dict(schema='leopard-auto-'+('gfni-boundary' if old else 'r19932')+
                '-check/v1',case=case,passed=True,timed=False))
            if case=='--guards': same(json.loads(lines[0]),guard_expected(old,int(parts[-1])))
    public_labels = []
    for p in ARCHIVES:
        for m in ((0,) if p=='native' else (0,1)):
            for c in range(9):
                for kind in (('check','exercise') if p=='sanitize' else ('check','exercise','plain')):
                    label = f'{p}-{m}-{c}-{kind}'; public_labels.append(label)
                    actual = read(root/'preparation'/(label+'.stdout'))
                    same(actual,public_expected(p,m,c,kind=='exercise'))
                    same(actual,state['public_records'][label])
    same(sorted(state['public_records']),sorted(public_labels))
    same(sha(REFERENCE/'SHA256SUMS'),'da76c223dacb2977699f39ab045150323b30ee020ab9b4cddcb19f5935495c46')
    refs = dict((name,digest) for digest,name in
                (line.split('  ',1) for line in (REFERENCE/'SHA256SUMS').read_text().splitlines()))
    older = {0:4,1:4,2:0,3:1,4:2,6:4,7:3,8:7}
    pairs = [(f'native-0-{c}',REFERENCE/'preparation'/f'native-native-{old}.parity',c)
             for c,old in older.items()]
    for _,ref,_ in pairs: same(sha(ref),refs[str(ref.relative_to(REFERENCE))])
    pairs += [(f'{p}-{m}-{c}',root/'preparation'/f'native-0-{c}.parity',c)
              for p in ('release','sanitize') for m in range(2) for c in range(9)]
    same(len(state['comparisons']),len(pairs))
    total = 0
    for row,(label,ref,c) in zip(state['comparisons'],pairs):
        path = root/'preparation'/(label+'.parity')
        same([str(locate(row['path'])),str(locate(row['original']))],[str(path),str(ref)])
        length = CELLS[c][1]*CELLS[c][2]
        same([row['bytes'],path.stat().st_size,ref.stat().st_size],[length]*3)
        same(sha(path),row['sha256']); same(sha(ref),row['original_sha256'])
        with path.open('rb') as a,ref.open('rb') as b:
            while True:
                data = a.read(65536)
                require(data==b.read(65536),'native Leopard1 parity bytes')
                if not data: break
                total += len(data)
            os.posix_fadvise(a.fileno(),0,0,os.POSIX_FADV_DONTNEED)
    same([total,state['parity_bytes']],[297765056]*2)
    hooks = read(root/'hooks/finished.json')
    for name,digest in hooks['files'].items(): same(sha(root/'hooks'/name),digest)
    same(hooks['source_sha256'],build['source_pins']['leopard2.cpp'])
    same(sha(Path(hooks['original'])),hooks['original_sha256'])
    before,after = hook_members(Path(hooks['original'])),hook_members(root/'hooks/candidate.a')
    same(sorted(before),sorted(after))
    same([name for name in before if before[name]!=after[name]],['leopard2.cpp.o'])
    same({name:d for name,d in after.items() if name!='leopard2.cpp.o'},hooks['unchanged_members'])
    same(after['leopard2.cpp.o'],sha(root/'hooks/leopard2.cpp.o'))
    checks = read(root/'hook-checks/checks.json')
    same([checks['completed'],checks['timed']],[True,False])
    hook_labels = [b+'-'+f for b in ('32','r199','r19932')
                   for f in ('kat','ff8-allocation','ff16-allocation')]
    hook_labels += ['disabled-inert','ineligible-inert','production','isa']
    same([r['label'] for r in checks['records']],hook_labels)
    for row in checks['records']:
        same(row['returncode'],0)
        name = row['label']
        for suffix in ('stdout','stderr'):
            same(sha(root/'hook-checks'/(name+'.'+suffix)),row[suffix+'_sha256'])
        same((root/'hook-checks'/(name+'.stderr')).read_text(),'')
        output = (root/'hook-checks'/(name+'.stdout')).read_text()
        require('skipped' not in output.lower() and bool(output),'hook/ISA output')
        if name.startswith(('32-','r199-','r19932-')):
            boundary,fault = name.split('-',1)
            r,size = (200,32768) if boundary=='32' else (199,65536) if boundary=='r199' else (199,32768)
            same(output,f'AUTO GF16 GFNI encode fallback passed: auto-gfni-encode-{fault}-fallback R={r} bytes={size}\n')
            same(row['argv'],[str(RAW/'hooks/check'),'auto-gfni-encode-'+fault+'-fallback','boundary-'+boundary])
        elif name=='isa':
            require('portable ISA check: PASS' in output,'ISA pass marker')
        else:
            same(row['argv'],[str(RAW/'hooks/check'),'auto-gfni-encode-'+name])
            require('passed' in output or 'stayed inert' in output,'backend pass marker')
    resources = {name:scope((root/(name+'.log')).read_text(),limit,events)
                 for name,limit,events in [('build',536870912,0),('preparation',268435456,612),
                     ('finish-hooks',536870912,0),('check-hooks',268435456,0),
                     ('cmake-configure',536870912,0),('ctest',268435456,0)]}
    ctest = read(root/'ctest/checks.json')
    same([ctest['completed'],ctest['timed'],ctest['full_cmake_library_build']],[True,False,False])
    same(ctest['executable_sha256'],hooks['files']['check'])
    same(sha(root/'cmake/leopard2_backend_failures_test'),hooks['files']['check'])
    cases = ET.parse(root/'ctest/results.xml').findall('.//testcase')
    names = ['leopard2_backend_auto_gfni_boundary_r19932_'+f for f in ('kat','ff8-allocation','ff16-allocation')]
    same(ctest['tests'],names)
    same(sorted(c.attrib['name'] for c in cases),sorted(names))
    for case in cases:
        require(case.find('failure') is None and case.find('skipped') is None,'CTest failure/skip')
        output = case.findtext('system-out','')
        require('fallback passed:' in output and 'R=199 bytes=32768' in output,'actual CTest output')
    return dict(bead=state['bead'],positive_native_records=223,total_records=259,
        public_records=117,clock_guards=5,cli_refusals=31,full_parity_comparisons=44,
        full_native_leopard1_parity_bytes=total,real_backend_checks=12,registered_ctests=3,release_isa_pass=True,
        archive_members=24,unchanged_members=23,hook_archive_members=23,hook_unchanged_members=22,
        default_enabled=False,timed=False,performance_qualified=False,production_promoted=False,
        resources=resources)


if __name__=='__main__':
    require(len(sys.argv)==2,'usage: verify_auto_r19932_checks.py ROOT')
    print(json.dumps(replay(Path(sys.argv[1]).resolve()),indent=2))

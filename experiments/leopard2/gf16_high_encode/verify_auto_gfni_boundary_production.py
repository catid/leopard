#!/usr/bin/env python3
"""Retained-only production integration audit; no codec execution or collector imports."""
import hashlib
import json
from pathlib import Path
import struct
import sys
import xml.etree.ElementTree as ET


def require(value, message):
    if not value: raise ValueError(message)


def same(a,b):
    require(json.dumps(a,sort_keys=True) == json.dumps(b,sort_keys=True),'typed identity')


def sha(path):
    with path.open('rb') as stream: return hashlib.file_digest(stream,'sha256').hexdigest()


def members(path):
    result, names = {}, b''
    with path.open('rb') as stream:
        require(stream.read(8) == b'!<arch>\n','ar magic')
        while True:
            header = stream.read(60)
            if not header: break
            require(len(header) == 60 and header[58:] == b'`\n','ar header')
            size = int(header[48:58])
            require(0 <= size < 64*1024*1024,'bounded member')
            payload = stream.read(size)
            require(len(payload) == size,'ar payload')
            if size % 2: require(stream.read(1) == b'\n','ar padding')
            name = header[:16].decode().strip()
            if name in ('/','/SYM64/'): continue
            if name == '//': names = payload; continue
            name = names[int(name[1:]):].split(b'/\n',1)[0].decode() if name.startswith('/') else name.rstrip('/')
            require(name not in result and Path(name).name == name,'member identity')
            result[name] = hashlib.sha256(payload).hexdigest()
    require(len(result) == 24,'complete dual-field archive')
    return result


def changed_bytes(before, after):
    require(before.stat().st_size == after.stat().st_size,'artifact length changed')
    changes,offset = [],0
    with before.open('rb') as a, after.open('rb') as b:
        while True:
            left,right = a.read(65536),b.read(65536)
            if not left: break
            if left != right:
                changes += [(offset+i,x,y) for i,(x,y) in enumerate(zip(left,right)) if x != y]
            require(len(changes) <= 1,'more than default-data byte changed')
            offset += len(left)
    return changes


def boundary_symbol_offset(path):
    data = path.read_bytes()
    require(data[:6] == b'\x7fELF\x02\x01','ELF64 little endian')
    offset, = struct.unpack_from('<Q',data,40)
    width,count,names_index = struct.unpack_from('<HHH',data,58)
    require(width == 64,'ELF section width')
    sections = [struct.unpack_from('<IIQQQQIIQQ',data,offset+i*width) for i in range(count)]
    section_names = sections[names_index]
    names = data[section_names[4]:section_names[4]+section_names[5]]
    matches = []
    for section in sections:
        if section[1] != 2: continue
        strings = sections[section[6]]
        strings = data[strings[4]:strings[4]+strings[5]]
        require(section[9] == 24,'ELF symbol width')
        for cursor in range(section[4],section[4]+section[5],24):
            name,info,other,index,value,size = struct.unpack_from('<IBBHQQ',data,cursor)
            label = strings[name:].split(b'\0',1)[0]
            if b'g_auto_gf16_gfni_boundary_mode' in label:
                target = sections[index]
                require(names[target[0]:].split(b'\0',1)[0] == b'.data' and size == 4,'boundary data symbol')
                matches.append(target[4]+value)
    require(len(matches) == 1,'one boundary symbol')
    return matches[0]


def resources(path, maximum):
    lines = path.read_text().splitlines()
    require(lines.count('memory.peak') == 1 and '\tExit status: 0' in lines,'scope completion')
    index = lines.index('memory.peak'); peak = int(lines[index+1])
    require(0 < peak <= maximum,'memory peak')
    same(lines[index:],['memory.peak',str(peak),'memory.max',str(maximum),'memory.events',
         'low 0','high 0','max 0','oom 0','oom_kill 0','oom_group_kill 0',
         'memory.swap.current','0','memory.swap.max','0'])
    return peak


def verify(root, candidate):
    for path,value in {
        root/'release/libleopard.a':'d8eb0985e7347e491ea069e2e39d64ec256a35808422d02c1df36641fb315c8a',
        candidate/'release/candidate.a':'80bffc9e873585d9a18fcf6a294413b8c2a76d84fa6e7f9a437fb96ca8458533',
        root/'sanitize/candidate.a':'e2d1bacb1142f71f53ab8aaf927dc0254c37fa9abfdea89727633faf634724a6',
        candidate/'sanitize/candidate.a':'9a87a6baa349195e02c810fab8131fe5d4edfaf86de5a0efd1b7cafedb207c38',
        root/'release/leopard2_auto_gf16_gfni_production_test':'5e3ab152301227e3c1ae5715e04af30062839b8929d43e82e10419e5f1345b2a',
        root/'release/leopard2_auto_encode_backend_test':'0fe2d27c573c69c9a076f27440bc9a4556c99769db0b86d9489496926c9cea38',
        root/'release/leopard2_backend_failures_test':'8fa3434828382d26de91d6af84a9328703d1972278a8c692909d5fffef947037',
        root/'release/check':'86c81a67d7aeb267b7d772acbbe6bd77152302cec48d6a0fd6e67f0946a9033b',
        root/'sanitize/check':'e136e1705797caec1c38260afaa16e97075ebd0540bec8fbaf89ae4164b3a7b4',
        root/'sanitize/production-check':'cc3612caf60ad366b0c02aca86c29297b30ac2a79c074e3bdac25a34bfb5a8d0',
        root/'focused-checks.log':'a83618ece3f1dbacc1040d3bf61959d0633e7b3b82afd5628bce46f0121dfad1',
        root/'final-tests.xml':'b97c6b9a8074ec3e6e52baea058335bb9b9711528c60b4e8b99b66ddb148aca2',
    }.items(): same(sha(path),value)
    changes = changed_bytes(candidate/'release/candidate.a',root/'release/libleopard.a')
    same(changes,[(285580,2,1)])
    core = root/'release/CMakeFiles/leopard.dir/leopard2.cpp.o'
    core_changes = changed_bytes(candidate/'release/leopard2.cpp.o',core)
    same(core_changes,[(boundary_symbol_offset(core),2,1)])
    for profile,path in (('release','release/libleopard.a'),('sanitize','sanitize/candidate.a')):
        before,after = members(candidate/profile/'candidate.a'),members(root/path)
        same(sorted(before),sorted(after))
        same([name for name in before if before[name] != after[name]],['leopard2.cpp.o'])
    count = 0
    for profile in ('release','sanitize'):
        expected = {'production':None,'routes':'--routes'}
        for cell in (0,1):
            for kind in ('api','concurrent'): expected[f'{kind}-{cell}'] = '--'+kind
            for fault in ('host','unavailable','oom','kat'): expected[f'fault-{cell}-{fault}'] = '--fault'
        for cell in range(8):
            for mode in (0,1): expected[f'guards-{cell}-{mode}'] = '--guards'
        for label,case in expected.items():
            prefix = root/'checks'/(profile+'-'+label)
            require(prefix.with_suffix('.stderr').stat().st_size == 0,'native stderr')
            lines = prefix.with_suffix('.stdout').read_text().splitlines()
            if case is None:
                same(lines,[f'Production AUTO GF16 GFNI route passed: R={r} bytes={size}'
                            for r,size in ((200,65536),(200,32768),(199,65536))])
            else:
                same(json.loads(lines[-1]),dict(schema='leopard-auto-gfni-boundary-check/v1',
                     case=case,passed=True,timed=False))
                require(len(lines) == (2 if case == '--guards' else 1),'record count')
                if case == '--guards':
                    # Compare every field to the already verified qualified candidate's exact same case.
                    old = (candidate/'checks'/(profile+'-'+label+'.stdout')).read_text().splitlines()
                    same(json.loads(lines[0]),json.loads(old[0]))
            count += 1
    require(len(list((root/'checks').glob('*.stdout'))) == count == 60,'native check inventory')
    xml = ET.parse(root/'final-tests.xml').getroot()
    same([xml.get(x) for x in ('tests','failures','disabled','skipped')],['13','0','0','0'])
    cases = xml.findall('testcase')
    require(len(cases) == 13 and all(x.get('status') == 'run' and
            not x.findall('failure') and 'skipped' not in x.findtext('system-out','').lower() for x in cases),
            'actual CTests, not skips')
    wanted = {f'leopard2_backend_auto_gfni_boundary_{boundary}_{fault}':
              f'auto-gfni-encode-{fault}-fallback R={r} bytes={size}'
              for boundary,r,size in (('32',200,32768),('r199',199,65536))
              for fault in ('kat','ff8-allocation','ff16-allocation')}
    outputs = {x.get('name'):x.findtext('system-out','') for x in cases}
    require(all(text in outputs[name] for name,text in wanted.items()),'real boundary fault cases')
    isa = ET.parse(root/'isa.xml').getroot()
    same([isa.get(x) for x in ('tests','failures','skipped')],['1','0','0'])
    require('oom-kill' in (root/'release-checks-oom-journal.log').read_text(),
            'failed initial scope must be retained')
    peaks = {name:resources(root/(name+'.log'),maximum) for name,maximum in (
        ('build-release',536870912),('build-sanitizer-checks',536870912),
        ('build-final-tests',536870912),('build-comment-refresh',536870912),
        ('focused-checks',268435456),('final-tests',268435456),('isa',268435456))}
    return dict(release_archive_changed_bytes=changes,boundary_symbol_offset=core_changes[0][0],
                release_executable_code_unchanged=True,archive_objects_per_build=24,
                unchanged_objects_per_build=23,native_checks=count,final_ctests=13,
                real_boundary_fault_cases=6,portable_isa=True,memory_peaks=peaks,
                initial_test_oom_retained=True,production_default_enabled=True,
                independent_model_converged=False,broader_goal_complete=False)


if __name__ == '__main__':
    require(len(sys.argv) == 3,'usage: verify_auto_gfni_boundary_production.py ROOT CANDIDATE_CHECKPOINT')
    print(json.dumps(verify(Path(sys.argv[1]),Path(sys.argv[2])),sort_keys=True))

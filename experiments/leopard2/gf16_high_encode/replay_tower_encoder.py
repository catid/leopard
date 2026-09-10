#!/usr/bin/env python3
"""Read-only, collector-free replay of the scoped tower encoder qualification."""
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    with path.open('rb') as stream:
        digest = hashlib.file_digest(stream, 'sha256').hexdigest()
        os.posix_fadvise(stream.fileno(),0,0,os.POSIX_FADV_DONTNEED)
        return digest


def equal(a, b):
    with a.open('rb') as left, b.open('rb') as right:
        while True:
            chunk = left.read(65536)
            require(chunk == right.read(65536), 'byte parity: '+str(a))
            if not chunk:
                os.posix_fadvise(left.fileno(),0,0,os.POSIX_FADV_DONTNEED)
                os.posix_fadvise(right.fileno(),0,0,os.POSIX_FADV_DONTNEED)
                return


def resource(lines):
    require(len(lines) == 15 and lines[0] == 'memory.peak', 'resource shape')
    peak = int(lines[1])
    require(0 < peak < 268435456 and lines[2:] == [
        'memory.max','268435456','memory.events','low 0','high 0','max 0','oom 0',
        'oom_kill 0','oom_group_kill 0','memory.swap.current','0','memory.swap.max','0'], 'resource envelope')
    return peak


def audit_isa(text, portable):
    instructions = 0
    for line in text.splitlines():
        match = re.match(r'^\s*[0-9a-f]+:\s+((?:[0-9a-f]{2}\s+)+)\s*(\S+)(.*)$', line)
        if not match:
            continue
        encoding, mnemonic, operands = match.groups()
        instructions += 1
        first = encoding.split()[0]
        require(first != '62' and not re.search(r'\b(?:zmm\d+|[xyz]mm(?:1[6-9]|2\d|3[01])|k[0-7])\b', operands),
                'unexpected wide ISA')
        require(not any(word in mnemonic for word in ('gf2p8','ternlog','pclmul')), 'excluded instruction family')
        if portable:
            require(first not in ('c4','c5') and not re.search(r'\bymm\d+\b',operands), 'AVX in baseline selector/control')
    require(instructions > 0, 'empty disassembly')
    return instructions


def replay(root, native):
    build = json.loads((root/'build.json').read_text())
    require(build['completed'] is True and build['timed'] is False, 'build incomplete')
    require(set(build['profiles']) == {'release','trace','sanitize'}, 'profiles')
    for name, digest in build['sources'].items():
        require(sha(root/'source'/name) == digest, 'source pin '+name)
    require(build['sources']['LeopardFF16.original.cpp'] ==
            'fd27e72cb09068d83bc32b23acf49c6d133382579f109d6cfa227361c43d8cbf', 'original field')
    for profile, data in build['profiles'].items():
        for name, digest in data['files'].items():
            require(sha(root/profile/name) == digest, 'artifact pin '+profile+'/'+name)
        require(len(data['unchanged_members']) == 23, 'unchanged archive members')
        require(data['original_sha256'] == ('c2aee8119a36bc647a450649106656437c55ca83fe08337a559edea6c533b0a9'
            if profile == 'sanitize' else '89f33d3d43f792de078c38469f9406dded2b69126af2c6fc8db6e2acb77cd334'), 'original archive')
        original = root/('original-sanitize.a' if profile == 'sanitize' else 'original-release.a')
        if not original.exists(): original = Path(data['original_archive'])
        require(sha(original) == data['original_sha256'], 'actual original archive')
        before = subprocess.check_output(['ar','t',str(original)],text=True,timeout=10).splitlines()
        archive = root/profile/'candidate.a'
        after = subprocess.check_output(['ar','t',str(archive)],text=True,timeout=10).splitlines()
        require(len(before) == len(set(before)) == 24 and after == before +
                ['tower_encoder.cpp.o','tower_butterfly_probe.cpp.o'], 'actual archive inventory')
        for name in before:
            if name == 'LeopardFF16.cpp.o': continue
            old = subprocess.check_output(['ar','p',str(original),name],timeout=10)
            require(hashlib.sha256(old).hexdigest() == data['unchanged_members'][name] and
                    old == subprocess.check_output(['ar','p',str(archive),name],timeout=10), 'actual unchanged member')
    require(sha(root/'release/tower_butterfly_probe.cpp.o') ==
            'fc363f48a317469117d04c5c68519335daffb9c1e46a5717aa3c10a25e435446', 'qualified exact kernel object')
    isa = {}
    for name in ('LeopardFF16.cpp.o','tower_encoder.cpp.o','tower_butterfly_probe.cpp.o'):
        saved = (root/'release'/(name+'.disassembly')).read_text()
        actual = subprocess.check_output(['objdump','-drwC','-Mintel',str(root/'release'/name)],text=True,timeout=10)
        normalize = lambda text: re.sub(r'^.*: +file format elf64-x86-64$', 'OBJECT: file format elf64-x86-64',text,flags=re.M)
        require(normalize(saved) == normalize(actual),'actual disassembly drift')
        isa[name] = audit_isa(actual,name != 'tower_butterfly_probe.cpp.o')
    checks = json.loads((root/'checks/checks.json').read_text())
    require(checks['completed'] is True and checks['timed'] is False and len(checks['results']) == 123 and
            checks['build_sha256'] == sha(root/'build.json'), 'public completion')
    by_key = {(r['profile'],r['shape']):r for r in checks['results']}
    require(len(by_key) == 123, 'duplicate check')
    peaks = []
    parity_bytes = 0
    for profile in ('release','trace','sanitize','original-release','original-sanitize'):
        original = profile.startswith('original-')
        for shape in range(21 if original else 27):
            result = by_key[profile,shape]
            prefix = root/'checks'/(profile+'-'+str(shape))
            require(result['exit_code'] == 0, 'public exit status')
            stderr = prefix.with_suffix('.stderr').read_text()
            require(re.fullmatch(r'Running as unit: [\w.-]+; invocation ID: [0-9a-f]+\n',stderr) is not None,
                    'unexpected public stderr')
            lines = prefix.with_suffix('.stdout').read_text().splitlines()
            if original:
                old = json.loads(lines.pop(0))
                require(old['cell'] == shape and old['subset_masks'] == 6 and old['timed'] is False, 'original guarded check')
            record = json.loads(lines.pop(0))
            require(record == result['native'] and record['timed'] is False and
                    record['trace'] is (not original and profile != 'release'), 'native record')
            peaks.append(resource(lines))
            require(peaks[-1] == result['memory_peak'], 'recorded peak')
            if shape < 21:
                require(record['shape'] == shape, 'shape identity')
                path = prefix.with_suffix('.parity')
                require(path.stat().st_size == record['r']*record['bytes'] and sha(path) == result['parity_sha256'], 'parity file')
                if profile != 'original-release':
                    equal(path,root/'checks'/('original-release-'+str(shape)+'.parity'))
                    parity_bytes += path.stat().st_size
                require(record['scratch_bytes'] == by_key['original-release',shape]['native']['scratch_bytes'], 'scratch unchanged')
                if not original:
                    require(record['masks'] == 6 and record['batch_calls'] == 2 and len(record['mask_counts']) == 6, 'public coverage')
                    c = record
                    prefix_bytes = c['bytes'] - c['bytes']%64
                    eligible_size = c['field'] == 2 and c['r'] > 128 and prefix_bytes > 16384
                    any_selected = False
                    prefixes = [c['r'],1,c['r']-1,c['r'],c['r'] if c['r']%2 else c['r']-1,0]
                    for mask, count in enumerate(c['mask_counts']):
                        selected = eligible_size and mask != 5 and (c['backend'] == 3 or (shape == 13 and mask != 0))
                        any_selected |= selected
                        if c['trace']:
                            require((count['passes'] > 0) == selected and
                                    count['source_bytes'] == (c['k']*prefix_bytes if selected else 0) and
                                    count['output_bytes'] == (prefixes[mask]*prefix_bytes if selected else 0), 'mask boundaries')
                        else: require(count == {'passes':0,'source_bytes':0,'output_bytes':0}, 'untraced counts')
                    require(c['initializations'] == int(any_selected), 'cache count')
            else:
                require(record['special'] == shape-21 and record['initializations'] == (0 if shape in (22,24) else 1), 'special case')
    # Use only the retained original native-L1 shapes with matching public
    # input seed/field/layout. Partial/tail cases instead have original-L2 controls.
    native_bytes = 0
    native_hashes = {}
    manifest = {line.split('  ',1)[1]:line.split('  ',1)[0] for line in (native/'SHA256SUMS').read_text().splitlines()}
    for shape, cell in ((0,0),(1,1),(2,2),(3,4),(12,3)):
        relative = 'checks/native-'+str(cell)+'-NNNN-1-plain.parity'
        path = native/relative
        require(sha(path) == manifest[relative], 'native parity manifest')
        stdout_relative = relative.removesuffix('.parity')+'.stdout'
        require(sha(native/stdout_relative) == manifest[stdout_relative], 'native metadata manifest')
        old = json.loads((native/stdout_relative).read_text())
        current = by_key['release',shape]['native']
        require(old['codec'] == 'native:3f83c55599cbfbeea5d15adf3627405bcde1cbed2ddf4a672fbc75b186f5c4c1' and
                old['timed'] is False and old['samples'] == [] and
                all(old[key] == current[key] for key in ('k','r','bytes')) and
                old['output_hash'] == current['parity_hash'], 'native workload identity')
        for profile in ('release','trace','sanitize','original-release','original-sanitize'):
            equal(path,root/'checks'/(profile+'-'+str(shape)+'.parity'))
            native_bytes += path.stat().st_size
        native_hashes[relative] = sha(path)
    kernels = json.loads((root/'kernel-checks/checks.json').read_text())
    require(kernels['completed'] is True and kernels['timed'] is False and len(kernels['results']) == 39 and
            kernels['build_sha256'] == sha(root/'build.json'), 'kernel completion')
    require([(r['profile'],r['case']) for r in kernels['results']] ==
            [(p,c) for p in ('release','trace','sanitize') for c in range(13)], 'kernel sequence')
    for result in kernels['results']:
        case, profile = result['case'], result['profile']
        require(result['exit_code'] == (0 if case == 0 else 134 if case < 10 else 1), 'kernel exit')
        prefix = root/'kernel-checks'/(profile+'-'+str(case))
        lines = prefix.with_suffix('.stdout').read_text().splitlines()
        if case == 0:
            r = json.loads(lines.pop(0))
            require(r['selector_cases'] == 9000 and r['basis_pairs'] == 2097152 and r['log_cases'] == 65536 and
                    r['initializations'] == 1 and r['timed'] is False, 'kernel counts')
        peaks.append(resource(lines))
        require(peaks[-1] == result['memory_peak'], 'kernel peak')
    return dict(tracker='leopard-79h.38.5.4.18.4',timed=False,public_processes=123,kernel_processes=39,
                cross_profile_parity_bytes=parity_bytes,native_parity_comparison_bytes=native_bytes,
                native_parity_pins=native_hashes,maximum_native_memory_peak=max(peaks),release_instruction_counts=isa,
                build_sha256=sha(root/'build.json'),checks_sha256=sha(root/'checks/checks.json'),
                kernel_checks_sha256=sha(root/'kernel-checks/checks.json'))


if __name__ == '__main__':
    require(len(sys.argv) == 3, 'usage: replay_tower_encoder.py BUILD NATIVE_BUNDLE')
    print(json.dumps(replay(Path(sys.argv[1]),Path(sys.argv[2])),sort_keys=True))

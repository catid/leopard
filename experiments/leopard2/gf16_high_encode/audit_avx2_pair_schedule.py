#!/usr/bin/env python3
"""Read-only codegen/native-evidence audit; never execute an input codec.

leopard-79h.38.5.4.18.2. Static counts are not dynamic work or timings.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess

FUNCTION = 'void leopard::backend::AVX2FF16Butterfly2<true>(void*, void*, unsigned short, unsigned long)'


def require(value, message):
    if not value:
        raise ValueError(message)


def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def instructions(text):
    result = []
    for line in text.splitlines():
        parts = line.split('\t')
        if len(parts) >= 3 and re.fullmatch(r'\s*[0-9a-f]+:', parts[0]):
            raw = bytes.fromhex(parts[1])
            assembly = parts[2].strip()
            require(raw and assembly, 'empty instruction')
            result.append((int(parts[0].strip()[:-1], 16), raw, assembly))
    require(result, 'no instructions')
    return result


def pair_loop(text):
    bodies = re.split(r'^([0-9a-f]+) <(.+)>:\n', text, flags=re.M)
    targets = [bodies[i + 2] for i in range(1, len(bodies), 3)
               if bodies[i + 1] == FUNCTION]
    require(len(targets) == 1, 'missing or duplicate target function')
    ops = instructions(targets[0])
    loops = []
    for address, raw, assembly in ops:
        branch = re.match(r'j(?!mp\b)[a-z]+\s+([0-9a-f]+)\s', assembly)
        if not branch or int(branch[1], 16) >= address:
            continue
        start = int(branch[1], 16)
        body = [op for op in ops if start <= op[0] <= address]
        if sum(op[2].split()[0] == 'vpshufb' for op in body) != 8:
            continue
        require(body[0][0] == start, 'branch target is not an instruction')
        stack = [op[2] for op in body if re.search(r'\([^)]*%(?:rsp|rbp)\b', op[2])]
        loops.append(dict(start=hex(start), end=hex(address), instructions=len(body),
                          byte_shuffles=8, stack_references=stack))
    require(len(loops) == 1, 'ambiguous or missing eight-shuffle loop')
    return loops[0]


def resource(path, limit):
    text = path.read_text()
    def value(name):
        matches = re.findall(r'^' + re.escape(name) + r'\n([0-9]+)$', text, re.M)
        require(len(matches) == 1, 'missing resource field: ' + name)
        return int(matches[0])
    require(value('memory.max') == limit and value('memory.peak') < limit, 'memory bound')
    require(value('memory.swap.current') == value('memory.swap.max') == 0, 'swap')
    require('memory.events\nlow 0\nhigh 0\nmax 0\noom 0\noom_kill 0\noom_group_kill 0\n' in text,
            'memory events')
    require('\tExit status: 0\n' in text, 'job failed')
    return value('memory.peak')


def audit(root):
    build = json.loads((root/'codegen/build.json').read_text())
    native = json.loads((root/'native/report.json').read_text())
    checks = json.loads((root/'checks-build/build.json').read_text())
    # Relocate only this experiment's recorded root; original oracle/archive
    # paths remain exact. This permits replay of a retained read-only copy.
    old_root = str(Path(next(iter(native['pins']))).parent.parent)
    def relocate(path):
        old = Path(path)
        return root/old.relative_to(old_root) if old.is_relative_to(old_root) else old
    for path, digest in native['pins'].items():
        require(sha(relocate(path)) == digest, 'native input pin')
    for name, digest in build['source_pins'].items():
        require(sha(root/'source'/name) == digest, 'codegen source pin')
    drivers = root/'checks-build/drivers'
    current = Path(__file__).resolve().parent
    for name in ('test_avx2_pair_schedule.cpp','test_gfni_boundary.cpp',
                 'gfni_boundary_screen.cpp','avx2_isa_screen.cpp'):
        require(sha(drivers/name) == sha(current/name), 'compiled driver source drift')
    result = {'bead':'leopard-79h.38.5.4.18.2', 'timed':False, 'loops':{}}
    for mode in range(3):
        directory = root/'codegen'/str(mode)
        obj = directory/'Leopard2BackendAVX2.cpp.o'
        require(sha(obj) == build['modes'][str(mode)]['object_sha256'], 'codegen object pin')
        fresh = subprocess.check_output(['objdump','-drwC',str(obj)], text=True)
        saved = (directory/'disassembly.txt').read_text()
        require(instructions(fresh) == instructions(saved), 'saved disassembly differs from actual object')
        for _, raw, assembly in instructions(fresh):
            require(raw[0] != 0x62 and not re.search(r'%ymm(?:1[6-9]|2[0-9]|3[01])\b|%zmm|\bvpternlog|\bvgf2p8', assembly),
                    'AVX2 ISA ceiling violated')
        loop = pair_loop(fresh)
        require(len(loop['stack_references']) == (0 if mode == 2 else 4), 'unexpected stack references')
        result['loops'][str(mode)] = loop
    require(build['modes']['0']['object_sha256'] ==
        'bf615cd0bb211ea3b8fd8a2c7e6c8a9651f5317857cd8ddcd910c03e74fe1e5d', 'production OFF identity')
    require(native['completed'] is True and native['timed'] is False, 'incomplete native checks')
    expected_labels = []
    for profile in ('release','sanitize'):
        labels = ['pairs','split',*map(str,range(8)),'roundtrip','concurrent',
                  *('parity-'+str(i) for i in range(4)),'clock-guard',
                  *('bad-'+str(i) for i in range(4))]
        expected_labels += [profile+'-'+label for label in labels]
        archive = root/'checks-build'/profile/'candidate.a'
        state = checks['profiles'][profile]
        require(sha(archive) == state['archive_sha256'], 'candidate archive drift')
        original = Path(state['original'])
        require(sha(original) == state['original_sha256'], 'original archive drift')
        original_members = subprocess.check_output(['ar','t',str(original)],text=True).splitlines()
        require(original_members == subprocess.check_output(['ar','t',str(archive)],text=True).splitlines(), 'archive inventory')
        require(len(original_members) == 24 and len(state['unchanged_members']) == 23, 'archive member coverage')
        candidate_member = subprocess.check_output(['ar','p',str(archive),'Leopard2BackendAVX2.cpp.o'])
        require(hashlib.sha256(candidate_member).hexdigest() == state['object_sha256'], 'candidate member identity')
        for member, digest in state['unchanged_members'].items():
            actual = subprocess.check_output(['ar','p',str(archive),member])
            before = subprocess.check_output(['ar','p',str(original),member])
            require(actual == before and hashlib.sha256(actual).hexdigest() == digest, 'unrelated member changed')
        prefix = root/'native'/profile
        for cell in range(8):
            guard = json.loads(Path(str(prefix)+'-'+str(cell)+'.stdout').read_text())
            require(guard['cell'] == cell and guard['subset_masks'] == 6 and
                    guard['timed'] is False and guard['field'] == (1 if cell == 6 else 2),
                    'public subset/field coverage')
        for cell in range(4):
            parity = json.loads(Path(str(prefix)+'-parity-'+str(cell)+'.stdout').read_text())
            require(parity['cell'] == cell and parity['encode_calls'] == 1 and
                    parity['samples_ns'] == [] and parity['outer_guards'] is True and
                    parity['codec_commit'] == 'pair-schedule-mode2', 'untimed parity record')
        require(Path(str(prefix)+'-pairs.stdout').read_text() == 'pair kernel cases: 66147\n', 'pair coverage')
        require(Path(str(prefix)+'-split.stdout').read_text() == 'split range cases: 256\n', 'range coverage')
        require(Path(str(prefix)+'-roundtrip.stdout').read_text() ==
                ''.join('roundtrip cell '+str(i)+' passed\n' for i in range(3)), 'roundtrip coverage')
        concurrent = Path(str(prefix)+'-concurrent.stdout').read_text().splitlines()
        require(concurrent.count('roundtrip cell 1 passed') == 9 and
                concurrent.count('roundtrip cell 2 passed') == 8 and
                concurrent.count('four-thread both-field roundtrips passed') == 1 and len(concurrent) == 18,
                'concurrent coverage')
    require([r['label'] for r in native['records']] == expected_labels, 'record inventory/order')
    for record in native['records']:
        label = record['label']
        expected = 86 if label.endswith('clock-guard') else 1 if '-bad-' in label else 0
        require(type(record['returncode']) is int and record['returncode'] == record['expected'] == expected, 'exit code')
        for ext in ('stdout','stderr'):
            path = root/'native'/(label+'.'+ext)
            require(sha(path) == record[ext+'_sha256'], 'raw output drift')
        stderr = (root/'native'/(label+'.stderr')).read_text()
        if expected == 0: require(not stderr, 'unexpected stderr')
        elif expected == 86: require(stderr == 'unexpected driver benchmark clock\n', 'clock guard')
    count = 0
    require(len(native['parity_comparisons']) == 8, 'parity inventory')
    for record in native['parity_comparisons']:
        actual, original = relocate(record['path']), Path(record['original'])
        require(sha(actual) == record['sha256'] == sha(original) == record['original_sha256'], 'parity hash')
        require(actual.stat().st_size == original.stat().st_size == record['bytes'], 'parity length')
        with actual.open('rb') as a, original.open('rb') as b:
            while True:
                chunk = a.read(65536)
                require(chunk == b.read(65536), 'full parity bytes')
                if not chunk: break
                count += len(chunk)
    require(count == native['parity_bytes'] == 69599232, 'parity coverage')
    result.update(native_records=42, positive_native_records=32, full_parity_bytes=count,
        memory_peaks={name:resource(root/name, limit) for name,limit in
            [('codegen-v2.log',512*1024**2),('checks-build.log',512*1024**2),('native.log',256*1024**2)]})
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('workspace', type=Path)
    args = parser.parse_args()
    print(json.dumps(audit(args.workspace.resolve(strict=True)),indent=2))

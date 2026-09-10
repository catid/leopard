#!/usr/bin/env python3
"""Untimed retained-input investigation; never execute a codec or select a winner.

Tracker: leopard-79h.38.5.4.19.1.4. Objdump text normalization is a diagnostic,
not a proof of semantic equivalence or of the cause of any timing difference.
"""
import argparse
import hashlib
from itertools import zip_longest
import json
from pathlib import Path
import re
import statistics
import subprocess

BEAD = 'leopard-79h.38.5.4.19.1.4'
PINS = {
    'tower/frozen/original': 'dd0c82b408e00744716f809c1c87a7f972eb5d7be010079b5154c44cbd290f68',
    'tower/frozen/current': '45ba62393647dcc64873eb592402e551b33de37719245540310ca0bf988920b9',
    'tower/frozen/original.a': '89f33d3d43f792de078c38469f9406dded2b69126af2c6fc8db6e2acb77cd334',
    'tower/frozen/current.a': '24a78109405a25a5164e656bcaabbcc9aff586849a7899492cbaf290b3018ba9',
    'tower/frozen/l2-driver.o': '99223acc12f2f616c0e4e5c0cb39d7a1a5207aa27570eedf3f08a1ec9290aa75',
    'tower/frozen/build.json': 'd304b72038eb460efd69ccd2526bd26c187fec083be7ba2f2724c7085a479818',
    'tower/frozen/leopard2.cpp': '93412028272454b2457f9928a6bc4c34684591f5bee2c882a3ad0000fd50cc70',
    'tower/frozen/avx2_adjacent_public.cpp': '3e098de4ac232fc58f2816908f09f7c882efd6b0f2f84605cec5a8835e47cccd',
    'tower/attempt1/attempt.json': '4b0040d9455698a22cbde47b0803cabfcd9c889a8695931671662d70475de163',
    'paired/frozen/current': '2a8cca6b4d4edfb945186506e3d170d868dc5ef9cbfef46dd1a3e9181b2d688c',
    'paired/frozen/current.a': '89f33d3d43f792de078c38469f9406dded2b69126af2c6fc8db6e2acb77cd334',
    'paired/frozen/paired_timer_r19932.cpp': '58890549f49632f657fd225b914f95cef15d7a88c6c3c940880db64587deb20c',
    'paired/attempt1/attempt.json': '40c2f5bcccfaf89efaaf83a939b2362e7e9c51af3ec6435207b0eb38bcd52409',
    'paired/SHA256SUMS': '2c1dd4857bc6dc35de7d6ce33421f6ef7de50e448f438d7743835b842df2d5c4',
}
SYMBOLS = (
    'main', 'leo2_encode',
    '_ZL14EncodeInternalPK10leo2_codecmPKPKvPKPvS6_mS3_mbb',
    '_ZN12_GLOBAL__N_1L24SelectTransformEncodeOpsEPK10leo2_codecmjjbb',
    '_ZN12_GLOBAL__N_1L26ExecuteTransformEncodePassEPK10leo2_codecRKN7leopard7backend3OpsEmmbjjPKPKvPPvSD_PKN17leopard2_internal26SparseForwardPlanBatchViewEbb',
    '_ZN7leopard3ff817ReedSolomonEncodeERKNS_7backend3OpsEmjjjjPKPKvPPvPKN17leopard2_internal26SparseForwardPlanBatchViewEbbb',
    '_ZN7leopard7backendL22AVX2FF8HighEncodeSmallEPKPKvjPKPvjPKhS9_m',
    '_ZN7leopard7backendL21AVX2FF8IFFTButterfly4EPvS1_S1_S1_tttm',
    '_ZN7leopard7backendL20AVX2FF8FFTButterfly4EPvS1_S1_S1_tttm',
    '_ZN12_GLOBAL__N_1L21ValidateEncodeBuffersEPK10leo2_codecmPKPKvPKPvS7_mS4_mPb',
    '_ZN7leopard7backendL13AVX2XorMemoryEPvPKvm',
    '_ZN7leopard7backendL20AVX2FF8FFTButterfly2EPvS1_tm',
    '_ZN7leopard7backendL21AVX2FF8IFFTButterfly2EPvS1_tm',
    '_ZN7leopard7backendL24AVX2FF8IFFTButterfly2XorEPKvS2_PvS3_tm',
    '_ZN7leopard7backendL24AVX2FF8IFFTButterfly4OutEPKvS2_S2_S2_PvS3_S3_S3_tttm',
    '_ZN7leopard7backendL25AVX2FF8FFTButterfly4RangeEPKPvjtttmb',
    '_ZN7leopard7backendL29AVX2FF8WeightedIFFTButterfly4EPKvS2_S2_S2_PvS3_S3_S3_tttthtttm',
)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(data):
    return hashlib.sha256(data).hexdigest()


def strict_json(data):
    def unique(pairs):
        result = {}
        for key, value in pairs:
            require(key not in result, 'duplicate JSON key')
            result[key] = value
        return result
    def bad_constant(value):
        raise ValueError('nonfinite JSON constant: ' + value)
    return json.loads(data, object_pairs_hook=unique, parse_constant=bad_constant)


def typed_equal(a, b):
    return json.dumps(a, sort_keys=True, allow_nan=False) == json.dumps(b, sort_keys=True, allow_nan=False)


def differences(a, b):
    return [dict(index=i, original=x, current=y)
            for i, (x, y) in enumerate(zip_longest(a, b)) if x != y]


def read(path):
    require(path.is_file() and path.stat().st_size <= 8 * 1024 * 1024, 'input bound')
    return path.read_bytes()


def members(data):
    """Bounded GNU ar member digests, including newly added tower objects."""
    require(data[:8] == b'!<arch>\n', 'regular ar required')
    offset, names, result = 8, b'', {}
    while offset < len(data):
        header = data[offset:offset + 60]
        require(len(header) == 60 and header[58:] == b'`\n', 'ar header')
        size = int(header[48:58].strip())
        require(0 <= size <= 4 * 1024 * 1024, 'ar size')
        start = offset + 60
        payload = data[start:start + size]
        require(len(payload) == size, 'ar payload')
        offset = start + size
        if size % 2:
            require(data[offset:offset + 1] == b'\n', 'ar padding')
            offset += 1
        name = header[:16].decode('ascii').strip()
        if name in ('/', '/SYM64/'):
            continue
        if name == '//':
            require(not names, 'duplicate name table')
            names = payload
            continue
        if name.startswith('/'):
            index = int(name[1:])
            require(0 <= index < len(names) and b'/\n' in names[index:], 'ar name offset')
            name = names[index:].split(b'/\n', 1)[0].decode('ascii')
        else:
            name = name.removesuffix('/')
        require(name and Path(name).name == name and name not in result, 'ar member name')
        result[name] = digest(payload)
    require(offset == len(data) and result, 'ar EOF')
    return result


def tool(argv):
    # Only these static inspection tools are permitted, never an input binary.
    require(argv[0] in ('/usr/bin/nm', '/usr/bin/objdump'), 'static tool only')
    run = subprocess.run(argv, capture_output=True, check=True, timeout=30,
                         env={'PATH': '/usr/bin:/bin', 'LC_ALL': 'C'})
    require(not run.stderr and len(run.stdout) < 2 * 1024 * 1024, 'tool output')
    return run.stdout.decode('ascii')


def symbols(path):
    result = {}
    for line in tool(['/usr/bin/nm', '-S', '--defined-only', str(path)]).splitlines():
        fields = line.split()
        if len(fields) == 4 and fields[2] in ('t', 'T'):
            address, size = int(fields[0], 16), int(fields[1], 16)
            result.setdefault(fields[3], []).append((address, size))
    for rows in result.values():
        rows.sort()
    return result


def normalize(assembly):
    """Remove ONLY symbol-annotated PC-relative address spelling.

    Keep literal constants, registers, widths, offsets and symbol names. RIP
    references without an objdump symbol annotation are deliberately unchanged.
    This does not normalize indirect targets or establish their runtime values.
    """
    if '#' in assembly:
        instruction, comment = assembly.split('#', 1)
        match = re.fullmatch(r'\s*[0-9a-f]+ <([^>]+)>\s*', comment)
        if match and re.search(r'-?0x[0-9a-f]+\(%rip\)', instruction):
            instruction = re.sub(r'-?0x[0-9a-f]+\(%rip\)',
                                 '<' + match[1] + '>(%rip)', instruction)
            return ' '.join(instruction.split())
    # Only direct PC-relative control transfers. Do not normalize an absolute
    # address or a comment on a non-RIP instruction merely because it has a name.
    transfer = re.fullmatch(r'(callq?|j[a-z]+|loop[a-z]*)\s+[0-9a-f]+ <([^>]+)>', assembly.strip())
    if transfer:
        assembly = transfer[1] + ' <' + transfer[2] + '>'
    return ' '.join(assembly.split())


def disassemble(path, address, size):
    raw = tool(['/usr/bin/objdump', '-d', '-w', '--show-raw-insn',
                '--start-address=' + str(address), '--stop-address=' + str(address + size), str(path)])
    record, instructions = parse_disassembly(raw, address, size)
    return raw, record, instructions


def parse_disassembly(raw, address, size):
    instructions, encoding, expected = [], bytearray(), address
    for line in raw.splitlines():
        match = re.fullmatch(r'\s*([0-9a-f]+):\s*\t([0-9a-f ]+)\t(.+)', line)
        if not match:
            continue
        pc = int(match[1], 16)
        require(pc == expected, 'disassembly discontinuity')
        data = bytes.fromhex(match[2])
        require(data, 'empty instruction')
        instructions.append((pc - address, normalize(match[3])))
        encoding.extend(data)
        expected += len(data)
    require(expected == address + size and instructions, 'complete symbol disassembly')
    record = dict(address=address, size=size, mod64=address % 64, mod4096=address % 4096,
                  instructions=len(instructions), raw_sha256=digest(encoding),
                  normalized_sha256=digest(json.dumps(instructions).encode()))
    return record, instructions


def grouped_loop(instructions, address):
    """Locate the public encode repetition backedge, without guessing alignment.

    Qualification's original/tower driver has two leo2_encode sites: preflight
    and grouped repetition. Only the latter is inside a backward main branch.
    """
    matches = []
    for i, (offset, opcode) in enumerate(instructions):
        if opcode != 'call <leo2_encode>':
            continue
        for j in range(i + 1, min(i + 10, len(instructions))):
            target = re.fullmatch(r'jne <main\+0x([0-9a-f]+)>', instructions[j][1])
            if target is None:
                continue
            start = int(target[1], 16)
            if not 0 < offset - start < 128:
                continue
            region = [(pc - start, op) for pc, op in instructions[:j + 1] if pc >= start]
            require(region and region[0][0] == 0, 'loop starts on instruction')
            matches.append(dict(start=address + start, start_mod64=(address + start) % 64,
                                call=address + offset, call_mod64=(address + offset) % 64,
                                backedge=address + instructions[j][0],
                                normalized_instructions=region))
    require(len(matches) == 1, 'one grouped public-call loop')
    return matches[0]


def process_summary(record):
    require(record['group'] == 256 and len(record['samples']) == 84, 'fixed GF8 group')
    samples = []
    for elapsed, per_call in record['samples']:
        require(type(elapsed) is int and elapsed > 0 and per_call == elapsed / 256, 'exact group average')
        samples.append(per_call)
    return dict(median=statistics.median(samples), minimum=min(samples), maximum=max(samples),
                slot_medians=[statistics.median(samples[slot::4]) for slot in range(4)],
                first_pass=statistics.median(samples[:4]), last_pass=statistics.median(samples[-4:]),
                retained_spans=len(samples), single_call_latency=False)


def audit(tower, paired, output):
    roots = {'tower': tower, 'paired': paired}
    paths = {key: roots[key.split('/')[0]] / key.split('/', 1)[1] for key in PINS}
    def verify():
        for key, path in paths.items():
            require(digest(read(path)) == PINS[key], 'pin: ' + key)
    verify()
    output.mkdir(mode=0o700, exist_ok=False)
    archives = {key: members(read(path)) for key, path in paths.items() if key.endswith('.a')}
    old = archives['tower/frozen/original.a']
    current = archives['tower/frozen/current.a']
    require(len(old) == 24 and len(current) == 26, 'complete archive inventory')
    require(archives['paired/frozen/current.a'] == old, 'paired/original archive identity')
    archive_comparison = dict(
        unchanged=[name for name in old if current.get(name) == old[name]],
        changed=[name for name in old if name in current and current[name] != old[name]],
        removed=sorted(set(old) - set(current)), added=sorted(set(current) - set(old)),
        hashes=archives)
    build = strict_json(read(paths['tower/frozen/build.json']))
    require(build['artifacts']['l2-release-objects/driver.o'] == PINS['tower/frozen/l2-driver.o'],
            'shared driver hash')
    links = [argv for argv in build['commands']
             if argv[-1].endswith(('/original/plain', '/release/plain'))]
    require(len(links) == 2, 'two original/tower links')
    drivers = [[arg for arg in argv if arg.endswith('/l2-release-objects/driver.o')] for argv in links]
    require(len(drivers[0]) == 1 and drivers[0] == drivers[1], 'same driver in both link records')
    functions, loops, original = [], {}, None
    for profile, key in (('original', 'tower/frozen/original'), ('tower', 'tower/frozen/current'),
                         ('paired', 'paired/frozen/current')):
        table = symbols(paths[key])
        profile_rows = {}
        for index, name in enumerate(SYMBOLS):
            require(name in table, 'missing symbol: ' + name)
            # Duplicate local symbols belong to separate backend variants.
            # Retain ALL occurrences; do not silently choose the first variant.
            for occurrence, (address, size) in enumerate(table[name]):
                raw, row, instructions = disassemble(paths[key], address, size)
                filename = f'{profile}-{index}-{occurrence}.asm'
                (output / filename).write_text(raw)
                row.update(profile=profile, symbol=name, occurrence=occurrence,
                           disassembly=filename, disassembly_sha256=digest(raw.encode()))
                identity = (name, occurrence)
                if original is not None:
                    baseline = original[identity]
                    row['normalized_equal_original'] = instructions == baseline
                    if instructions != baseline:
                        changes = differences(baseline, instructions)
                        row['difference_count'] = len(changes)
                        row['instruction_count_delta'] = len(instructions) - len(baseline)
                        row['first_differences'] = changes[:8]
                profile_rows[identity] = instructions
                if name == 'main' and profile in ('original', 'tower'):
                    loops[profile] = grouped_loop(instructions, address)
                functions.append(row)
        if original is None:
            original = profile_rows
        else:
            require(set(profile_rows) == set(original), 'symbol occurrence inventory')
    # Diagnostic projection of ALL GF8 same-OFF processes, no filtering or new
    # performance decision. Read retained samples; never execute the codec.
    manifest = {}
    for line in read(paths['paired/SHA256SUMS']).decode('ascii').splitlines():
        value, name = line.split('  ', 1)
        require(name not in manifest and re.fullmatch('[0-9a-f]{64}', value), 'manifest entry')
        manifest[name] = value
    attempt = strict_json(read(paths['paired/attempt1/attempt.json']))
    require(attempt['complete'] is True and len(attempt['invocations']) == 318, 'complete paired attempt')
    gf8 = [row for row in attempt['invocations'] if row['cell'] == 8 and row['comparison'] == 'same_off']
    require([(r['round'], r['slot']) for r in gf8] == [(r, s) for r in range(3) for s in range(4)],
            'all twelve same-OFF processes in original order')
    samples = []
    for row in gf8:
        filename = f"cell-8-round-{row['round']}-same_off-slot-{row['slot']}.stdout"
        raw = read(paired / 'attempt1' / filename)
        require(digest(raw) == manifest['attempt1/' + filename], 'raw/manifest identity')
        require(typed_equal(strict_json(raw), row['record']), 'typed raw/journal equality')
        require(row['order'] == '0000' and row['sibling_delta'] == 0, 'same-OFF identity')
        samples.append(dict(round=row['round'], process_slot=row['slot'], raw_sha256=digest(raw),
                            **process_summary(row['record'])))
    verify()
    result = dict(schema='leopard-gf8-retained-shift-audit/v1', bead=BEAD,
                  codec_executions=0, new_timings=0, performance_decision=None,
                  pins=PINS, archive_comparison=archive_comparison, functions=functions,
                  original_tower_link_records=links, grouped_public_loop=loops,
                  grouped_loop_normalized_equal=(loops['original']['normalized_instructions'] ==
                                                 loops['tower']['normalized_instructions']),
                  requested_backend_from_pinned_source={'paired_GF8_cell8': 'AUTO',
                                                        'tower_GF8_cell7': 'AVX2'},
                  paired_gf8_same_off=samples,
                  limitation='Static placement and normalized disassembly are not a causal performance attribution. '
                  'Historical process mappings, allocation addresses and hardware counters are absent.')
    (output / 'result.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps(dict(changed=archive_comparison['changed'], added=archive_comparison['added'],
                         functions=len(functions), differing=[(r['profile'], r['symbol'], r['occurrence'])
                         for r in functions if r.get('normalized_equal_original') is False]), indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('tower', type=Path)
    parser.add_argument('paired', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    audit(args.tower, args.paired, args.output)

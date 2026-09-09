#!/usr/bin/env python3
"""Untimed retained full-butterfly evidence replay; leopard-79h.18.20.3.

Uses stdlib and objdump only; does not execute a codec or import its collector.
Checks native records and actual object code, not independent native execution.
"""
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def manifest(root, name):
    entries = {}
    for line in (root / name).read_text().splitlines():
        match = re.fullmatch(r'([0-9a-f]{64})  ([^\n]+)', line)
        require(match is not None, 'manifest syntax')
        digest, relative = match.groups()
        path = Path(relative)
        require(not path.is_absolute() and '..' not in path.parts and relative not in entries,
                'manifest path or duplicate')
        require(sha(root / path) == digest, 'manifest drift: ' + relative)
        entries[relative] = digest
    require(entries, 'empty manifest')
    return entries


def stack(row):
    _, mnemonic, operands = row
    return bool(re.search(r'\[[^\]]*\b(?:rsp|rbp|esp|ebp)\b', operands) or
                re.fullmatch(r'(?:push|pop)[a-z]*|(?:l?call)[a-z]*|enter|leave', mnemonic))


def codegen(assembly):
    functions = {}
    current = None
    for line in assembly.splitlines():
        head = re.fullmatch(r'[0-9a-f]+ <(.+)>:', line)
        if head:
            current = head[1]
            require(current not in functions, 'duplicate function')
            functions[current] = []
            continue
        match = re.match(r'\s*([0-9a-f]+):\s+((?:[0-9a-f]{2} )+)\s*(\S+)(.*)', line)
        if match:
            require(current is not None, 'instruction outside function')
            address, raw, mnemonic, operands = match.groups()
            require(raw.split()[0] != '62', 'EVEX instruction')
            require(not re.search(r'\b(?:[xyz]mm(?:1[6-9]|[2-9][0-9])|zmm\d+|k[0-7])\b', operands),
                    'higher ISA register')
            require(not re.search(r'gf2|ternlog|pclmul', mnemonic), 'excluded ISA')
            functions[current].append((int(address, 16), mnemonic, operands))
    required = {'tower_ifft_pair': 6, 'tower_ifft_accumulate': 6,
                'tower_convert_involution': 2, 'forward_shared': 6}
    result = {}
    for name, rows in functions.items():
        key = 'forward_shared' if '::pair<false, false, false>(' in name else name
        for end, mnemonic, operands in rows:
            target = re.match(r'\s*([0-9a-f]+) <', operands)
            if mnemonic not in ('ja', 'jne') or not target or int(target[1], 16) >= end:
                continue
            start = int(target[1], 16)
            loop = [row for row in rows if start <= row[0] <= end]
            shuffles = sum(row[1] == 'vpshufb' for row in loop)
            if not shuffles:
                continue
            require(key in required and key not in result, 'unexpected/duplicate shuffle loop')
            require(shuffles == required[key], 'shuffle count')
            require(not any(stack(row) for row in loop), 'vector loop stack access')
            result[key] = dict(loop_start=hex(start), loop_end=hex(end), instructions=len(loop),
                               shuffles=shuffles, masks=sum(row[1] == 'vpand' for row in loop),
                               shifts=sum(row[1] == 'vpsrlw' for row in loop),
                               xors=sum(row[1] == 'vpxor' for row in loop),
                               vector_loop_stack_accesses=0,
                               whole_function_stack_accesses=sum(stack(row) for row in rows))
    require(set(result) == set(required), 'missing shuffle loop')
    # Both ordinary forward entrypoints must actually jump to the counted helper.
    for name in ('tower_fft_pair', 'tower_fft_out'):
        require(name in functions and any(row[1] == 'jmp' and '::pair<false, false, false>(' in row[2]
                                         for row in functions[name]), 'forward helper dispatch')
    return result


def conversion_proof(basis):
    require(len(basis) == 16, 'basis length')
    def poly(a, b):
        result = 0
        while b:
            if b & 1:
                result ^= a
            a <<= 1
            if a & 65536:
                a ^= 0x1002D
            b >>= 1
        return result
    polynomial = [0] * 65536
    for x in range(1, 65536):
        bit = x & -x
        polynomial[x] = polynomial[x ^ bit] ^ basis[bit.bit_length() - 1]
    require(len(set(polynomial)) == 65536, 'noninvertible basis')
    canonical = {p: x for x, p in enumerate(polynomial)}
    products = [canonical[poly(polynomial[256], polynomial[b])] for b in range(256)]
    require(all(x >> 8 == b for b, x in enumerate(products)), 'high map not identity')
    require(all(products[b] == products[b & 15] ^ products[b & 240] for b in range(256)),
            'nibble map decomposition')
    def convert(x):
        return x ^ (products[x >> 8] & 255)
    require(all(convert(convert(x)) == x for x in range(65536)), 'not involutive')
    return dict(high_identity_values=256, involution_symbols=65536, table_bytes=32,
                u_times_basis=[products[1 << i] for i in range(8)])


def native_record(data):
    expected = dict(tracker='leopard-79h.18.20.3', exhaustive_butterfly_cases=65536*4,
                    boundary_butterfly_cases=14*3*7*4 + 32*3*4,
                    input_pair_basis_symbols_per_log=32, logs_including_zero_skew=65536,
                    conversion_symbols=65536, high_identity_values=256, zero_byte_null_checks=9,
                    timed=False, public_codec_qualified=False)
    require(data == expected, 'native record counts/claims')
    return expected


def resource(path, cap):
    lines = path.read_text().splitlines()
    require(lines.count('memory.peak') == 1 and lines.count('\tExit status: 0') == 1, 'resource status')
    index = lines.index('memory.peak')
    peak = int(lines[index + 1])
    require(0 < peak <= cap, 'resource cap')
    require(lines[index+2:] == ['memory.max', str(cap), 'memory.events', 'low 0', 'high 0', 'max 0',
                               'oom 0', 'oom_kill 0', 'oom_group_kill 0', 'memory.swap.current', '0',
                               'memory.swap.max', '0'], 'resource events')
    return dict(peak=peak, maximum=cap, all_events_zero=True, swap=0, exit_status=0)


def replay(root):
    inputs = manifest(root, 'inputs.sha256')
    artifacts = manifest(root, 'artifacts.sha256')
    expected = {'release/codec.a': '89f33d3d43f792de078c38469f9406dded2b69126af2c6fc8db6e2acb77cd334',
                'sanitize/codec.a': '501d5b0375455b137c8c24befb2e152a0ea2d773f8b094aa1388d0b5e6cb4410'}
    for name, digest in expected.items():
        require(inputs[name] == artifacts[name] == digest, 'oracle archive identity')
    require(inputs['probe.o'] == artifacts['probe.o'], 'codegen/native object identity')
    source = (root / 'source/test_tower_algebra.cpp').read_text()
    basis = [int(x, 16) for x in re.findall(r'0x[0-9a-fA-F]+', re.search(r'basis\[16\]\s*=\s*\{([^}]+)\}', source)[1])]
    records = []
    for profile in ('release', 'sanitize'):
        require((root / (profile + '.status')).read_text() == '0\n', 'native exit')
        require((root / (profile + '.stderr')).read_bytes() == b'', 'native stderr')
        records.append(native_record(json.loads((root / (profile + '.stdout')).read_text())))
        for argument in ('--measure', '--help', '1'):
            stem = profile + '-cli-' + argument
            require((root / (stem + '.status')).read_text() == '2\n', 'CLI exit')
            require((root / (stem + '.stdout')).read_bytes() == b'', 'CLI stdout')
            require((root / (stem + '.stderr')).read_text() == 'No options or timing mode supported\n', 'CLI stderr')
    assembly = (root / 'codegen.txt').read_text()
    actual = subprocess.check_output(['objdump', '-drwC', '-Mintel', str(root / 'probe.o')], text=True, timeout=10)
    normalize = lambda value: re.sub(r'^.*: +file format elf64-x86-64$', 'OBJECT: file format elf64-x86-64', value, flags=re.M)
    require(normalize(actual) == normalize(assembly), 'actual disassembly drift')
    return dict(tracker='leopard-79h.18.20.3', native_records=records, cli_rejections=6,
                conversion=conversion_proof(basis), codegen=codegen(assembly),
                resources={name: resource(root / name, cap) for name, cap in
                           [('codegen-build.log', 536870912), ('native-build.log', 536870912),
                            ('native-check.log', 268435456)]},
                object_sha256=artifacts['probe.o'], codec_archives=expected,
                timed=False, codec_integrated=False, independent_log_derivation=False)


if __name__ == '__main__':
    require(len(sys.argv) == 2, 'usage: replay_tower_butterflies.py ROOT')
    print(json.dumps(replay(Path(sys.argv[1]).resolve()), sort_keys=True, indent=2))

#!/usr/bin/env python3
"""Collector-free untimed algebra/object replay, leopard-79h.18.20.2.

No native execution or production-module imports. Proof by linearity checks
all 65,536 fixed multipliers on all 16 canonical input basis vectors.
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


def field_proof(basis):
    # Shift/reduce multiplication: independent of C++ full convolution.
    def poly(a, b):
        value = 0
        while b:
            if b & 1:
                value ^= a
            a <<= 1
            if a & 65536:
                a ^= 0x1002D
            b >>= 1
        return value
    polynomial = [0] * 65536
    canonical = [0] * 65536
    for x in range(1, 65536):
        bit = x & -x
        polynomial[x] = polynomial[x ^ bit] ^ basis[bit.bit_length() - 1]
    require(len(set(polynomial)) == 65536, 'basis not invertible')
    for x, p in enumerate(polynomial):
        canonical[p] = x
    def product(a, b):
        return canonical[poly(polynomial[a], polynomial[b])]
    subfield = [product(a, b) for a in range(256) for b in range(256)]
    require(max(subfield) < 256, 'subfield closure')
    for b in range(256):
        linear = [0] * 256
        for a in range(1, 256):
            bit = a & -a
            linear[a] = linear[a ^ bit] ^ subfield[bit * 256 + b]
            require(linear[a] == subfield[a * 256 + b], 'subfield linearity')
    delta = product(256, 256) ^ 256
    require(delta == 128, 'quadratic relation')
    require(all((subfield[a * 256 + a] ^ a) != delta for a in range(256)), 'quadratic reducible')
    u_times = [product(256, b) for b in range(256)]
    require(len({x >> 8 for x in u_times}) == 256, 'high map not invertible')
    high_inverse = [0] * 256
    for b, value in enumerate(u_times):
        high_inverse[value >> 8] = b
    def to_tower(x):
        b = high_inverse[x >> 8]
        return ((x & 255) ^ (u_times[b] & 255)) | (b << 8)
    def from_tower(x):
        return (x & 255) ^ u_times[x >> 8]
    tower = [to_tower(x) for x in range(65536)]
    for x in range(65536):
        require(from_tower(tower[x]) == x and to_tower(from_tower(x)) == x, 'conversion roundtrip')
        if x:
            bit = x & -x
            require(tower[x] == tower[x ^ bit] ^ tower[bit], 'conversion linearity')
    for input_bit in range(16):
        a = tower[1 << input_bit] & 255
        b = tower[1 << input_bit] >> 8
        basis_product = [product(1 << input_bit, 1 << i) for i in range(16)]
        expected = [0] * 65536
        for coefficient in range(65536):
            if coefficient:
                bit = coefficient & -coefficient
                expected[coefficient] = expected[coefficient ^ bit] ^ basis_product[bit.bit_length() - 1]
            c = tower[coefficient] & 255
            d = tower[coefficient] >> 8
            ac = subfield[a * 256 + c]
            bd = subfield[b * 256 + subfield[delta * 256 + d]]
            cross = subfield[(a ^ b) * 256 + (c ^ d)]
            observed = from_tower((ac ^ bd) | ((cross ^ ac) << 8))
            require(observed == expected[coefficient], 'three-product identity')
    return dict(delta=delta, conversion_symbols=65536, subfield_products=65536,
                constant_basis_products=1048576, u_times_basis=[u_times[1 << i] for i in range(8)])


def codegen(text):
    functions = {}
    current = None
    for line in text.splitlines():
        head = re.fullmatch(r'[0-9a-f]+ <(tower_\w+)>:', line)
        if head:
            current = head[1]
            functions[current] = []
        match = re.match(r'\s*([0-9a-f]+):\s+((?:[0-9a-f]{2} )+)\s*(\S+)(.*)', line)
        if match and current:
            address, raw, mnemonic, operands = match.groups()
            require(raw.split()[0] != '62', 'EVEX instruction')
            require(not re.search(r'\b(?:[xyz]mm(?:1[6-9]|[2-9][0-9])|zmm\d+|k[0-7])\b', operands), 'higher ISA register')
            require(not re.search(r'gf2|ternlog|pclmul', mnemonic), 'excluded ISA')
            functions[current].append((int(address, 16), mnemonic, operands))
    require(set(functions) == {'tower_product_blocks', 'tower_convert_blocks'}, 'function set')
    result = {}
    for name, rows in functions.items():
        branches = [(address, int(m[1], 16)) for address, mnemonic, operands in rows
                    if mnemonic == 'jne' and (m := re.match(r'\s*([0-9a-f]+) <', operands))
                    and int(m[1], 16) < address]
        require(len(branches) == 1, 'unique loop')
        end, start = branches[0]
        loop = [row for row in rows if start <= row[0] <= end]
        stack = [row for row in rows if re.search(r'\[(?:r|e)(?:sp|bp)\b', row[2])]
        shuffles = sum(row[1] == 'vpshufb' for row in loop)
        require(not stack, 'stack references')
        require(shuffles == (6 if name == 'tower_product_blocks' else 4), 'shuffle count')
        result[name] = dict(loop_start=hex(start), loop_end=hex(end), instructions=len(loop),
                            shuffles=shuffles, shifts=sum(row[1] == 'vpsrlw' for row in loop),
                            masks=sum(row[1] == 'vpand' for row in loop),
                            xors=sum(row[1] == 'vpxor' for row in loop), stack_references=0)
    return result


def resource(path, maximum, expected_exit):
    lines = path.read_text().splitlines()
    require(lines.count('memory.peak') == 1, 'resource peak')
    require(lines.count('\tExit status: ' + str(expected_exit)) == 1, 'resource exit')
    i = lines.index('memory.peak')
    require(0 < int(lines[i + 1]) <= maximum, 'memory cap')
    require(lines[i + 2:] == ['memory.max', str(maximum), 'memory.events', 'low 0', 'high 0', 'max 0',
                            'oom 0', 'oom_kill 0', 'oom_group_kill 0',
                            'memory.swap.current', '0', 'memory.swap.max', '0'], 'resource counters')
    return dict(peak=int(lines[i + 1]), maximum=maximum, exit_status=expected_exit, all_events_zero=True, swap=0)


def replay(root):
    build = json.loads((root / 'build.json').read_text())
    checks = json.loads((root / 'checks.json').read_text())
    for name, digest in build['sources'].items():
        require(sha(root / 'source' / name) == digest, 'source drift')
    for name, digest in build['artifacts'].items():
        require(sha(root / name) == digest, 'artifact drift')
    require(checks['artifact_sha256'] == build['artifacts'], 'execution input mismatch')
    require(build['artifacts']['release/codec.a'] == '89f33d3d43f792de078c38469f9406dded2b69126af2c6fc8db6e2acb77cd334', 'release archive')
    require(build['artifacts']['sanitize/codec.a'] == '501d5b0375455b137c8c24befb2e152a0ea2d773f8b094aa1388d0b5e6cb4410', 'sanitizer archive')
    source = (root / 'source/LeopardFF16.cpp').read_text()
    basis = [int(x, 16) for x in re.findall(r'0x[0-9a-fA-F]+', re.search(r'kCantorBasis\[kBits\]\s*=\s*\{([^}]+)\}', source)[1])]
    require(len(basis) == 16, 'basis length')
    proof = field_proof(basis)
    require(len(checks['results']) == 2 and checks['cli_rejections'] == 6, 'native record count')
    for index, profile in enumerate(('release', 'sanitize')):
        observed = json.loads((root / (profile + '.stdout')).read_text())
        require(checks['results'][index] == dict(profile=profile, result=observed, exit_status=0), 'native record')
        require((root / (profile + '.stderr')).read_bytes() == b'', 'native stderr')
        expected = dict(proof, tracker='leopard-79h.18.20.2', additional_products=1048576,
                        vector_product_checks=2097152, mutated_table_rejections=6, polynomial=65581, u=256,
                        table_bytes_per_multiplier=96, conversion_table_bytes_each=64,
                        timed=False, production_integrated=False)
        require(observed == expected, 'native/algebra mismatch')
        for i, argument in enumerate(('--measure', '--help', '1')):
            require((root / f'{profile}-cli-{i}.stdout').read_bytes() == b'', 'CLI stdout')
            require((root / f'{profile}-cli-{i}.stderr').read_text() == 'No timing/CLI options supported: ' + argument + '\n', 'CLI rejection')
    assembly = (root / 'release/codegen.txt').read_text()
    actual_assembly = subprocess.check_output(['objdump', '-drwC', '-Mintel', str(root / 'release/tower_avx2_probe.o')], text=True)
    # Only objdump's filename banner changes in a read-only relocated copy.
    normalize = lambda value: re.sub(r'^.*: +file format elf64-x86-64$', 'OBJECT: file format elf64-x86-64', value, flags=re.M)
    require(normalize(actual_assembly) == normalize(assembly), 'actual disassembly differs')
    return dict(tracker='leopard-79h.18.20.2', proof=proof, codegen=codegen(assembly),
                resources={name: resource(root / path, cap, status) for name, path, cap, status in
                           [('initial_build', 'build.log', 536870912, 1),
                            ('accepted_build', 'finish-build.log', 536870912, 0),
                            ('native_checks', 'check.log', 268435456, 0)]},
                release_object_sha256=sha(root / 'release/tower_avx2_probe.o'), timed=False,
                codec_integrated=False, independent_model_review=False)


if __name__ == '__main__':
    require(len(sys.argv) == 2, 'usage: replay_tower_algebra.py ROOT')
    print(json.dumps(replay(Path(sys.argv[1]).resolve()), indent=2, sort_keys=True))

#!/usr/bin/env python3
"""Actual-object, untimed audit; leopard-79h.38.5.4.18.3.

This is codegen evidence only, not correctness, dynamic attribution or timing.
"""
import argparse
import json
from pathlib import Path
import re
import subprocess
from audit_avx2_pair_schedule import FUNCTION, instructions, require, resource, sha

FUNCTIONS = {
    'inverse': FUNCTION,
    'forward': 'void leopard::backend::AVX2FF16Butterfly2<false>(void*, void*, unsigned short, unsigned long)',
    'forward_range': 'void leopard::backend::AVX2FF16Butterfly2RangePrepared<false>(void* const*, unsigned int, unsigned int, unsigned int, unsigned short, unsigned long)',
    'accumulating': 'leopard::backend::AVX2FF16IFFTButterfly2Xor(void const*, void const*, void*, void*, unsigned short, unsigned long)',
}


def selected_loop(text, function):
    parts = re.split(r'^([0-9a-f]+) <(.+)>:\n', text, flags=re.M)
    bodies = [parts[i + 2] for i in range(1, len(parts), 3) if parts[i + 1] == function]
    require(len(bodies) == 1, 'missing or duplicate selected function')
    ops = instructions(bodies[0])
    candidates = []
    for end, _, assembly in ops:
        branch = re.match(r'j(?!mp\b)[a-z]+\s+([0-9a-f]+)\s', assembly)
        if not branch or int(branch[1], 16) >= end:
            continue
        start = int(branch[1], 16)
        body = [op for op in ops if start <= op[0] <= end]
        if sum(op[2].split()[0] == 'vpshufb' for op in body) == 8:
            require(body[0][0] == start, 'branch target is not an instruction')
            candidates.append((start, end, body))
    # Prepared ranges nest the byte loop inside a shard-pair loop. Only the
    # innermost eight-shuffle loop denotes one 64-byte iteration. Two distinct
    # innermost loops are ambiguous and remain a hard failure.
    innermost = [(start, end, body) for start, end, body in candidates
                 if not any(start <= a and b <= end and (start, end) != (a, b)
                            for a, b, _ in candidates)]
    require(len(innermost) == 1, 'ambiguous or missing innermost eight-shuffle loop')
    start, end, body = innermost[0]
    return dict(start=hex(start), end=hex(end), instructions=len(body), byte_shuffles=8,
                stack_references=[op[2] for op in body if re.search(r'\([^)]*%(?:rsp|rbp)\b', op[2])])


def isa_ceiling(text):
    for _, raw, assembly in instructions(text):
        require(raw[0] != 0x62 and not re.search(
            r'%ymm(?:1[6-9]|2[0-9]|3[01])\b|%zmm|\bvpternlog|\bvgf2p8', assembly),
            'AVX2 ISA ceiling violated')


def audit(root):
    build = json.loads((root / 'codegen/build.json').read_text())
    require(build['bead'] == 'leopard-79h.38.5.4.18.3' and build['timed'] is False, 'build scope')
    require(set(build['modes']) == {'0', '1', '2', '3'}, 'mode inventory')
    require({p.name: sha(p) for p in (root / 'source').iterdir()} == build['source_pins'], 'source pins')
    result = dict(bead=build['bead'], timed=False, codec_executed=False,
                  correctness_qualified=False, dynamic_path_counts_qualified=False, modes={})
    for mode in range(4):
        directory = root / 'codegen' / str(mode)
        obj = directory / 'Leopard2BackendAVX2.cpp.o'
        require(sha(obj) == build['modes'][str(mode)]['object_sha256'], 'object drift')
        actual = subprocess.check_output(['objdump', '-drwC', str(obj)], text=True)
        require(instructions(actual) == instructions((directory / 'disassembly.txt').read_text()), 'disassembly drift')
        isa_ceiling(actual)
        result['modes'][str(mode)] = dict(object_sha256=sha(obj),
            loops={key: selected_loop(actual, name) for key, name in FUNCTIONS.items()})
    require(result['modes']['0']['object_sha256'] ==
            'bf615cd0bb211ea3b8fd8a2c7e6c8a9651f5317857cd8ddcd910c03e74fe1e5d', 'production baseline')
    result['build_memory_peak'] = resource(root / 'codegen.log', 512 * 1024**2)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('workspace', type=Path)
    print(json.dumps(audit(parser.parse_args().workspace.resolve(strict=True)), indent=2))

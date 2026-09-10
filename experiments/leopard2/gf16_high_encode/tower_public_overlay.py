"""Pinned adaptation of the qualified public driver; no historical source edits.

Only the link interface, diagnostics and schema change. Public calls, workload,
clock boundary, parity comparisons and grouped schedule remain identical.
Tracker: leopard-79h.38.5.4.18.4.2. This is not timing preregistration.
"""
import hashlib

BASES = {
    'driver': '3e098de4ac232fc58f2816908f09f7c882efd6b0f2f84605cec5a8835e47cccd',
    'witness': '7c5e563c397cce8c12ab36a0e8458c2e4a6c0b66722edba069e0b0d7b6b2f10a',
}


def adapt(source, kind):
    if kind not in BASES or hashlib.sha256(source.encode()).hexdigest() != BASES[kind]:
        raise ValueError('unqualified public frontend source')
    source = source.replace('avx2_adjacent_public_link.h', 'tower_public_link.h')
    source = source.replace('LeoAdjacent', 'LeoTower')
    if kind == 'witness':
        return source
    source = source.replace('leopard-79h.38.5.4.18.3.1', 'leopard-79h.38.5.4.18.4.2')
    source = source.replace('adjacent-scheduling', 'tower-encoder')
    source = source.replace('select adjacent state', 'select tower state')
    source = source.replace('leopard-adjacent-public/v1', 'leopard-tower-public/v1')
    begin = source.index('void PrintCounts(')
    end = source.index('\n}\n}', begin) + 2
    source = source[:begin] + '''void PrintCounts(const LeoTowerCounts& counts)
{
    std::printf("{\\\"values\\\":[");
    for (unsigned i = 0; i < 8; ++i)
        std::printf("%s%llu", i ? "," : "",
            static_cast<unsigned long long>(counts.values[i]));
    std::printf("],\\\"initializations\\\":%u}", counts.initializations);
}''' + source[end:]
    source = source.replace('pair_probes', 'tower_probes').replace('pair_totals', 'tower_totals')
    source = source.replace('LeoTowerCounts tower_probes[3]', 'LeoTowerCounts tower_probes[4]')
    source = source.replace('''            // First public call may initialize caches. Only the next three are
            // isolated pair probes; initialization is explicitly excluded.''',
        '''            // Each preflight records tower work and lifetime cache initialization.
            // Cold setup is outside later warm-throughput spans, never hidden as free.''')
    source = source.replace('if (slot) tower_probes[slot - 1] = LeoTowerPublicCounts();',
                            'tower_probes[slot] = LeoTowerPublicCounts();')
    source = source.replace('''        Require(!tower_totals.overflow, "pair counter overflow");
        for (unsigned i = 0; i < 3; ++i) Require(!tower_probes[i].overflow, "probe overflow");
''', '')
    source = source.replace('i < 3; ++i) { if (i) std::printf(","); PrintCounts(tower_probes[i]);',
                            'i < 4; ++i) { if (i) std::printf(","); PrintCounts(tower_probes[i]);')
    # A fresh process starts without the tower cache, including untraced Release.
    anchor = '        Require(LeoTowerPublicDefault(), "runtime must default OFF");'
    if source.count(anchor) != 1:
        raise ValueError('initialization anchor')
    source = source.replace(anchor, anchor + '''
        Require(LeoTowerPublicCounts().initializations == 0, "tower cache must start cold");''')
    if any(word in source for word in ('LeoAdjacent', 'pair_probes', '.overflow', 'pair_totals')):
        raise ValueError('incomplete tower adaptation')
    return source

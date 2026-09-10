"""Separate timing-capable diagnostic adapter; no clocks or codec execution.

Only pinned, completed three-epoch sources are accepted. Qualification and
preregistration are still required before the actual steady binary is timed.
"""
import hashlib

BEAD = 'leopard-79h.38.5.4.19.1.4.4'
DRIVER_SHA = '9cb3c3a1d550d6f06f551f5aff5d31e2f5365afa8a44bbc535300440e4facd16'
CLOCK_SHA = '746905101d048ae2023e75205a7a1325514a19320af57d194664d41201d447b8'
SCHEMA = 'leopard-paired-epoch-diagnostic/v1'
PLAIN_MARK = '''// No public-call observer; epoch marks are outside sample boundaries.
extern "C" void LeoPairedWitnessMark(unsigned) {}
'''


def pinned(source, digest):
    if hashlib.sha256(source.encode()).hexdigest() != digest:
        raise ValueError('unqualified timing adapter input')
    return source


def once(source, old, new):
    if source.count(old) != 1:
        raise ValueError('timing adapter anchor not unique: ' + old[:80])
    return source.replace(old, new)


def driver(source):
    original = pinned(source, DRIVER_SHA)
    source = once(source,
        '        Require(argc < 2 || std::strcmp(argv[1], "--measure"), "metadata refuses real timing");\n', '')
    source = once(source,
        '''        Require(!std::strcmp(clock_kind, "abort") || !std::strcmp(clock_kind, "synthetic"),
            "metadata requires an abort or synthetic clock binding");''',
        '''        Require(!std::strcmp(clock_kind, "abort") || !std::strcmp(clock_kind, "synthetic") ||
                !std::strcmp(clock_kind, "steady"), "unknown diagnostic clock binding");''')
    source = once(source, 'leopard-paired-epoch-r19932/v1', SCHEMA)
    for part in original.split('        const auto encode = [&]() {')[1:]:
        body = '        const auto encode = [&]() {' + part.split('\n        };', 1)[0] + '\n        };'
        if source.count(body) != 1:
            raise ValueError('public encode lambda changed')
    body = original.split('                const bool sampled =', 1)[1].split('                check_buffers();', 1)[0]
    if source.count('                const bool sampled =' + body) != 1:
        raise ValueError('group timing boundary changed')
    return source


def fake_steady_clock(source):
    # This test artifact deliberately labels itself steady so --measure takes
    # its actual branch. Only its link/hash proves that no real clock is read.
    return once(pinned(source, CLOCK_SHA),
        'extern "C" const char* LeoPairedClockKind() { return "synthetic"; }',
        'extern "C" const char* LeoPairedClockKind() { return "steady"; }')

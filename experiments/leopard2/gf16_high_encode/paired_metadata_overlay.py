"""Untimed-only adaptation; leopard-79h.38.5.4.19.1.4.1.

The consumed timer source remains untouched. Refuse any other source version.
"""
import hashlib

BASE_SHA = '58890549f49632f657fd225b914f95cef15d7a88c6c3c940880db64587deb20c'
BEAD = 'leopard-79h.38.5.4.19.1.4.1'


def adapt(source):
    if hashlib.sha256(source.encode()).hexdigest() != BASE_SHA:
        raise ValueError('unqualified paired timer source')
    original = source

    def replace(old, new):
        nonlocal source
        if source.count(old) != 1:
            raise ValueError('overlay anchor not unique: ' + old[:60])
        source = source.replace(old, new)

    replace('int main(int argc, char** argv)',
            '#include "PairedRuntimeMetadata.h"\n\nint main(int argc, char** argv)')
    replace('    try {\n', '''    try {
        Require(argc < 2 || std::strcmp(argv[1], "--measure"), "metadata refuses real timing");
''')
    replace('        const char* const clock_kind = LeoPairedClockKind();', '''        const char* const clock_kind = LeoPairedClockKind();
        Require(!std::strcmp(clock_kind, "abort") || !std::strcmp(clock_kind, "synthetic"),
            "metadata requires an abort or synthetic clock binding");''')
    replace('        const auto select = [&](unsigned) { ++selections; };', '''        const auto select = [&](unsigned slot) {
            paired_metadata::Select(selections, slot, -1, -1, -1, -1);
            ++selections;
        };''')
    replace('''            Require(diag::AutoGF16GFNIEncodeSelectedForDiagnostics(codec.get(), cell.bytes) == selected_gfni(slot),
                "selected route");''', '''            const bool selected = diag::AutoGF16GFNIEncodeSelectedForDiagnostics(codec.get(), cell.bytes);
            Require(selected == selected_gfni(slot), "selected route");
            paired_metadata::Select(selections, slot, static_cast<int>(options.backend),
                static_cast<int>(leo2_context_backend(context.get())),
                diag::AutoGF16GFNIR19932EnabledForDiagnostics() ? 1 : 0, selected ? 1 : 0);''')
    replace('''        for (unsigned slot = 0; slot < 4; ++slot) {
            select(slot);''', '''        const auto snapshot = [&](unsigned endpoint) {
            paired_metadata::Capture(endpoint, source, reference, scratch,
#ifdef LEO_PAIRED_NATIVE
                NULL,
#else
                &output,
#endif
                parity, output_bytes, inputs, outputs, cell.bytes);
        };
        snapshot(0);
        for (unsigned slot = 0; slot < 4; ++slot) {
            select(slot);''')
    replace('        const unsigned expected = exercise ? 1 + 25 * group : 1;',
            '        snapshot(1);\n        const unsigned expected = exercise ? 1 + 25 * group : 1;')
    replace('        return 0;\n', '        paired_metadata::Print(selections, probes);\n        return 0;\n')
    # Explicitly guard the complete public call lambdas and RunGroup boundary.
    for block in original.split('        const auto encode = [&]() {')[1:]:
        body = '        const auto encode = [&]() {' + block.split('\n        };', 1)[0] + '\n        };'
        if source.count(body) != 1:
            raise ValueError('public call boundary changed')
    span = original.split('                const bool sampled =', 1)[1].split('                check_buffers();', 1)[0]
    if source.count('                const bool sampled =' + span) != 1:
        raise ValueError('group boundary changed')
    return source

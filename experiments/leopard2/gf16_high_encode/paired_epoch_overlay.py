"""Separate clock-free three-epoch adapters; never edit consumed frontend sources."""
import hashlib
from paired_metadata_overlay import adapt as metadata_adapt

BEAD = 'leopard-79h.38.5.4.19.1.4.3'
PINS = {
    'PairedRuntimeMetadata.h': '20f9325ab7238a39e06810181e05f59cf6b646135c7b7a39cc2f089efd5e6da4',
    'paired_timer_clock.cpp': 'e4fa103e6d943444629bd794e7d7166a66c44b010d0676b3212b73081f422105',
    'paired_timer_witness.cpp': '36229212de2775bd160959f8baa44cdd8972e858f3d7eb48dff69ab2c662425a',
}

def pinned(name, source):
    if hashlib.sha256(source.encode()).hexdigest() != PINS[name]:
        raise ValueError('unqualified epoch input: ' + name)
    return source

def once(source, old, new):
    if source.count(old) != 1:
        raise ValueError('epoch anchor not unique: ' + old[:80])
    return source.replace(old, new)

def driver(source):
    original = source
    source = metadata_adapt(source)
    source = once(source, 'extern "C" const char* LeoPairedClockKind();',
                  'extern "C" const char* LeoPairedClockKind();\nextern "C" void LeoPairedWitnessMark(unsigned);')
    source = once(source, 'samples.reserve(84);', 'samples.reserve(252);')
    source = once(source, '        snapshot(0);', '''        unsigned epoch_probes[3][4] = {};
        unsigned epoch_calls[3] = {}, epoch_selections[3] = {};
        unsigned epoch_slots[3][4] = {};
        for (unsigned epoch = 0; epoch < 3; ++epoch) {
        const unsigned calls_before = calls, selections_before = selections;
        unsigned slots_before[4];
        for (unsigned slot = 0; slot < 4; ++slot) slots_before[slot] = per_slot[slot];
        LeoPairedWitnessMark(2 * epoch);
        snapshot(2 * epoch);''')
    source = once(source, 'if (slot == 0) std::memcpy(reference.data, parity, output_bytes);',
                  'if (epoch == 0 && slot == 0) std::memcpy(reference.data, parity, output_bytes);')
    source = once(source, '        snapshot(1);', '''        epoch_calls[epoch] = calls - calls_before;
        epoch_selections[epoch] = selections - selections_before;
        for (unsigned slot = 0; slot < 4; ++slot) {
            epoch_slots[epoch][slot] = per_slot[slot] - slots_before[slot];
            epoch_probes[epoch][slot] = probes[slot];
            Require(epoch_slots[epoch][slot] == (exercise ? 1 + 25 * group : 1), "epoch slot calls");
        }
        Require(epoch_calls[epoch] == 4 * (exercise ? 1 + 25 * group : 1) &&
            epoch_selections[epoch] == (exercise ? 104U : 4U), "epoch accounting");
        Require(Hash(source.data, input_bytes) == input_hash, "epoch input changed");
#ifndef LEO_PAIRED_NATIVE
        // Check before the next preflight can reset the real probe counter.
        paired_metadata::RequireQuiescentProbe(probes[3]);
#endif
        snapshot(2 * epoch + 1);
        LeoPairedWitnessMark(2 * epoch + 1);
        }
''')
    source = once(source, 'const unsigned expected = exercise ? 1 + 25 * group : 1;',
                  'const unsigned expected = 3 * (exercise ? 1 + 25 * group : 1);')
    source = once(source, 'selections == (exercise ? 104U : 4U)',
                  'selections == (exercise ? 312U : 12U)')
    source = once(source, 'samples.size() == (clocks_enabled ? 84U : 0U)',
                  'samples.size() == (clocks_enabled ? 252U : 0U)')
    source = source.replace('leopard-paired-timer-r19932/v1', 'leopard-paired-epoch-r19932/v1')
    source = source.replace('\\"warmup_passes\\"', '\\"warmup_passes_per_epoch\\"')
    source = source.replace('\\"exercise_passes\\"', '\\"exercise_passes_per_epoch\\"')
    source = once(source, '        paired_metadata::Print(selections, probes);', '''        std::printf("{\\"schema\\":\\"paired-epoch-accounting/v1\\",\\"epochs\\":[");
        for (unsigned epoch = 0; epoch < 3; ++epoch) {
            std::printf("%s{\\"epoch\\":%u,\\"calls\\":%u,\\"selections\\":%u,"
                "\\"per_slot_calls\\":[%u,%u,%u,%u],\\"probes\\":[%u,%u,%u,%u],"
                "\\"sample_begin\\":%u,\\"sample_count\\":%u}", epoch ? "," : "", epoch,
                epoch_calls[epoch], epoch_selections[epoch],
                epoch_slots[epoch][0], epoch_slots[epoch][1], epoch_slots[epoch][2], epoch_slots[epoch][3],
                epoch_probes[epoch][0], epoch_probes[epoch][1], epoch_probes[epoch][2], epoch_probes[epoch][3],
                clocks_enabled ? epoch * 84 : 0, clocks_enabled ? 84 : 0);
        }
        std::puts("],\\"timed\\":false}");
        paired_metadata::Print(selections, epoch_probes);''')
    # Complete encode lambdas and the grouped boundary must remain source-identical.
    for part in original.split('        const auto encode = [&]() {')[1:]:
        block = '        const auto encode = [&]() {' + part.split('\n        };',1)[0] + '\n        };'
        if source.count(block) != 1: raise ValueError('epoch public lambda changed')
    block = original.split('                const bool sampled =',1)[1].split('                check_buffers();',1)[0]
    if source.count('                const bool sampled =' + block) != 1:
        raise ValueError('epoch grouped boundary changed')
    return source

def header(source):
    source = pinned('PairedRuntimeMetadata.h', source)
    source = source.replace('leopard-79h.38.5.4.19.1.4.1', BEAD)
    source = once(source, 'kSelections = 104', 'kSelections = 312')
    source = once(source, 'Snapshot snapshots[2];', 'Snapshot snapshots[6];')
    source = once(source, 'static Store records = {};', '''static Store records = {};
// This destructor also runs on the clock wrapper's std::exit path. It reports
// only progress, never a complete metadata/success record. Store is trivial
// static storage and remains readable throughout process teardown.
struct Progress {
    ~Progress() {
        std::printf("{\\\"schema\\\":\\\"paired-epoch-progress/v1\\\",\\\"selection_count\\\":%u,"
            "\\\"snapshot_count\\\":%u,\\\"timed\\\":false}\\n", records.selection_count, records.snapshot_count);
    }
};
static Progress progress;''')
    source = once(source, 'inline uintptr_t End(uintptr_t address, uint64_t bytes)', '''#ifndef LEO_PAIRED_NATIVE
inline void RequireQuiescentProbe(unsigned expected)
{
    namespace d = leopard2_internal;
    // The mode getter maps armed raw3 and normalized raw1 to the same value.
    // Finish returns false without writes only for an already-normalized mode.
    // If it normalizes a leaked raw3 probe, fail immediately; never continue.
    Require(d::AutoGF16GFNIEncodeCallCountForDiagnostics() == expected &&
        d::AutoGF16GFNIEncodeModeForDiagnostics() == 1 &&
        !d::FinishAutoGF16GFNIEncodeRouteProbeForDiagnostics(), "epoch probe leaked into exercise");
}
#endif

inline uintptr_t End(uintptr_t address, uint64_t bytes)''')
    source = once(source, 'endpoint < 2,', 'endpoint < 6,')
    source = once(source, 'if (endpoint == 1) Require(Same(records.snapshots[0], s)',
                  'if (endpoint != 0) Require(Same(records.snapshots[0], s)')
    source = once(source, 'const unsigned* probes)', 'const unsigned probes[3][4])')
    source = once(source, 'records.snapshot_count == 2 && records.selection_count == selections',
                  'records.snapshot_count == 6 && records.selection_count == selections && (selections == 12 || selections == 312)')
    source = source.replace('paired-runtime-metadata/v1', 'paired-epoch-runtime-metadata/v1')
    source = once(source, '\\"observation\\":\\"new_frontend_endpoints\\"',
                  '\\"observation\\":\\"new_three_epoch_endpoints\\"')
    source = once(source, '''        "\\"preflight_gfni_counts\\":[%u,%u,%u,%u],\\"selections\\":[", probes[0],probes[1],probes[2],probes[3]);''',
        '''        "\\"preflight_gfni_counts\\":[[%u,%u,%u,%u],[%u,%u,%u,%u],[%u,%u,%u,%u]],\\"selections\\":[",
        probes[0][0],probes[0][1],probes[0][2],probes[0][3],probes[1][0],probes[1][1],probes[1][2],probes[1][3],
        probes[2][0],probes[2][1],probes[2][2],probes[2][3]);''')
    source = once(source, '        const Selection& s = records.selections[i];',
                  '        const Selection& s = records.selections[i];\n        const unsigned stride = selections / 3, local = i % stride;')
    source = once(source, '\\"index\\":%u,\\"phase\\":', '\\"index\\":%u,\\"epoch\\":%u,\\"phase\\":')
    source = once(source, 'i ? "," : "", i, i < 4 ? "preflight" : "exercise", i < 4 ? -1 : static_cast<int>((i-4)/4),',
                  'i ? "," : "", i, i / stride, local < 4 ? "preflight" : "exercise", local < 4 ? -1 : static_cast<int>((local-4)/4),')
    source = once(source, 'for (unsigned i = 0; i < 2; ++i)', 'for (unsigned i = 0; i < 6; ++i)')
    source = once(source, '''std::printf("%s{\\"endpoint\\":%u,\\"allocations\\":{", i ? "," : "", i);''',
                  '''std::printf("%s{\\"endpoint\\":%u,\\"epoch\\":%u,\\"phase\\":\\"%s\\",\\"allocations\\":{",
            i ? "," : "", i, i / 2, i % 2 ? "after" : "before");''')
    return source

def clock(source):
    source = pinned('paired_timer_clock.cpp', source)
    source = source.replace('168', '504').replace('paired-synthetic-clock/v1','paired-epoch-synthetic-clock/v1')
    source = once(source, '    if (fault && *fault) {', '''    const char* fault_epoch = std::getenv("LEO_PAIRED_TEST_CLOCK_FAULT_EPOCH");
    if (fault_epoch && (std::strlen(fault_epoch) != 1 || *fault_epoch < '0' || *fault_epoch > '2')) {
        std::fputs("invalid synthetic fault epoch\\n", stderr); std::exit(89);
    }
    const unsigned target_epoch = fault_epoch ? static_cast<unsigned>(*fault_epoch - '0') : 0;
    if (fault && *fault && index / 168 == target_epoch) {''')
    return source

def witness(source):
    source = pinned('paired_timer_witness.cpp', source)
    return source + '''
// Fixed-capacity boundary marks; the public witness itself is never reset.
namespace {
struct EpochMarks {
    struct Mark { unsigned calls, states[3], apis[3]; uint64_t order; };
    Mark marks[6] = {};
    unsigned count = 0;
    ~EpochMarks() {
        std::printf("{\\"schema\\":\\"paired-epoch-public-marks/v1\\",\\"marks\\":[");
        for (unsigned i = 0; i < count; ++i) {
            const Mark& m = marks[i];
            std::printf("%s{\\"endpoint\\":%u,\\"calls\\":%u,\\"states\\":[%u,%u,%u],"
                "\\"apis\\":[%u,%u,%u],\\"order_hash\\":\\"%016llx\\"}", i ? "," : "", i,
                m.calls,m.states[0],m.states[1],m.states[2],m.apis[0],m.apis[1],m.apis[2],
                static_cast<unsigned long long>(m.order));
        }
        std::puts("],\\"timed\\":false}");
    }
} epoch_marks;
}
extern "C" void LeoPairedWitnessMark(unsigned endpoint)
{
    if (endpoint != epoch_marks.count || endpoint >= 6) {
        std::fputs("epoch witness mark order/capacity\\n", stderr); std::exit(87);
    }
    EpochMarks::Mark& mark = epoch_marks.marks[endpoint];
    mark.calls = witness.calls; mark.order = witness.order;
    for (unsigned i = 0; i < 3; ++i) { mark.states[i] = witness.states[i]; mark.apis[i] = witness.apis[i]; }
    ++epoch_marks.count;
}
'''

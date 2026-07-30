# Leopard2 C API

`leopard2.h` adds an object-based API without changing `leopard.h` or the old
`leo_*` behavior.  The old entry points retain their count and 64-byte-size rules.

## Object lifetimes

Create one `leo2_context`, then any number of immutable codecs.  A codec fixes
`K`, `R`, profile, field, shard layout, parent length, and coordinate map.  A
decode plan copies one presence/erasure pattern and precomputes its locator
values, profile-specific normalization factors, and deterministic
received-coordinate selection.  A context must outlive every codec created from
it, and a codec must outlive every plan created from it.

Codec and plan execution is read-only and may be called concurrently when every
call has distinct outputs and scratch.  Setup may allocate; encode and reusable
plan execution do not.

Setup functions clear a non-null output handle before validating the remaining
arguments, so a failed create never leaves a stale object pointer.  Scratch-size
queries similarly clear a non-null size output to zero before validation, after
first proving that the complete scalar output object span is representable.
The wire-size query follows the same rule.  Null introspection arguments return
their documented zero, AUTO, or native-layout sentinel.

`leo2_codec_create` accepts profile and field selectors through integer ABI
parameters and validates them before converting to the internal enum types.
Applications still pass the `leo2_profile` and `leo2_field` constants normally.
This representation keeps negative and unknown values well-defined under C
short-enum modes and C++ enum sanitizers instead of invoking undefined behavior
while trying to reject them.  Profile, field, layout, and backend introspection
remain enum-typed because the library returns only values it has validated.

Option structures must be zero-initialized, `struct_size` set, and every
`reserved` field left zero.  A nonzero reserved field is rejected.  A
`struct_size` larger than the current definition is accepted and only the known
prefix is read; unknown trailing bytes are ignored.  The `backend` and
`shard_layout` option members use fixed-width `uint32_t` ABI storage even when a
C caller enables short enums.  Batch-item structures are fixed V1 layouts; a
future incompatible extension will use new suffixed types and entry points.

## Encoding

Call `leo2_encode_scratch_size`, allocate at least that many bytes aligned to
`leo2_scratch_alignment()`, then call `leo2_encode`.  Every original pointer is
required.  A null recovery output requests that parity shard be skipped.  Inputs
may alias other inputs, but outputs and scratch must be disjoint from all shard
ranges; unsupported overlap returns `LEO2_OVERLAP`.

The high profile calls the existing block-IFFT accumulation and truncated final
FFT and matches old Leopard in the compatibility test.  The low profile performs
one padded-`P` IFFT followed by only the shifted parity blocks needed by the
requested output prefix/subset.

A ragged transform encode runs its aligned prefix directly from application
buffers, then reuses the same `2P` or `2T` work slots for one staged 64-byte
tail.  Source staging is fixed at `64*K` bytes; Low V1 additionally reserves
`64*R` bytes for compact parity-tail scatter.  These fixed terms do not grow
with shard length.  GF16 uses the same compact/ALTMAP tail conversion as the
single-pass encoder, and requested null parity outputs remain unmaterialized.

For bounded `K,R <= 16` codecs, setup can precompute exact full-parent
systematic generator rows. AUTO uses this allocation-free direct evaluator only
for the measured Low V1, one-output, regular-shard region; scalar requires
`K >= 3`, SSSE3/AVX2/explicit AVX512 requires `K >= 2`, and physical shards
must be at least 1 KiB and a multiple of 64 bytes. All other cells retain the transform encoder.
This is an internal kernel decision and does not change profile identity or
parity bytes. See `leopard2_direct_encode.md` for the derivation and evidence.

## Decoding

Create a plan from byte presence arrays with `leo2_decode_plan_create`.  Execution
must provide pointers that match that pattern.  Missing original inputs are null,
and their corresponding `restored_original` entries are writable.  Present
original output entries are ignored.  Missing parity is never rebuilt implicitly.

`leo2_decode_plan_execute` restores missing originals only.  A pattern with no
missing originals has zero scratch and is a true no-op, including when parity is
missing.  `R=1` high-profile repair and `K=1` low-profile repair use direct paths.
Other patterns use the specialized IT2026 low- or high-rate decoder selected by
the codec profile.  Locator construction, profile normalization, pruning inputs,
and output selection depend only on the erasure pattern and are performed during
plan creation; execution is the byte-heavy reusable step.

The plan stores a compact fused-transform output bitmap, generic and per-block
input prefixes, sorted requested coordinates, and sparse high-rate output-block
descriptors. It neither rebuilds the former field-sized error bitfield nor
scans shard pointers for transform prefixes during execution. See
`leopard2_decode_plan_schedules.md` for representation, bounds, differential
tests, and pinned measurements.

For a physical shard size that is a multiple of 64 bytes, selected transform
inputs already use the kernels' native tile layout and remain read-only.
Execution therefore references those caller shards directly.  The specialized
low decoder uses at most `min(N,2P)` full-shard data slots: one accumulator and
one reusable parent-block tile when that is smaller.  The specialized high
decoder similarly uses `min(N,2T+L)` slots, where `L` is the number of requested
missing originals; its last `L` tiled slots retain outputs until
application-layout scatter.  The retained materialized specialized kernel is
used when its regular `N` slots are no larger.  One offline-calibrated
legacy-high GF8 region also reserves and uses `N` slots because that regular
traversal is faster for a single stripe; measured AVX2 and explicit AVX512
multi-item batches execute tiled using the same conservative allocation.  SSSE3
batches retain the materialized traversal pending dedicated evidence.  Forced generic decoding
retains `N` work slots.  A ragged final
tile adds `64*(K+R)` bytes for GF8 zero padding or GF16 compact-to-ALTMAP
scatter.  The public overlap checks are unchanged.

`leo2_decode_plan_scratch_size` uses the exact `L` for its immutable pattern.
The pattern-independent one-shot query uses `min(N,2T+R)` for a high profile,
so it is safe for every valid pattern accepted later by `leo2_decode`.  A
ragged shard is executed as an aligned prefix plus one padded 64-byte tail.
Its `K+R` staging term is therefore fixed at 64 bytes per public coordinate;
it no longer scales with the full shard length.  Pointer and range-validation
metadata are additional small terms; neither is a full-shard data slot.  See
`leopard2_decode_tiled_workspace.md` for the execution equivalence and scratch
formulas.  In the measured GF8 generic-fallback region, aligned SSSE3/AVX2/AVX512
prefixes of at least 4 KiB also reveal directly from the `N`-slot workspace into
caller outputs, eliminating the separate in-scratch reveal and public scatter
passes without changing scratch size or wire identity.

The plan uses exactly `K` received public coordinates.  It keeps every surviving
systematic shard, then keeps the lowest-index received parity shards needed to
reach `K`; any surplus received parity is treated as a deterministic virtual
erasure.  This selection affects only the work schedule, not the decoded message.
Applications still pass the original presence pattern and may leave pointers for
surplus received parity populated.

For differential testing and diagnosis, set at most one of
`LEO2_CODEC_FORCE_GENERIC_DECODE` or
`LEO2_CODEC_FORCE_SPECIALIZED_DECODE` in `leo2_codec_options.flags`.  The first
selects the retained full `O(N log N)` active-parent decoder; the second selects
the profile-specific transform decoder and bypasses bounded direct-matrix repair
and automatic crossover dispatch.  The flags govern nontrivial transform or
direct-matrix repairs; no-loss, single-XOR-parity, and single-original copy
fast paths remain direct.  Supplying both flags is invalid.  These diagnostic
choices are not wire profiles and do not change encoded data.

Within a profile-specific decoder,
`LEO2_CODEC_FORCE_TILED_DECODE` and
`LEO2_CODEC_FORCE_MATERIALIZED_DECODE` select the side-sized or retained
`N`-slot traversal for offline comparison.  They imply specialized decoding,
may be combined with `LEO2_CODEC_FORCE_SPECIALIZED_DECODE`, and are mutually
exclusive.  Either workspace flag combined with forced generic is invalid.
Normal applications should leave all four flags clear and use deterministic
AUTO dispatch.  See `leopard2_decode_tiled_workspace.md` for the calibrated
region and its authenticated evidence.

### Dispatch inputs and plan reuse

The V1 API deliberately has no caller-supplied expected-plan-reuse hint and no
persistent batch-count hint.  A batch call already supplies its exact
`item_count`; current byte-kernel policy distinguishes a single item from a
multi-item batch only where isolated evidence supports that stable class.  Each
decode item continues to use the storage returned by
`leo2_decode_plan_scratch_size`.  A batch-specific path must fit that
conservative single-item allocation, so callers never need a second per-item
scratch query merely because the next call is a larger batch.

Plan construction is complete before byte-heavy execution.  The same immutable
plan may safely be shared by callers whose actual reuse differs, so a mutable or
guessed reuse count would make otherwise identical calls depend on application
forecasting without changing the underlying arithmetic.  Benchmarks therefore
report setup separately and amortize it at multiple reuse counts, while normal
execution selects a deterministic kernel and never benchmarks online.  A
future measured setup strategy that genuinely needs a reuse promise must use a
new versioned plan-options/create entry point; it must not reinterpret a
reserved field or alter an existing plan after construction.

This decision changes neither ABI nor defaults.  Backend and execution-layout
choices do not alter the selected profile, field, coordinates, parity ordering
or wire bytes.  Exact batch count remains available internally to the
scheduler, and a future internal exact-count rule is compatible only if it
preserves the documented scratch bound and deterministic results.

`leo2_decode` is a convenience wrapper that allocates and destroys a plan.  Use a
reusable plan when setup amortization matters.

## Profiles

Legacy high V1:

    T = ceil_pow2(R)
    N = ceil_pow2(K + T)
    D = N - T
    [ parity 0..T-1 ][ message T..T+K-1 ][ shortened zero tail ]

Only parity `0..R-1` is transmitted; `R..T-1` is punctured.

Low V1:

    P = ceil_pow2(K)
    N = ceil_pow2(P + R)
    [ message 0..K-1 ][ shortened zeros to P ]
    [ parity P..P+R-1 ][ punctured parity tail ]

`AUTO` deterministically uses high when `R <= K`, otherwise low.  Field auto uses
GF8 for a parent of at most 256 coordinates and GF16 otherwise.  Profile and field
are mathematical code identity; backend selection never changes them.

`LEO2_PROFILE_EXACT_EXPERIMENTAL_V1` reserves a distinct research code identity.
The production constructor currently returns `LEO2_UNSUPPORTED` for it; it never
falls back to or aliases either production V1 profile.

Shard layout is also code identity.  `LEO2_SHARD_LAYOUT_NATIVE_V1` is zero, so a
zero-initialized options structure and the exact API-version-1 prefix retain
native behavior.  API version 2 appends `shard_layout` to
`leo2_codec_options`.  Its storage is fixed-width `uint32_t` even when a C
caller uses short-enum compiler options; set it to a `leo2_shard_layout` value.
An older `struct_size` ending at that field's offset is accepted, while a size
containing only part of the field is rejected.
`LEO2_SHARD_LAYOUT_GF16_PADDED_ODD_V1` must be requested explicitly together
with `LEO2_FIELD_GF16`; field AUTO never chooses it.

## Byte lengths, physical storage, and batch calls

GF8 supports every positive byte length and internally handles SIMD tails through
caller scratch.  Native GF16 supports every positive even byte length.  Complete
64-byte tiles remain unchanged and use zero-copy encode input and decode input
paths.  A partial final GF16 tile contains `q` complete symbols in `2q`
application bytes:
the first `q` bytes are scattered to ALTMAP low lanes `0..q-1`, the next `q`
bytes to high lanes `0..q-1`, and unused lanes are zero.  Outputs use the inverse
gather.  Odd GF16 byte lengths return `LEO2_UNSUPPORTED` because an unpaired byte
is not a complete GF16 symbol.  A no-loss decode plan remains a true no-op with
zero scratch even for an otherwise unsupported byte length.  See
`leopard2_gf16_tails.md` for the construction and the odd-length limitation.

Core encode/decode `shard_bytes` always names the physical buffer size.  For an
explicit padded-odd codec and an odd application payload `B`, call
`leo2_codec_wire_shard_bytes` to obtain `W=B+1`, then use
`leo2_pack_systematic_shard` to form every physical systematic coordinate as
`payload || 0`.  Encode and decode use `W`.  Every parity coordinate stores all
`W` bytes: parity byte `B` may be nonzero and cannot be omitted.  Recovered
systematic wires can be converted back with `leo2_unpack_systematic_shard`,
which verifies their final zero before returning the first `B` bytes.

Pack and unpack allocate nothing and have `memmove` overlap semantics, including
in-place operation.  Their wire-size argument must exactly match the size query.
An unrepresentable payload or wire address span returns
`LEO2_INVALID_ARGUMENT` before either range is accessed.  Padded codecs reject
even application-size queries, odd physical core lengths, and nonzero final
bytes in systematic inputs.  Native layout remains the only layout for even
payloads, so all legacy-valid and existing compact-even parity bytes remain
unchanged.  The no-loss decoder intentionally performs no buffer or pad
validation because its contract is a true no-op.

The batch entry points execute independent items using each item's own buffers and
scratch.  `item_count` may not exceed `UINT32_MAX`; larger counts return
`LEO2_INVALID_ARGUMENT` before any item is read.  A context with an effective
`thread_count` greater than one owns a persistent pool control object, but starts
no worker threads during context construction.  Empty and single-item batches
remain serial and do not start workers.  The first larger batch starts
up to `min(thread_count - 1, item_count - 1)` workers.  A later, larger batch
can grow the pool to the additional parallelism it needs; started workers remain
persistent and are reused.  A lazy start can return `LEO2_OUT_OF_MEMORY` and a
failed pool remains failed deterministically.
One-time field/table initialization and immutable codec/decode-plan setup are
deliberately serial even in an OpenMP build, so `leo_init`, context creation,
GF16 normalization setup, and sparse or dense locator-plan construction do not
leave an OpenMP worker team behind.  This does not disable OpenMP in the
existing byte-heavy encode/decode execution paths.
The calling thread participates.  A count of one always executes serially.

On Linux, `thread_count = 0` counts the CPUs in `sched_getaffinity`, so cpuset
and container restrictions are honored.  On Windows it counts the process
affinity mask.  Other platforms, or a failed affinity query, use
`std::thread::hardware_concurrency()`.  The default falls back to one and is
capped at 128.  An explicit count from 1 through 128 is retained even when it
exceeds the current affinity budget; a larger count returns
`LEO2_INVALID_ARGUMENT`.
`leo2_context_thread_count` reports the effective value.  This option controls
parallelism across batch items, not the wire profile or the result bytes.

On Linux, a context that explicitly requests AVX2 also snapshots the creating
thread's affinity for GF16 cache-tiling policy.  Leopard2 reads the level-three
data cache for every allowed CPU and uses the smallest capacity, so an affinity
spanning asymmetric cache dies remains conservative while a context created
after pinning to a larger-cache die can defer tiling until the larger crossover.
AUTO, GFNI, and other backends retain the established 32-MiB calibration until
independently qualified.  Missing, malformed, or changing topology does the
same.  This snapshot affects scratch and execution geometry only, never field
arithmetic or wire identity.  Assuming every later execution CPU supports the
context's selected ISA, changing or widening affinity after context creation
remains correct but may make the cached performance policy suboptimal;
applications managing affinity should set it before creating the context and
retain the same mask for ordinary execution and every batch call.  The pool
starts lazily and may grow on a later batch, so newly created workers inherit
the affinity of that batch caller; Leopard2 does not pin them itself.

An explicit context `backend` option selects one immutable execution table.
`AUTO` reports the production-default table that passed the startup capability
checks and known-answer tests as its immutable baseline.  In one ordinary
production x86 binary, explicit `SCALAR`,
`SSSE3`, `AVX2`, and `AVX512` requests select that table when it
was compiled and supported by the
host.  Lower tables are allocated and known-answer-tested once, on the first
explicit context request, so legacy and `AUTO`-only applications do not pay
their setup or memory cost.  Qualification is serialized and its result is
cached: an unavailable ISA returns `LEO2_UNSUPPORTED`, allocation failure
returns `LEO2_OUT_OF_MEMORY`, and a known-answer-test failure returns
`LEO2_INTERNAL_ERROR`.  A context never changes the process-default table
used by the legacy `leo_*` API, and independent contexts selecting different
tables can execute concurrently.
`AUTO` deliberately reports AVX2 as its baseline on qualifying x86 hosts.
Isolated complete-codec evidence supports AVX-512VL for bounded GF16
legacy-high encode regions, but did not justify selecting it universally.
An AUTO codec may therefore widen one full-output transform encode to a
startup-qualified immutable AVX-512 table only on the calibrated AMD family
1Ah/model 44h host class and in a deterministic offline-calibrated region.
That region is legacy-high GF16 with `K >= 8`, `N >= 16`,
`2 <= R <= 4096`, all parity outputs requested, and a shard length that is an
exact 64-byte multiple from 64 bytes through 4 MiB inclusive.  Unknown CPU
models and cells outside those bounds retain AVX2.  Explicit backend requests
never widen, and a
failed or unavailable optional qualification retains the baseline.  The
reported context backend remains the baseline because the kernel choice is a
property of the codec, byte length, and requested outputs rather than a new
context or wire identity.  See
`docs/leopard2_avx512_codec.md`.  The explicit AVX-512 request
requires AVX2 plus AVX-512F/BW/VL and operating-system support for opmask and
ZMM state; its current kernels retain 256-bit data operations and use the
expanded architectural register file.  Unsupported hosts return
`LEO2_UNSUPPORTED` without executing the backend.

Backend choice affects execution-path thresholds and kernels only, never the
profile, field, coordinate map, or parity bytes.  Diagnostic builds still cap
`AUTO` and explicit selection at their forced backend.
Legacy custom builds that raise the entire core translation unit to SSSE3 or
AVX2 retain `AUTO` compatibility but reject explicit table requests, because
their inline kernels can bypass an exact lower-table selection.  The default
CMake build uses isolated ISA objects and supports the explicit x86 choices.

Byte-heavy calls select their immutable table once and pass it explicitly
through every transform and OpenMP worker; selection is not recovered from
mutable global or thread-local state.  Context-pool workers likewise receive
it with each batch item, so independent contexts cannot contaminate one
another.  On ARM, explicit `NEON`
remains an exact request for the existing native-NEON or SSE2NEON transform
path.  Selecting a lower scalar table on an active NEON path is rejected until
the native-NEON backend is fully represented by the same ops-table boundary;
introspection continues to describe effective execution rather than a scalar
tail or fallback implementation detail.

Only one batch call at a time uses a given context pool.  Concurrent batch calls
sharing that context are serialized at the scheduler, while calls using separate
contexts may run independently.  Ordinary encode or plan-execute calls remain
safe to invoke concurrently with distinct outputs and scratch.  Within a batch,
the implementation reports the result from the lowest-index failing item
deterministically; completion order is otherwise unspecified.

API version 6 provides an opt-in scalable alias preflight for repeated and large
batches. Query caller-owned storage with
`leo2_encode_batch_preflight_scratch_size` or
`leo2_decode_plan_batch_preflight_scratch_size`, then pass that aligned storage
to the corresponding `*_with_preflight_scratch` entry point. General codecs
select this path from two items; one item and legacy-high K=1,R=1 through eight
items retain their measured specialized compatibility paths. A flattened,
allocation-free interval sort and sweep replaces the
compatibility path's pairwise item comparisons without weakening atomic
rejection or immutable-input sharing.

When every buffer and per-item scratch address remains stable across many
encodes, `leo2_encode_batch_binding_create` deep-copies and validates that
metadata once. `leo2_encode_batch_binding_execute` then performs only the
captured byte-heavy work without allocation or repeated alias sorting. Source
bytes may change, but the codec/context and every captured buffer must outlive
the binding; one binding is not concurrently executable because it owns no
separate parity or scratch storage. The full algorithms, alias treatment,
tests, lifetime rules, and diagnostic crossover evidence are in
`docs/leopard2_batch_preflight.md`.

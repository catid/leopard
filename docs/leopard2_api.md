# Leopard2 C API

`leopard2.h` adds an object-based API without changing `leopard.h` or the old
`leo_*` behavior.  The old entry points retain their count and 64-byte-size rules.

## Object lifetimes

Create one `leo2_context`, then any number of immutable codecs.  A codec fixes
`K`, `R`, profile, field, parent length, and coordinate map.  A decode plan copies
one presence/erasure pattern and precomputes its locator values.  The context and
codec must outlive their codecs and plans respectively.

Codec and plan execution is read-only and may be called concurrently when every
call has distinct outputs and scratch.  Setup may allocate; encode and reusable
plan execution do not.

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

## Decoding

Create a plan from byte presence arrays with `leo2_decode_plan_create`.  Execution
must provide pointers that match that pattern.  Missing original inputs are null,
and their corresponding `restored_original` entries are writable.  Present
original output entries are ignored.  Missing parity is never rebuilt implicitly.

`leo2_decode_plan_execute` restores missing originals only.  A pattern with no
missing originals has zero scratch and is a true no-op, including when parity is
missing.  `R=1` high-profile repair and `K=1` low-profile repair use direct paths.
Other patterns currently use the generic active-parent LCH derivative decoder with
locator setup cached in the plan.  Specialized IT2026 low/high execution remains
behind its proof and differential gates.

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

## Byte lengths and batch calls

GF8 supports every positive byte length and internally handles SIMD tails through
caller scratch.  Legacy GF16 ALTMAP currently requires a multiple of 64 bytes;
partial tiles return `LEO2_UNSUPPORTED` because truncating their high halves is
not decodable.  See `leopard2_math_and_sources.md` for the proof obligation.

The initial batch entry points execute independent items using per-item scratch.
They provide a stable interface for later pool/executor and SIMD-across-stripes
dispatch; the current implementation preserves deterministic serial ordering.

# GF16 tower arithmetic: untimed algebra and isolated AVX2 probe

Bead `leopard-79h.18.20.2`, a bounded prerequisite under existing Experiment T.
**Algebra and isolated code generation passed; no measured speedup or codec
integration.** This is a materially different direction from the completed
[negative adjacent-scheduling screen](avx2_adjacent_screen.md), not a retry.

## Actual field representation

The first eight coordinates of Leopard's actual GF16 Cantor basis form a
GF256 subfield. With `u = β8` (canonical symbol `0x0100`), the relation is
`u² = u + δ`, where `δ = β7 = 0x80` in that subfield. The quadratic is
irreducible: none of its 256 subfield elements satisfies `t² + t = δ`.
This uses GF16's polynomial `0x1002D`, **not** an assumed reuse of GF8 tables.

Represent a symbol as `a + b*u`, with subfield bytes `a,b`. For a fixed
multiplier `c + d*u`, three GF256 products suffice:

```
A = a*c
B = b*(δ*d)
C = (a+b)*(c+d)
product = (A+B) + (C+A)*u
```

All additions here are XOR. Canonical-to-tower conversion is derived from
the complete map `u*b`, not by assuming the canonical high byte already is b.
`high_byte(u*b)` is a permutation. Its inverse recovers b from the canonical
high byte; then `a = canonical_low XOR low_byte(u*b)`. The reverse conversion
XORs a into the low byte of `u*b`.

The eight canonical `u * (1<<i)` basis images are
`0100, 02cf, 04ab, 0821, 108a, 209d, 4027, 801f` (hexadecimal).
Both conversions are GF(2)-linear and reversible over all 65,536 symbols.

## Isolated code-generation result and costs

The probe operates on complete 64-byte ALTMAP blocks, with 32 low and 32 high
bytes. Product tables occupy 96 bytes per coefficient (six 16-byte rows),
versus the existing eight-row, 128-byte GF16 nibble representation. A nominal
65,536-slot table would occupy 6 MiB rather than 8 MiB; no full table cache or
initialization policy is implemented here. Each conversion uses 64 bytes of
fixed tables.

Actual GCC Release object, `-march=x86-64 -mavx2 -mno-avx512f -mno-gfni`:

| Loop, per 64-byte block | Instructions | Shuffles | Masks | Shifts | XORs | Stack references |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Standalone tower product | 30 | 6 | 6 | 3 | 6 | 0 |
| One conversion direction | 18 | 4 | 2 | 1 | 3 | 0 |

The object has no EVEX, high-vector-register, ZMM, GFNI or ternary-logic
instructions. Six product table broadcasts are outside its loop. This is a
standalone multiply, **not** a full forward/inverse/accumulating butterfly;
its total instruction count cannot be compared directly to those kernels.
Two fewer product shuffles do not imply a throughput gain: the third subfield
input adds masks/shifts, and conversion, loads, table setup and register live
ranges in actual butterflies still matter. The conversion loop even rereads
the source high vector through a memory operand. Converting around each
multiply would erase the intended shuffle saving; a useful design would keep
the representation internal across enough FFT work to amortize conversion.

## Validation actually completed

Per Release and full ASan/UBSan/LSan probe:

- All 65,536 GF256 input/coefficient products agree with independent GF16
  polynomial arithmetic and the existing original-Leopard2 scalar field API.
- All 65,536 symbols pass scalar canonical/tower roundtrips and vector
  conversion, including in-place two-block conversion. Linearity is checked.
- All 65,536 GF16 constants times all 16 canonical input basis vectors
  (1,048,576 products), plus 1,048,576 zero/all-one/deterministic extra cases,
  agree across polynomial oracle, original scalar field, scalar tower formula
  and actual AVX2 product with both boundary conversions.
- The product runs on unaligned one-block buffers with write canaries, repeats
  in place, accepts zero blocks without accessing null pointers, and detects
  six separately corrupted product-table rows. Six total CLI invocations reject
  options, including `--measure`, before initialization.

For each fixed coefficient the construction is linear in the input; exhaustive
input-basis agreement and the checked linear maps establish the scalar product
identity for all symbols. This is not a claim that 2^32 native products were
executed. The existing field API comparison is **Leopard2 scalar arithmetic**,
not a new full-byte Leopard1 encoder parity comparison.

The collector-free Python replay uses a separate shift/reduce polynomial
algorithm, rederives the subfield, conversion and all 1,048,576 constant/basis
products, verifies native records and binary/source hashes, and compares saved
disassembly to the actual object. Eight pure adversarial tests pass normally
and with Python `-O`, including changed shuffles, EVEX/high registers, stack
references, missing loop/function, and a singular basis. Raw and read-only-copy
replays pass in both modes with the identical SHA-256:
`4bd667a35d8ca0284b4cb3b8c0d23d08c93619841bf68bbcceec8ad9099814f2`.

Original production sources/archives were not rebuilt or edited. The original
scalar field source is byte-identical to `45e2eff` (`fd27e72c…`). The Release
oracle archive is `89f33d3d…`; the reused full-sanitizer archive is `501d5b03…`
from adjacent qualification, with the same field implementation. The new
Release product/conversion object is
`f2adbdc0dd192d1be74fe4bee4d4c1bf90a53752ed4ad4d241d1198d5d73a35e`.

| Scope | Peak bytes | Cap | Outcome |
| --- | ---: | --- | --- |
| Initial build | 180,490,240 | 512 MiB | exit1: sanitizer non-PIE archive/PIE link mismatch |
| Corrected sanitizer probe build | 144,424,960 | 512 MiB | exit0 |
| Both native profiles and CLI checks | 28,479,488 | 256 MiB | exit0 |
| Raw normal/optimized replay | 25,108,480 | 256 MiB | exit0 |
| Retention, both sealed replays and tests | 73,494,528 | 256 MiB | exit0 |

Every listed scope has all six memory events zero and swap zero. The initial
link failure and its exact sources/objects/logs remain retained. It was fixed
by using the existing sanitizer non-PIE/non-recovery policy for the two new
probe objects and link, without changing or rebuilding the codec archive.

Raw: `/tmp/leopard-tower-algebra.T3NxrT`.
Read-only: `.research/leopard-79h/tower-algebra-qualified.jWRjOr`, 81 files /
40,981,118 bytes plus outer manifest. Manifest SHA-256:
`b0a35bb9701d62b645d4e0145adcae6aa615721dea94765c2735d31d9aa60835`.
Delivery: `/tmp/leopard-tower-delivery.eJaDrq`.
See the [machine-readable result](results/tower_algebra_20260909.json).

## Remaining gates

Next is an **untimed** integration/cost audit and isolated full-butterfly
prototype: enumerate every canonical/tower boundary and all constant/special
skew handling, including copied sources, directly bound parity, partial outputs,
decode and other backends. This probe does not qualify any public codec,
compact tail, arbitrary overlap, scratch contract, multi-block product path or
concurrent codec behavior. Whole-transform tests and original Leopard1 parity
must precede any public performance experiment. Experiment T's eventual
10% end-to-end promotion threshold is unchanged; timing would need a new,
separately reviewed, committed-and-pushed preregistration with fixed controls.

Everything ran locally, serial under the canonical lock, with capped builds
and checks, no swap, no Claude, subagents, SSH workers or unrelated host changes.
Codex self-review and deterministic/adversarial checks follow the user's Claude
opt-out; no independent-model `CONVERGED`. The performance goal stays open.

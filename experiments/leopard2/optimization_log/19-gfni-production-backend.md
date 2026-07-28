# 19 — GFNI promoted to a production backend member

**Disposition: SHIPPED — default builds now carry the 1.33x-1.73x GFNI codec as
`LEO2_BACKEND_GFNI`, explicit opt-in, AUTO unchanged**

## Why this was the top-ranked performance item
The paper-conformance audit (report 18's sibling work) surfaced a scoping fact:
`LEO2_GFNI_VARIANT` was defined by **no CMake target**, so the campaign's
largest kernel wins — GFNI affine multiplication (report 01), radix-eight
staging (report 03), and their interaction with GF16 fusion (report 02) —
compiled out of every default build.  Reports 01/03/04's measured 1.15x-1.67x
delivered nothing to an actual user.  Promoting the existing, audited kernels
was worth more than any new optimization.

The repo's own gate said exactly what promotion required
(`docs/leopard2_gfni_codec.md`, requirement 1): an own translation unit, own
`leo2_backend` identity, CPUID/XCR0 qualification, startup known-answer test,
and own allowed-instruction list — which is precisely what
`leopard2_portable_isa` was enforcing every time it rejected `vgf2p8affineqb`
in the AVX2 object.

## What shipped
- **`LEO2_BACKEND_GFNI = 6`** appended to the public enum (additive ABI, same
  pattern as AVX512's addition at 5).  Ops name `"avx2-gfni"`.
- **`Leopard2BackendGFNI.cpp`**: re-includes the AVX2 source with
  `LEO2_GFNI_MEMBER` + `LEO2_GFNI_VARIANT`, compiled `-mavx2 -mgfni
  -mno-avx512f` (GNU/Clang only, mirroring AVX512; MSVC has no GFNI switch
  below /arch:AVX512).  The TU owns its affine tables outright — the
  double-8 MiB hazard the conformance audit flagged is avoided by construction,
  and the shared-table accessors are compiled out of the member TU.
- **Qualification**: `X86Features.gfni` = AVX2 gate + CPUID.(7,0):ECX bit 8;
  the pure classifier grew a `leaf7_ecx` parameter with synthetic tests.
  Explicit-only: `AUTO` never selects GFNI.
- **KAT**: `TestFF8Butterfly8Out` — both fused radix-eight kernels verified at
  startup against a reference composed from the audited two-point butterfly,
  closing the audit's "no KAT names butterfly8" gap.
- **Portable-ISA `gfni` class**: AVX2 VEX set + `vgf2p8affineqb` + four
  plain-AVX2 mnemonics only the GFNI-only kernels emit; EVEX/ZMM rejected;
  `vgf2p8mulb` explicitly rejected; leak fixtures mutation-tested.
- **Policy mapping — the subtle part.**  Thirty `kind == AVX2`-family predicates
  (tiling from reports 13-16, R=1 XOR coarse/pairwise/fused-final, direct
  encode/repair, sink metadata, syndrome fusion, forward fusion, decode
  dispatch) were extended to the new identity site-by-site, because the
  experimental configuration ran *under the AVX2 identity* and every measured
  number depended on those policies being active.  Two deliberate exclusions:
  the AUTO-AVX512 encode exception (unreachable for explicit contexts) and the
  AVX512-only `ff8_high_encode_one_block` gate (the GFNI table does not publish
  that op).  Nine pinned source-text contracts in
  `tools/leopard2_operation_counts.py` were updated to the extended predicates.

## Validation
- **Full suite: 103/103** (grew from 100: three GFNI fault-injection tests) —
  including `leopard2_portable_isa` **passing on a GFNI-bearing archive for the
  first time**, `leopard2_field_options_matrix`, `leopard2_test_hook_isolation`,
  the vcxproj contract (134/134 after one documented exemption extension), and
  the operation-counts model (397 checks).
- **Randomized backend differential: 200 shapes, 0 mismatches** — avx2 vs gfni
  in one binary, digests + round-trip + acceptance; 42 shapes rejected
  identically.
- **Same-binary performance**, one pinned core, default build flags:

| Cell | AVX2 | GFNI | speedup |
| --- | ---: | ---: | ---: |
| GF16 K=200 R=50, 64 KiB | 1333.1 us | 769.6 us | **1.732x** |
| GF16 K=1000 R=200, 64 KiB | 9106.0 us | 6016.4 us | **1.514x** |
| GF16 K=2000 R=500, 256 KiB | 92582.1 us | 63561.2 us | **1.457x** |
| GF8 K=192 R=64, 64 KiB | 725.1 us | 505.9 us | **1.433x** |
| GF8 K=224 R=32, 64 KiB | 803.0 us | 603.4 us | **1.331x** |

GF8 K=224 R=32 reproduces the experimental 1.33x to three digits.  The 256 KiB
row doubles as the policy-mapping witness: both members tile there, and GFNI
wins by the kernel margin — had tiling silently detached from the new identity,
that row would have shown GFNI losing.

## What the KAT caught on its first run
The composed reference initially used the branchless two-point butterfly at
sentinel skew 255, which multiplies by **one** (`AddMod` wraps log 255 to the
identity).  The fused radix-eight contract instead **absorbs the skip branch**
that two-point callers perform before dispatching — at the sentinel the multiply
half is skipped outright.  The KAT failed, the divergence was traced to the
documented contract, and the reference now encodes the skip.  A KAT that fails
informatively on its first contact with the kernel is doing its job.

## Versus leopard1
Per the standing method (report 18): in a default-flags build the in-process
leopard1 runs the plain AVX2 kernels — leopard1 at its best available code —
so this measures the promoted member against a leopard1 that does NOT get
GFNI.  `legacy_comparison` reported `matched` for both backends in all 11
cells.  Ratios are leopard1 time / Leopard2 time; decode includes leopard1's
per-call locator setup against Leopard2's amortized plan (reuse 8).

| Cell | AVX2 enc | GFNI enc | AVX2 dec | GFNI dec |
| --- | ---: | ---: | ---: | ---: |
| GF16 K=1000 R=200, 256 KiB | 2.218 | **3.302** | 6.45 | **9.41** |
| GF16 K=300 R=100, 1 MiB | 2.295 | **3.365** | 4.53 | **6.78** |
| GF16 K=2000 R=500, 64 KiB | 1.737 | **2.139** | 5.12 | **7.49** |
| GF8 K=224 R=32, 64 KiB | 1.071 | **1.877** | 2.08 | **2.50** |
| GF16 K=255 R=129, 64 KiB | 0.991 | **1.705** | 2.41 | **3.28** |
| GF16 K=200 R=50, 64 KiB | 0.988 | **1.668** | 3.42 | **5.00** |
| GF16 K=1000 R=200, 64 KiB | 1.096 | **1.662** | 5.36 | **8.35** |
| GF8 K=128 R=128, 64 KiB | 1.111 | **1.514** | 1.49 | **1.75** |
| GF8 K=192 R=64, 64 KiB | 1.057 | **1.398** | 1.72 | **2.09** |
| GF8 K=240 R=16, 64 KiB | 1.048 | **1.190** | 2.78 | **3.14** |
| GF16 K=2000 R=500, 256 KiB | 3.924* | 3.344* | 4.99 | **8.38** |

**Medians: AVX2 enc 1.096 / dec 3.42 — GFNI enc 1.705 / dec 5.00.**

\* The K=2000 R=500 256 KiB encode pair inverts only because each column
divides by its own child's leopard1 run, and leopard1's ~250 ms encode at that
cell is noisy at 5 iterations; the same-binary comparison at this exact cell is
1.457x in GFNI's favour (table above).  The decode pair, measured in the same
children, shows GFNI ahead as expected.

The headline: cells where Leopard2's AVX2 member sits at **parity** with
leopard1 on encode (0.99-1.10) run **1.4-1.9x** on the GFNI member, and the
best cells reach **3.3-3.4x encode / 9.4x decode** — in a default build, with
byte-identical output.  Screen evidence, not promotion evidence: the isolated
`run_abba.py` campaign remains required before these enter
`docs/leopard2_vs_main_benchmark.md`.

## Left explicitly open (requirements 2-6 of the GFNI doc)
Packed GF16 tables (8 MiB -> 2 MiB, footprint-only), AUTO selector policy
(needs the same same-binary screen the AVX-512VL callback received), the
isolated `run_abba.py` exact-main campaign, a second microarchitecture, and the
512-bit width follow-up.

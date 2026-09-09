# GF16 tower butterflies: untimed qualification and boundary audit

Tracker `leopard-79h.18.20.3`, review child `.3.1`, 2026-09-09.
Experimental kernels only. No production routing, codec integration, or timing
change. The previous [algebra milestone](tower_algebra.md) remains unchanged.

## Native result and actual code generation

Release and full ASan/UBSan with leak detection each pass 262,144 butterfly
cases (65,536 log values times four forms), 1,560 boundary/tail cases, 65,536
two-shuffle conversion values and nine zero-byte/null calls. Six total CLI
requests refuse timing/help/numeric arguments before initialization.

The four forms are in-place inverse, in-place forward, disjoint forward-out,
and disjoint inverse accumulation. Original scalar and AVX2 operations agree
with independent polynomial multiplication and converted tower results.
Accumulating twice restores outputs; out-form inputs remain unchanged.
All32 symbol-pair basis inputs are distributed across32 lanes for each log,
not every basis at every lane or every possible block. Boundary lengths are
2,30,32,62,66,126,128,130,8190,8192,32768,32770,65536,65598, at offsets0/1/17;
every even compact tail after128 bytes is covered separately. ASan cannot
poison up to7 bytes immediately before misaligned data; canaries cover writes
there, and aligned prefixes provide exact poison. No stronger claim is made.
The log-to-coefficient mapping comes from the original field API, not an
independently reconstructed log table.

The actual pure-AVX2 object SHA-256 is
`fc363f48a317469117d04c5c68519335daffb9c1e46a5717aa3c10a25e435446`.
GCC13.3 uses generic x86-64/AVX2, no AVX-512/GFNI, strict warnings. Ordinary
full64-byte loop counts, including loop control:

| Operation | Loop | Instructions | Shuffles | Masks / shifts / XORs | Loop stack accesses |
| --- | --- | ---: | ---: | --- | ---: |
| Shared forward pair/out | 0x150–0x1fa | 36 | 6 | 6 / 3 / 10 | 0 |
| Inverse pair | 0x480–0x523 | 37 | 6 | 6 / 3 / 10 | 0 |
| Inverse accumulation | 0x760–0x820 | 40 | 6 | 6 / 3 / 14 | 0 |
| Conversion involution | 0x338–0x379 | 15 | 2 | 2 / 1 / 2 | 0 |

Whole functions DO use stack: callee saves, stack arguments, and inverse
scalar-tail storage. Their respective syntactic stack-access counts are
12/22/19/4. These are not regular vector-loop spills. Both forward entrypoints
actually jump to the counted helper. The object contains no EVEX, high SIMD
registers, ZMM, opmask, GFNI, ternary logic or carryless multiplication.

The proper GF16 fixed-multiply comparison is eight baseline nibble shuffles
versus six tower shuffles per64-byte block, with all butterfly operands and
extra XOR/dependency work included. A GF8 multiply's four shuffles is not the
same field workload. Static instruction reductions are not measured speedups;
prologues, table preparation, tails and conversions are additional costs.

## Conversion and table costs

Independent shift/reduce polynomial replay verifies all256 values of
`high_byte(u*b)=b`. Therefore canonical/tower conversion is one involution:
keep the high byte, XOR the low byte with `low_byte(u*high)`. The linear map
needs two16-byte nibble rows, two shuffles, two masks, one shift and two XORs
per64-byte block. Both directions use the same32-byte table.

Product tables remain six16-byte rows (96 bytes) per fixed coefficient.
The prototype constructs96 subfield entries per coefficient, after mapping
the coefficient and forming c, delta*d, and c+d. A complete65536-slot table
would occupy6 MiB versus8 MiB for128-byte canonical entries. This is nominal
table storage only: if canonical fallback tables stay resident it is not a
2 MiB reduction in total memory. The test's64-KiB subfield table is an oracle
convenience, not an integrated initialization/cache strategy. Initialization
latency, actual cache residency and whole-codec memory use are unmeasured.

## Actual integration boundaries

Pinned sources: `LeopardFF16.cpp` fd27e72c…, `LeopardFF16.h` feed12b1…,
`Leopard2BackendAVX2.cpp` e1546599…, `Leopard2Backend.h` 7d4535a1…,
`leopard2.cpp` 93412028…. Full hashes and source copies are retained.

- `IFFT_DIT_Encoder_Impl` in LeopardFF16.cpp:1226 directly memcpy-copies
  large high-rate sources; swapping Ops cannot convert those inputs.
  The compact source-fused `ff16_ifft_butterfly4_out` path and partial source
  groups also need explicit canonical-to-tower boundaries. Zero padding stays
  zero. Conversion can potentially replace an existing source copy, but that
  has not been implemented or measured.
- `ReedSolomonEncodeWithSourcePolicy` at2263 keeps the first inverse result
  in the accumulator, then transforms later source blocks in temporary space
  and accumulates before the forward transform. All internal buffers and all
  corresponding callbacks must consistently use tower coordinates. Do not
  preconvert all K sources outside the existing2T workspace.
- Public high outputs in leopard2.cpp:22073 bind directly to caller recovery
  buffers, including byte-tiled calls. A final tower-to-canonical conversion
  is genuinely new work: there is no existing final output copy to remove.
  The staged64-byte compact tail near22150 needs conversion before
  `GatherShardFromKernel`. EncodeLayout already rejects odd GF16 byte counts.
- Ordinary multiplier log0 means one. Butterfly log65535 means zero skew;
  raw `MultiplyLogElement(1,65535)` instead returns one because the exponent
  table duplicates its endpoint. Preserve this contract distinction. Pair
  and accumulating baseline callers specialize the sentinel; forward-out
  accepts it. The new experimental wrappers accept it explicitly, with null
  tables allowed only for zero skew or zero bytes.
- Decode is not qualified by this encode audit: direct source staging,
  reveal multiplication, in-place multiply, direct recovery outputs,
  block-zero XORs and generic/low paths need their own boundary audit.
  Other backends and every public buffer must remain canonical.

No source audit presently rules out an experiment, but the two-shuffle saving
does not establish that new boundary work, cache effects or dependencies pay
off. Whole-codec/native-Leopard1 parity, decode, dispatch, public bounds,
scratch limits, both fields and concurrency remain separate gates. Any later
timing needs reviewed, committed AND pushed preregistration; Experiment T's
10% end-to-end promotion gate and2% controls/neighbors remain unchanged.

## User-requested independent static review

Local read-only Codex subworker `/root/tower_bug_review` performed three bounded
passes, without file edits, native jobs, Claude or remote hosts. It found one
low-severity bug in the prior algebra replayer: implicit stack instructions
were missed. Main fixed this and tested push/pop/pushfq/popfq/call/lcall/enter/
leave mutations. Nine verifier tests pass in both Python modes; fixed replay
matches the old result exactly. The original object has no such instructions,
so prior qualification is still valid and its eight-test historical report is
not rewritten. The reviewer also identified the two-shuffle opportunity.

No additional static correctness defect was found in the new kernels/harness;
the known nibble-helper name typo was corrected before compilation. Native
validation is main's execution evidence. This is not a Claude fixed-point
gate or independent-model `CONVERGED` claim.

Ten new adversarial replay tests pass normally and under Python `-O`, from
source and retained copies. Collector-free replays verify source/artifact
manifests, actual disassembly, native counts/exits, CLI refusals, resource
logs and independent conversion algebra. Four outputs match SHA-256
`6ea504f706d7c811c6e6f8fa1c89708d61a6fd5885fcb4de15b8bbac1afdf50f`.
They do not execute a codec. Prior product qualification is referenced,
not rerun; linking the prior standalone probe only resolves the unused
renamed previous test main.

All jobs were serial/local under the canonical lock, no swap, bounded CPU
and wall limits. Codegen build96,677,888/512MiB; native build178,720,768/512MiB;
native checks49,225,728/256MiB; new replay23,318,528/256MiB; prior-verifier
review65,683,456/256MiB; sealing70,250,496/256MiB. All exited0, all six memory
events zero. Existing original Release/full-sanitizer oracle archives are
unchanged (89f33d3d… /501d5b03…); no monolithic suite was run.

Raw: `/tmp/leopard-tower-butterfly.aQABl4`; verifier fix:
`/tmp/leopard-tower-review.ViXf9s`; sealing log:
`/tmp/leopard-tower-butterfly-seal.log`.
Readonly local bundle `.research/leopard-79h/tower-butterflies-qualified.YkqbBN`
contains89 files /42,984,375 bytes including its manifest. Manifest SHA-256:
`cb91cae127db3e0aa5bc872b5bacd5ac7c4d2a0de941e3fd3f09f3a6b679d65f`.

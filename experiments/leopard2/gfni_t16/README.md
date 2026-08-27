# Native-Cantor AVX-512/GFNI T16 encoder

This directory contains the bounded experiment and qualification tooling for
the exact legacy-high GF8 `K=16, R=16, T=16, B=64` encode operation.

The production candidate keeps the sixteen 64-byte shard rows in sixteen ZMM
registers for the complete inverse-plus-forward schedule.  It uses
`VGF2P8AFFINEQB` directly in Leopard's existing `0x11d` Cantor representation;
there is no conversion to AES representation and no change to public field,
coordinate, or wire semantics.  The leaf is compiled in a separate raised-ISA
object and is published only after AVX-512F/BW/VL, GFNI, XCR0, calibrated-host,
and startup-known-answer gates pass.  Explicit backends and every neighboring
shape retain the established implementation.

Files:

- `screen.cpp` is the direct-leaf selection screen.
- `public_screen.cpp` compares the enabled and disabled operation leaf through
  the public ordinary and one-item-batch APIs.
- `run_k16r16_b64_avx512_gfni_abba.py` is the frozen same-binary acquisition
  controller.  It gates the target on both ordinary and one-shot 95% lower
  confidence bounds of at least 1.05 and treats inactive-cell timings as
  descriptive only.
- `audit_k16r16_b64_avx512_gfni_abba.py` is an independent, no-timing replay
  auditor.  It does not import the acquisition controller and reconstructs the
  complete journal, raw launch identities, route evidence, statistics, gates,
  and source closure.  The authoritative wrapper uses its archive-only source
  mode so durable replay remains valid after the live branch advances; the
  acquisition controller and wrapper still check the live repository before
  and after the measured run.
- `run_authoritative_k16r16_b64_avx512_gfni.sh` creates fresh live and detached
  replay builds in an allowlisted environment, proves byte-identical build
  closure, requires the exact seven-test Release census, runs the frozen
  campaign, and publishes a sealed completion envelope only after the
  independent post-seal replay succeeds.  The completion record is written
  through a pre-opened descriptor only after every file and directory is
  read-only, and its precomputed digest is already in the envelope manifest.
  Its `--verify` mode rechecks the sealed envelope without collecting timing.
- `scheduler-quiet-selection.json` retains the first scheduler-only CPU-pair
  selection, while `scheduler-quiet-selection-v2.json` retains the final
  full-node selection after pair 10/74 later failed its contamination gate.

The short screens are diagnostic selection evidence, not promotion evidence.
Production default-on status is conditional on a completed, immutable
live-versus-replay campaign and a successful independent replay audit.

## Authoritative result

That qualification completed in the sealed V6 envelope at
`.research/leopard-79h/f14aeead-t16-final-v6`, bound to source commit
`f14aeead4ad3c65897c9990e47164a65ecd608a9` and tree
`797d51815d8f5bf237fbcf712a41758cf518a7ac`.  Against the same binary with
only this operation leaf disabled, the exact `K=16, R=16, T=16, B=64`,
legacy-high GF8 AUTO workload measured:

- ordinary encode: `3.3800831324506677x`, 95% confidence interval
  `[3.3303699132991933, 3.430538432578953]`;
- one-shot encode: `3.7003885797783203x`, 95% confidence interval
  `[3.6607625548212575, 3.74044353773237]`.

All 49 ABBA rounds and 196 launches were accepted with zero contamination
retries.  All route, call-count, digest, isolation, inactive-neighbor, and
timer-floor gates passed.  The frozen verifier and two independent no-timing
replays reproduced every statistic and the complete source/build/evidence
closure.  The retained campaign, binary, and completion SHA-256 values are,
respectively, `7ad53caedfd1fd9668d240b3549e4b7ea8476a1cef9b68164fd4c0db4a98b3be`,
`67f2b27961b9914ca1493ad1854d57b8e83a77b2c4cef81383467bfbe41284ea`,
and `30b293d3c4fdfe814afee31c64f8b41acdd8226dc9cde9bc521c80fd2dcae69d`.

This promotes only the optional runtime-qualified leaf on calibrated AMD
family `1Ah`, model `08h` hosts with AVX-512F/BW/VL, GFNI, and OS ZMM state.
Explicit backends and all other shapes retain their prior path.  The result is
not an exact-main Leopard1 comparison and does not establish that Leopard2 is
faster for every configuration.

The final authoritative campaign is pinned to CPU 52 with SMT sibling 116.
CPU 13/sibling 77 first exhausted the frozen contamination budget, a longer
prereacquisition check rejected CPU 9/sibling 73, and the first retained scan
selected CPU 10/sibling 74.  A complete campaign on 10/74 passed its measured
gates but stopped at a subsequently repaired archive-auditor defect; the next
fresh campaign then exhausted the unchanged contamination budget on that pair.

The final selection used a 240-by-250-ms `/proc/stat` scan of every remaining
online sibling pair.  Because the first scan's relative-minimum pair later
failed the contamination gate, the final rule was deliberately tightened to
require zero nonidle ticks on both logical CPUs.  Prior failed/rejected pairs
and the scanner's own 63/127 pair were excluded before the scan.  The fixed
rule selected the lowest primary CPU among eligible pairs: 52/116 and 54/118
qualified, so 52/116 was selected.  Full counters are retained in
`scheduler-quiet-selection-v2.json`.  No candidate timing was collected during
either CPU-selection scan, and no retry, timing, or statistical gate changed.

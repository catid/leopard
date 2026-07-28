# Session resume state — 2026-07-28, branch `codex/leopard2` at `2942f35`

Everything below is also tracked in Beads (`leopard-79h.43` through `.48`).
**Use the 1.x binary explicitly**: `/home/catid/.nvm/versions/node/v24.1.0/bin/bd`
(`~/.local/bin/bd` is a forbidden 0.47.0 legacy install; `type -a bd` first,
per CLAUDE.md).

## What is done and durable

Commit `2942f35` (pushed to `origin/codex/leopard2`): the whole optimization
campaign. Narrative + evidence live in
`experiments/leopard2/optimization_log/README.md` (index of reports 00-23) and
`docs/leopard2_gfni_codec.md` / `docs/leopard2_open_optimization_findings.md`.
Highlights: GF16/GF8 live-set workspace tiling (decode to 3.38x, encode to
2.92x, scratch to 126x smaller), production `LEO2_BACKEND_GFNI = 6` member
(vs leopard1 medians enc 1.10->1.71, dec 3.42->5.00, byte-matched), dead
block-0 schedule removal (plan setup 1.06-1.15x), writable-partitioned batch
preflight, IT-2026 conformance audit clean. Validation at the commit: full
suite 103/103 four consecutive times, `leopard2_portable_isa` passing on the
GFNI-bearing archive, vcxproj 134/134, operation-counts 397 checks, ~2,200
randomized differential shapes with zero mismatches.

## In-flight: the authoritative exact-main campaign (leopard-79h.43, P1)

FINAL STATE: three model-drift layers were fixed and pushed (`a972a27`): the
TU closure (AVX2Xor + GFNI), the library `LEO2_HAVE_GFNI_BACKEND` define, and
the bench source-attestation defines (pre-existing since `a306d95`). All 22
candidate compile entries verified byte-for-byte against the model. The
campaign then PASSED every identity check and reached cell execution (first
baseline child output retained under the evidence dir) — and was SIGTERMed by
the affinity supervisor: an interactive Claude session's own tool/monitor
process churn under the same UID trips the fail-closed containment ("new
process has an already restricted mask but no safe uniform creator
provenance"). **The campaign must be run from a bare terminal with no Claude
session active on this host**: delete `/tmp/leopard2-vs-main-evidence-2942f35`
and `/tmp/campaign-2942f35.status`, refresh the report-dir suffix inside
`/tmp/launch-campaign-2942f35.sh` (the sed line pattern is in the script
history), run it, wait for `/tmp/campaign-2942f35.status`, then
`run_abba.py verify` the bundle. Candidate = production build at `a972a27`
(rebuilt so the baked bench attestation matches `--candidate-commit`).
The original drift diagnosis is retained below for context:

- Error: `candidate compile-command entry closure differs`.
- Cause (verified against `/tmp/leopard2-production/compile_commands.json`,
  22 actual entries vs 20 modeled): `CANDIDATE_LIBRARY_SOURCES` in
  `experiments/leopard2/main_compare/run_abba.py` (~line 175) is missing
  `Leopard2BackendAVX2Xor.cpp` (pre-existing drift) and
  `Leopard2BackendGFNI.cpp` (new this session).
- Exact fix (all in `run_abba.py`): add both to `CANDIDATE_LIBRARY_SOURCES`;
  `backend_targets` (~1528): `AVX2Xor -> leopard2_backend_avx2.dir`,
  `GFNI -> leopard2_backend_gfni.dir`; `isolated_flags` (~1567):
  AVX2Xor `['-mavx2','-mno-avx512f','-falign-functions=64']`, GFNI
  `['-mavx2','-mgfni','-mno-avx512f','-falign-functions=64']`; a `definitions`
  elif mirroring AVX512's (~1609): GFNI needs
  `['-DLEO2_HAVE_AVX2_BACKEND=1','-DLEO2_HAVE_GFNI_BACKEND=1']`, AVX2Xor none.
  Cross-check each constructed argv against the real compile_commands entries
  before rerunning (`expected_compile_argv` orders tokens as: definitions,
  includes, -Wall -Wextra -fopenmp -O3 -DNDEBUG -O3 -std=gnu++11, isa flags,
  propagated openmp).

Ready-made campaign artifacts (all still on disk):

| Piece | Path |
| --- | --- |
| Baseline (adapter-built, pinned `6e5725e`) | `/tmp/leopard-main-compare-2942f35/` (`leopard_main_benchmark`, `libleopard_main_exact.a`) |
| Baseline source root | `/tmp/leopard-main-exact` (clean at `6e5725e`) |
| Candidate (production, tests OFF) | `/tmp/leopard2-production/` at `2942f35` |
| Reservation (cpu 14 / sibling 30) | `/tmp/leopard2-main-reservation.json` |
| Launcher template | `/tmp/launch-campaign-2942f35.sh` |

Launch gotchas learned the hard way: the runner refuses a pre-existing
`--output` dir; the affinity supervisor refuses a pre-existing `--report`
path; the baseline must be adapter-built from *this* repo's
`experiments/leopard2/main_compare` (retained July-22 artifacts are rejected);
monitor completion via a status file the wrapper writes, never via `pgrep`
(self-matching). Keep the whole machine quiet during timing — the codec is
DRAM-bound, and the sibling-idleness gate only polices CPUs 14/30.

## Then (dependency order)

1. **leopard-79h.44 (P1)**: explicit-GFNI campaign pass; if campaign-grade
   evidence matches the screens, flip `LEO2_EXPERIMENT_AUTO_GFNI` on by
   default (mechanism + calibrated Zen 5 pin already shipped default-off) and
   update the `leopard2.h` AUTO doc text.
2. **leopard-79h.48 (P3)**: transcribe verified campaign numbers into
   `docs/leopard2_vs_main_benchmark.md` (blocked on .43/.44). Mind report 18:
   two different leopard1 baselines exist; the campaign's external pinned
   baseline is the correct one for that doc.
3. **leopard-79h.45 (P2)**: GF16 affine table packing 8 MiB -> 2 MiB.
4. **leopard-79h.46 (P2)**: derive tiling live-set constants from detected L3
   (constants are calibrated for the 32 MB CCD1; this host's CCD0 is 96 MB).
   Do NOT thread-scale the target — refuted, report 17.
5. **leopard-79h.47 (P3)**: second-microarchitecture (Intel GFNI) validation.

## Standing cautions for the next session

- Screen evidence never enters `docs/leopard2_vs_main_benchmark.md`; only
  verified `run_abba` bundles do.
- The method notes in `experiments/leopard2/optimization_log/README.md` are
  the distilled failure modes of this campaign (stale stored evidence,
  unchanged-witness non-results, silent guard widening, component-share
  overstatement, self-matching monitors). Read them before measuring.
- `/tmp/l2-*` build trees are disposable screen infrastructure; the only
  campaign-grade artifacts are the table above.

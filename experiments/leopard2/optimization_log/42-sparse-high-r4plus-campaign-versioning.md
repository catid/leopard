# 42 — Sparse-high R4+ production campaign versioning

## Result

The fresh ordinary-archive qualification runner is versioned and ready for a
new acquisition. It preserves every v10 evidence identity and adds a distinct
v11 dialect for the mechanism-bounded R in `{4,8,16}` candidate. This change
does not reuse the passing cells from the earlier campaign, run the fresh
timing matrix, enable the compiled AUTO default, or authorize promotion.

The current 88-cell grid contains:

- 56 candidates: 24 one-shot tuples, eight parity-row samples, and 24
  batch/binding API lanes;
- 18 padded-side-two structural controls: 12 one-shot tuples and six
  `K=16,R=2,4096` batch/binding lanes;
- the original 14 requested-AVX2 and AUTO/thread-four performance controls.

The v11 cells retain the v10 list order. Only the former R=2 candidate region
labels change, to `control_padded_side2_tuple` and
`control_padded_side2_public_api`. The frozen grid digest is
`5cb8acb9cdb2fc00de22753c9789f59bed4b515d0f62a663d45c2934c3e8b627`.

## Fail-closed evidence contract

| Evidence layer | Retained campaign | Fresh campaign |
| --- | --- | --- |
| outer manifest | v10 | v11 |
| job | v10 | v11 |
| analysis | v6 | v7 |
| cell | v3 | v4 |
| raw benchmark | v1 | v2 |
| campaign | qualification-v1 | padded-side-ge4-qualification-v2 |
| ordered grid digest | `ff6ab5...b6d2` | `5cb8ac...b627` |

New acquisition selects only the fresh dialect and a new default result
directory. Retained v10 manifests dispatch to their original 74-candidate,
14-control grid, raw route semantics, schemas, campaign name, and digest.
The v10 outer contract is pinned explicitly to build-attestation v14,
configuration-file v16, and the v7 launcher/environment family rather than
mutable current aliases. Cross-generation cell labels and raw schemas are
rejected.

For v11, the raw validator mirrors the production selector rather than
assuming that every table-enabled shape is direct-capable. Every R=2 policy
invocation must report zero prepared generator rows, transform selection, zero
direct witness calls, and the exact expected transform-call count. The
requested-AVX2 and thread-four controls remain at R=16, so table-enabled arms
still report R prepared rows while retaining transform routing.

The 56 candidates retain the preregistered route and setup-inclusive net lower
confidence-bound gates at five percent. The original 14 performance controls
retain the minus-two-percent net point guard. R=2 timings remain in the raw,
job, regional, and structural summaries, but do not act as performance evidence
or a noise-sensitive performance veto; all 18 structural records must instead
be present, valid, correct, and route-consistent for the v11 gate to pass.

## Verification

The deterministic runner self-test passed under normal Python and `python -O`,
including exact v10/v11 grids, cross-schema rejection, v10 R=2 replay, v11 R=2
transform truth tables, frozen v10 configuration/launcher/environment dispatch,
analysis classification, threshold immutability, and the existing
provenance/process-containment mutation suite.

Two independent read-only Codex reviews reached convergence. One review found
that the v16 configuration-variable tuple was still accepted only through the
mutable current alias; the corrected runner now names v16 explicitly and its
self-test advances the current aliases synthetically before replaying the
v14/v16 dialect. The final reviews independently confirmed byte-identical v10
contracts and the exact v11 role counts, grid digest, and R=2 route truth table.

On `foureyes.lan`, a source-only worktree copy passed six focused GCC 13 Release
tests: high-direct production, sparse-high AUTO production, direct encode,
benchmark-only registration, and both normal and assertion-disabled raw-schema
smokes. A separate ASan+UBSan build passed the same five executable tests (the
nested registration build excluded), with no sanitizer diagnostics. These are
correctness gates, not timing evidence.

Report 41 supplies the mechanism attribution behind the R>=4 boundary. The
remaining work is to acquire and replay all 88 fresh cells from sealed binaries
under the canonical timing and CPU-pair locks, then make a decision from that
new corpus only.

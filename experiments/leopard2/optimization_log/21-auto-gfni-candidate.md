# 21 — AUTO-GFNI: screen passed, candidate shipped default-off

**Disposition: MECHANISM SHIPPED (default-off), SCREEN EVIDENCE RECORDED.
Flipping the default awaits the isolated exact-main campaign, per the repo's
own AVX-512 precedent.**

## Why AUTO is the highest-reach item
Report 19's explicit `LEO2_BACKEND_GFNI` member only helps callers who select a
backend.  Most callers use AUTO, which deliberately reports the AVX2 baseline.
`docs/leopard2_gfni_codec.md` requirement 3 defines the promotion bar: "a
same-binary screen across the parent/redundancy/shard-size space", the same
treatment the AVX-512VL encode callback received before its calibrated-host
exception.

## The screen
44 cells, one binary, explicit `avx2` vs explicit `gfni` (identical dispatch
except the table), both fields, K = 2..4096, low and high rate (including
R > K), 1 KiB - 1 MiB, best-of-3 x 5 iterations, pinned core:

- **0 digest mismatches** across all 44 cells.
- **0 encode regressions** below 0.98x; encode gains up to 1.81x.
- **1 flagged decode cell** (GF8 K=129 R=1 at 64 KiB, 0.831x) — re-measured
  over four clean trials: 1.253, 1.001, 0.987, 1.001.  Noise.  The mechanism
  check agrees: R=1 decode is pure XOR reconstruction, and the two members
  share identical XOR kernels and identical R=1 dispatch policies, so a real
  17% delta had no available cause.

Verdict: GFNI >= AVX2 across the entire measured space on this host.

## What shipped
- `IsCalibratedAutoGFNIProcessor/Host` — the same Zen 5 Granite Ridge pin
  (family 0x1a, model 0x44) as the AVX-512 encode exception, with the screen
  recorded at the predicate.
- A `SelectBackend` arm behind **`LEO2_EXPERIMENT_AUTO_GFNI`** (default-off):
  on a calibrated, GFNI-qualified host the process default becomes the GFNI
  table, so AUTO contexts inherit it; unknown models retain AVX2 even with the
  experiment enabled.  Verified live: the default build's AUTO resolves to
  `avx2`, the experiment build's to `avx2-gfni`.
- A real bug the experiment caught before it could matter: `Initialize`'s
  startup if-chain had no GFNI arm (only the lazy explicit-request path did),
  so a GFNI *process default* failed with OUT_OF_MEMORY.  Explicit-backend
  usage never hits that path, which is why the promotion gate could be
  103/103 while the future AUTO state was broken.  Fixed and covered by the
  gate.

## Deliberately not done
The default is not flipped.  The repo's evidence bar for changing AUTO
(requirements 3+4) demands the isolated `run_abba.py` exact-main campaign with
the CPU-pair lease, build closure, and sibling-idleness gate — my screen is
explicitly screen-grade.  The AVX-512 precedent went through exactly this
sequence, and it widened only one operation on one cell; this candidate, with
its broader evidence, still owes the same campaign.  One `-D` flag and one
campaign stand between the current state and AUTO users getting 1.3-1.8x.

## Interaction noted for the campaign
With AUTO defaulting to GFNI, `CodecMayUseAutoAVX512Encode` goes inert (it
requires the AVX2 baseline).  That is the correct outcome — GFNI at the
AVX-512-widened cell measures 1.514x over AVX2, superseding the ~10-13%
widening — but the campaign should confirm it rather than assume it.

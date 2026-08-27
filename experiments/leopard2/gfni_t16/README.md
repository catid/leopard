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

The short screens are diagnostic selection evidence, not promotion evidence.
Production default-on status is conditional on a completed, immutable
live-versus-replay campaign and a successful independent replay audit.

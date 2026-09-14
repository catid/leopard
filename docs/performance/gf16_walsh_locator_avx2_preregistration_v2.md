# GF16 active-parent Walsh locator AVX2 same-process screen

This successor protocol follows the retained inconclusive campaign described
by `gf16_walsh_locator_avx2_preregistration.md`. That campaign's separate
scalar/AVX2 processes produced noisy direct-path controls and one non-idle
SMT tick; its results are not pooled or promoted. This protocol changes the
measurement boundary materially: both active-path implementations are timed in
one process with process-internal ABBA ordering.

## Frozen workloads and method

The six cells remain `(n, erasures) = (32,32), (64,40), (256,64),
(1024,64), (4096,256), (65536,8192)`, with seed `17`, 8 calls per sample,
31 measured samples, and 5 warmup samples. The `--backend compare` benchmark
initializes one qualified AVX2 backend, verifies scalar and AVX2 outputs against
the independent full-field reference, then alternates scalar-active and
AVX2-active setup in ABBA order within the same process. Direct locator setup
is not a timing acceptance gate here: the candidate callback is reachable only
from the dense active setup path, while direct setup and its existing tests are
unchanged.

## Resource and identity controls

- Build from a clean checkout of this pushed commit; copy the executable into a
  lane-owned immutable artifact directory and hash it before and after timing.
- Hold `${ATLAS_TMP}/leopard-gf8-authoritative.lock` directly with a nonblocking
  exclusive `flock`; do not add an outer lock around a runner that owns it.
- Pin to physical CPU 0 and reserve SMT sibling CPU 64. For every cell record
  sibling **non-idle** jiffies from `/proc/stat` fields
  `user+nice+system+irq+softirq+steal` (fields 2,3,4,7,8,9); idle and iowait
  are excluded. Any nonzero sibling delta invalidates that cell.
- Use `OMP_NUM_THREADS=1`, `OMP_DYNAMIC=FALSE`, record maximum RSS and exit
  status, and retain raw JSON and runtime logs for every cell.

## Acceptance

All six cells must have correct scalar and AVX2 outputs, clean source and
stable executable identity, zero reserved-sibling non-idle activity, and
successful resource/exit checks. The AVX2 active median must be at least 1.10x
faster than the scalar active median in every cell. Report raw medians, MADs,
RSS, sibling deltas, and the complete command/provenance. This setup result,
even if accepted, does not claim an end-to-end Leopard1 speedup and does not
change any AUTO route.

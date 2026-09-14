# GF16 active-parent Walsh locator AVX2 screen

This is the preregistration for the setup-only comparison of the qualified
pure-AVX2 `FF16WalshLocator` callback added in `b536ea7`. It is not a claim of
promotion or of an end-to-end codec speedup.

## Comparison

The candidate is the exact same pushed source and benchmark executable run
with `--backend avx2`; the control is `--backend scalar`, which calls the
existing scalar active-parent Walsh implementation. Both paths are checked
against the independent full-field/reference locator before timing. Direct
locator setup is measured in the same process as a fixed-work control, but is
not pooled with the active-path result.

## Frozen workloads

| parent `n` | erasures |
|---:|---:|
| 32 | 32 |
| 64 | 40 |
| 256 | 64 |
| 1,024 | 64 |
| 4,096 | 256 |
| 65,536 | 8,192 |

Every cell uses seed `17`, 8 calls per sample, 31 measured samples, and 5
warmup samples. The benchmark's per-process order is ABBA (direct/active), and
the process-level backend order is scalar, AVX2, AVX2, scalar. No samples may
be pooled across cells or across source/executable identities.

## Resource and identity controls

- Build from a clean checkout of the pushed commit; copy the executable into a
  lane-owned immutable artifact directory and hash it before and after timing.
- Hold `${ATLAS_TMP}/leopard-gf8-authoritative.lock` directly with a nonblocking
  exclusive `flock`; do not wrap a runner that acquires this lock in another
  lock.
- Pin the benchmark to physical CPU 0 and reserve its SMT sibling CPU 64.
  Record sibling **non-idle** jiffies before and after each cell from
  `/proc/stat` fields `user+nice+system+irq+softirq+steal` (fields 2, 3, 4,
  7, 8, and 9); idle and iowait fields are excluded. Reject the cell if this
  reserved-sibling activity delta is nonzero.
- Use `OMP_NUM_THREADS=1`, `OMP_DYNAMIC=FALSE`, and no other benchmark jobs.
  Record maximum resident set size and process exit status for every run.

## Acceptance and reporting

Accept the setup candidate only if every run is correct and provenance-stable,
the artifact hash is unchanged, the reserved sibling is idle, and the AVX2
active median is at least 1.10x faster than the scalar active median in at
least four of six cells without a direct-path regression beyond 2%. Report
all six cells, raw medians/MADs, RSS, control status, and limitations. A
failure or inconclusive control is retained and does not authorize pooling,
retrying, or relaxing the gates. Even an accepted setup result does not imply
an end-to-end Leopard1 comparison or change any AUTO route.

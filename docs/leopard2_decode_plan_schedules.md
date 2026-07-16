# Leopard2 reusable decode-plan schedules

Leopard2 decode execution no longer derives transform dependencies or input
prefixes from shard pointers. The immutable decode plan computes those values
after deterministic survivor and virtual-erasure selection, at the same time as
locator setup. The retained `ReedSolomonDecode*Prepared` functions still build
their schedules from masks and pointers and remain differential oracles; new
plan execution calls the `*Planned` companions.

## Compact output dependencies

The final fused-two-layer FFT asks whether an output is needed in groups of
`N`, `N/4`, `N/16`, and so on. When `log2(N)` is odd, it finishes with groups
of two. The plan concatenates one bit per queried group without per-level
padding. For level `j`, the bit offset is `(4^j - 1) / 3`; the coordinate's
group index is `coordinate >> mip_level`.

This replaces the field-sized stack hierarchy previously rebuilt on every
execution:

| Transform | Old per-call schedule | Immutable compact schedule |
|---|---:|---:|
| GF8 maximum, `N=256` | 224 bytes | 85 bits, stored in 16 bytes |
| GF16 maximum, `N=65536` | 41,760 bytes | 21,845 bits, stored in 2,736 bytes |

The GF16 old size is the exact sum of its five 1,024-word levels, six 16-word
levels, and four final words. The old object was stack-local, so these bytes
were initialized and prepared on every execution. A plan in AUTO mode owns
both generic and specialized metadata because byte size and backend policy may
select either path; forced diagnostic modes own only their selected path.

`test_decode_plan_schedule.cpp` exhaustively enumerates every requested-output
mask for `N <= 16`, then checks deterministic random and sparse masks through
`N=65536` against an independent interval-any predicate. It also rejects bad
sizes, coordinates, storage lengths, mip levels, and traversal parity. The
current deterministic test performs 898,788 dependency queries.

## Input and scatter metadata

The plan caches:

- the generic decoder's selected-input prefix;
- one 16-bit selected-input prefix per `P` or `T` block;
- sorted requested parent coordinates;
- sparse high-rate output block descriptors, including the requested prefix
  and coordinate-list range;
- missing-original scatter indices; and
- the selected recovery index for the `K=1` direct-copy path.

`P` and `T` are at most 32,768 for a nontrivial parent, so a 16-bit local
prefix is exact. In the worst legal `P/T=2`, `N=65536` case, the prefix array
is 65,536 bytes. No-loss, one-parity XOR, one-original copy, and bounded direct
repair plans do not allocate transform schedules.

Only public coordinates that are present and not virtual-erased contribute to
input prefixes. Shortened zeros are never counted, punctures are permanently
erased, and surplus received parity remains excluded by deterministic virtual
selection. Empty blocks are already zero after staging, so planned low/high
execution skips their inverse transforms and byte-heavy reduction/XOR work.

Invocation pointer, padding, scratch, and alias validation remains per call;
those checks depend on caller buffers rather than only on the erasure pattern.
The final output gather walks the cached missing-original list rather than
rescanning all `K` positions.

## Differential and concurrency evidence

The dedicated test compares the complete work arrays from old prepared and new
planned generic, low-rate, and high-rate kernels in both GF8 and GF16. Its
patterns include holes, completely empty input blocks, and sparse outputs in
multiple high-rate blocks. It compares 3,840 work slots.

The same test warms a public low-rate plan, enables a process-wide C++
allocation counter, and observes zero allocations during execution. Eight
threads then share that immutable plan for 128 executions with distinct output
and scratch buffers and reproduce every original byte.

## Pinned scalar measurement

Measurements used GCC 13.3, the portable scalar backend, one OpenMP thread,
64-byte GF16 shards, four missing originals, and CPU 16 (`taskset -c 16`). Old
and new binaries were run in an ABBA order, three repetitions each, with five
timing samples, two warmups, and reuse 256. The table reports the upper middle
of six run medians; the neighboring medians are retained under the ignored
`build/plan-abba/` directory.

| Decoder cell | Old execute (us) | Planned execute (us) | Change | Old setup (us) | Planned setup (us) |
|---|---:|---:|---:|---:|---:|
| generic, `K=1000 R=200` | 293.210 | 292.707 | -0.2% | 193.623 | 194.212 |
| high, `K=1000 R=200` | 233.134 | 229.612 | -1.5% | 193.462 | 195.213 |
| low, `K=32 R=1000` | 116.313 | 15.841 | -86.4% | 193.642 | 193.873 |

The generic result is neutral within run-to-run noise. High gains slightly
from cached prefixes and sparse output descriptors. The low-rate cell has 64
`P=32` blocks but only four parity blocks selected in addition to surviving
systematic data; skipping the remaining zero blocks explains its large,
parameter-specific improvement. This does not imply the same gain for dense
low-rate patterns.

A separate three-run sweep measured execution and setup-amortized latency at
the required reuse counts. Each entry is `execution / (execution + setup/reuse)`
in microseconds; setup medians stayed between 193.7 and 194.9 microseconds.

| Decoder cell | reuse 1 | reuse 8 | reuse 64 | reuse 1024 |
|---|---:|---:|---:|---:|
| generic, `K=1000 R=200` | 291.733 / 485.795 | 293.332 / 317.612 | 291.927 / 294.954 | 290.807 / 290.996 |
| high, `K=1000 R=200` | 228.603 / 423.495 | 229.506 / 253.835 | 229.145 / 232.184 | 229.387 / 229.578 |
| low, `K=32 R=1000` | 15.800 / 209.482 | 16.015 / 40.224 | 15.763 / 18.791 | 15.840 / 16.030 |

The machine-readable runs are ignored build artifacts under
`build/plan-reuse-stable/`. Representative reuse-1024 SHA-256 values are
`a3ab598d934f1da77201aa5f4cf48bd77d77f1147ff35c6553986b4d9eb6ea06`
(generic),
`438ff0a113e9a8fb6db7cc2c0e902c5bc5ce7948a899c63428b6880293c1ca91`
(high), and
`1bd8b0d76656f9ed1394d34e6b3aa2ff0117a53cb0f2dbd37096150d04e32c85`
(low).

Reproduce one side of the paired run with:

    taskset -c 16 env OMP_NUM_THREADS=1 \
      ./build/plan-schedule-dev/bench_leopard2 \
      --k 1000 --r 200 --profile high --field gf16 --bytes 64 --loss 4 \
      --batch 1 --reuse 256 --iterations 5 --warmup 2 --threads 1 --seed 1 \
      --force-specialized --json build/plan-result.json

The validation build and test commands are:

    cmake -S . -B build/plan-schedule -DCMAKE_BUILD_TYPE=Release \
      -DLEO2_BUILD_TESTS=ON -DLEO2_BUILD_BENCHMARKS=ON
    cmake --build build/plan-schedule -j32
    ctest --test-dir build/plan-schedule -j32 --output-on-failure

The forced-backend matrix also runs this schedule test under AUTO, scalar,
SSSE3, and AVX2. All four variants passed with identical normalized output in
the final run. Its source fingerprint was
`052487dffcad285fab43dbe3c1b57b2e2f39ff83dfbeb4f941cec5f1723bc214`;
the merged `matrix.json` SHA-256 was
`b463b1a95cf1f04b4382f0ae24c5cec06b87c7b8d21b8a69ad2f62f44fdbe764`.
The matrix fingerprints `Leopard2Plan.cpp`, `Leopard2Plan.h`, and the schedule
test itself so resume cannot accept an artifact that predates this code.

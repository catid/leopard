# GF16 AVX2 range-table evidence

Commit 7921998 hoists the eight fixed GF16 nibble-table broadcasts out of each
pair in a non-fused radix-four AVX2 transform range. The scheduler contract
already guarantees disjoint shards in the range, so the implementation runs
each original layer across all pairs while preserving the forward or inverse
layer order at every coordinate. It does not change the wire format, public
API, AUTO backend policy, fused schedule, or AVX-512 schedule.

The accepted same-source comparison used the production public encode path at
K=4096, R=512, 4096-byte shards, explicit AVX2, one thread, CPU 12, and a
reserved SMT sibling that recorded zero non-idle jiffies. Three independent
ABBA rounds produced a 1.05645x speedup with a Student-t 95% confidence interval
of [1.05030, 1.06263]. All 12 processes verified public parity byte-for-byte
against both split-output and contiguous-work direct production encodes. A
K=1000, R=200, 64-KiB one-round neighbor screen was 1.00396x and is directional
only.

The measured candidate initially duplicated the prepared butterfly body.
Commit 7921998 factors it through the ordinary butterfly helper. The measured
and factored AVX2 objects have identical text/data/bss sizes, and both prepared
range symbols have identical addresses, sizes, and disassembly. The remaining
ordinary-tail differences are equivalent commutative LEA operand ordering and
NOP encodings. Thus the accepted A/B result transfers specifically to the
unchanged range kernels in the reviewable implementation.

Exact Leopard main comparison remains open. Three complete 36-invocation
campaigns passed parity, round-trip, and workload-digest checks but were
rejected because the reserved SMT sibling accumulated 30, 11, and 9 non-idle
jiffies. The last attempt followed a clean 5.25-second 0/0 non-idle presample;
all 24 temporarily shielded processes were restored with no errors or vanished
PIDs. The failure bundles pass the runner's failure verifier, but none of their
timings may be analyzed or reused. No exact-main speedup claim is made. Bead
leopard-79h.38.5.4 must stay open until the unchanged isolation policy passes on
an exclusive host.

Validation included the effective 83-test Release suite, fresh AUTO/scalar/
SSSE3/AVX2/AVX-512 backend matrices, a 72-test Clang ASan+UBSan build, strict
GCC and Clang warning builds, focused concurrent/API/oracle/GF16 tests, and the
portable-ISA gate. The AVX2 object contains no ZMM, opmask, or YMM16-31 code.

The retained process medians are in raw_same_source_abba.json. The concise
machine-readable implementation, validation, exact-main rejection, digest,
restoration, and reproduction record is in summary.json. Large failure bundles
remain under the ignored .research/leopard2/gf16-main-gap directory.

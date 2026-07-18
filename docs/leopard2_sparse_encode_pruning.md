# Leopard2 sparse parity-output pruning

Leopard2 permits a caller to request any parity subset by passing null output
pointers to `leo2_encode`.  The legacy high encoder and the first production
low encoder trimmed only the unused suffix: a request for parity 0 and parity
R-1 still executed the same forward-transform prefix as a dense request.

An experimental encoder path now compiles an exact forward dependency schedule
into caller-provided encode scratch.  The schedule has a two-bit requested-row
mask per radix-2 butterfly in the ordinary parent-preserving LCH traversal.
This distinguishes x-only, y-only, and two-row operations so a one-row
dependency does not write its dead peer. Compilation walks
that traversal backwards from the exact output mask.  For a forward butterfly

    x' = x + m y
    y' = x + (m + 1) y

the required input rows are

    need(x) = need(x') or need(y')
    need(y) = (need(x') and m != 0) or
              (need(y') and m != 1).

The tests exercise the `m = 0` and `m = 1` special cases as well as ordinary
fixed multipliers.  Complete retained two-layer groups use the selected
backend's fused radix-four operation.  Ragged boundaries use its radix-two
operation, and zero multipliers reduce to XOR.  Low-profile evaluation keeps
the coefficient block immutable: only retained root butterflies are written
out of place, then the remaining exact schedule runs in place.  High-profile
evaluation applies the same schedule directly to its accumulated coefficient
block.

The schedule is compiled once per forced experimental `leo2_encode` call and
reused for both the aligned prefix and a staged final 64-byte tile.  It does not
mutate the codec, does not allocate, and is therefore safe for concurrent
executions that use independent caller scratch.  The worst legal GF16 geometry
adds at most 65,536 bytes of schedule storage: one packed operation-mask array
per parity block and one reusable packed P/T-bit dependency workspace. Mid-size
low-profile geometries whose many parity-coset masks would exceed that bound
fall back to the mature prefix transform pending a compact measured schedule.

Normal AUTO execution retains the mature fused-prefix kernel for dense and
sparse requests and does not pay schedule compilation or extra scratch. Exact
mask execution is currently selected only by the existing test-only
FORCE_TRANSFORM diagnostic. Promotion requires authoritative
pinned, same-source crossover measurements that include schedule setup for
one-shot and reused calls, followed by deterministic size/mask dispatcher cells.
The current evidence is structural and diagnostic, not a throughput claim.

Correctness is checked in three independent ways:

- high/low GF8/GF16 output subsets through the forced path are compared with
  the direct systematic generator-matrix oracle;
- the schedule executor is compared with the full padded LCH transform for
  scalar, SSSE3, and AVX2 backends, shifted cosets, arbitrary byte tails, and
  both in-place and immutable-source forms;
- small transforms are compared directly with polynomial evaluation over the
  legacy GF8 and GF16 representations.

The test hooks count exact blocks, prefix butterfly equivalents, retained
butterfly equivalents, and requested output copies; plan statistics separately
count one-row operations. They are diagnostic instrumentation only and are
absent from the public API.

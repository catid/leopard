# Leopard2 static batch scheduler

The candidate context pool assigns regular encode and decode batch items with
deterministic static ranges instead of a shared atomic fetch-add queue.  A batch
is regular only when every item has the same `shard_bytes` and, for encode, the
same requested parity mask.  A difference in either cost signal keeps the
mature dynamic queue.  The calling thread is participant zero.  Persistent
workers keep the increasing participant slots assigned when they are created.
For `Q` regular items and `W` current participants, participant `i` receives

    base = Q / W
    remainder = Q % W
    begin = i * base + min(i, remainder)
    end = begin + base + (i < remainder ? 1 : 0)

These half-open ranges are contiguous, cover `[0,Q)` exactly once, and differ
by at most one item.  The division form avoids multiplying `Q` by `i`, so its
address-independent scheduling arithmetic remains defined at `SIZE_MAX`.

Each public batch item has its own caller-provided execution scratch.  Batch
preflight rejects cross-item writable overlap before execution.  Consequently,
a participant has exclusive ownership of every item scratch span in its static
range for the duration of the call.  The scheduler requires no per-call queue,
range table, or internal scratch allocation.  Lazy creation and growth of the
persistent worker set are unchanged; once grown, excess workers receive empty
ranges when a later batch is smaller.

This is a scheduling change only.  It does not change profile or field
selection, codec bytes, shard tails, alias validation, error codes, the
lowest-index failure rule, the caller's thread-count budget, or OpenMP's
existing single-thread treatment inside pool participants.  Multi-item calls
that use one context pool remain serialized.  Empty and single-item batches
bypass the pool and may execute concurrently under the ordinary distinct-buffer
rules; a one-thread context has no pool.  Separate contexts and ordinary
one-shot calls remain concurrently usable under their documented buffer rules.

## Current topology boundary

The default thread-count query honors the process affinity mask where supported,
but the pool intentionally does not pin or migrate threads.  Static ranges are
the first production step toward topology-aware scheduling because they provide
stable work and scratch ownership on which a later placement policy can rely.
The following work remains separate and opt-in:

- enumerate allowed CPU identifiers and physical-core/SMT relationships rather
  than retaining only their count;
- expose a caller executor or explicit placement policy without changing the
  default affinity state;
- associate stable worker slots with caller-selected CPUs or NUMA nodes;
- first-touch caller- or library-owned per-worker scratch on its selected node;
- partition stripes by NUMA node and measure cross-node reads and writes;
- refine the conservative heterogeneous-work model if measured output-mask or
  codec-specific weights justify more static regions; and
- define safe fork, context-destruction, and affinity-restoration behavior for
  any future placement layer.

No NUMA, affinity, huge-page, or topology-discovery behavior is enabled by this
milestone.  Such behavior must be separately selectable and validated on each
supported platform before production use.

Static ranges remove shared queue contention but trade away dynamic balancing.
No throughput improvement is claimed by this correctness milestone.  The
production guard therefore selects them only for the exact uniform byte/mask
case.  Heterogeneous byte sizes and mixed encode output masks retain dynamic
fetch-add scheduling.  A broader static cost model still requires isolated
neighbor evidence before promotion.

## Correctness evidence

The public API contract test exhaustively checks the range invariants across
representative task counts, configured scheduling scales through 128
participants, `UINT32_MAX`, and `SIZE_MAX`.  It also checks deterministic replay,
invalid range queries, lazy worker growth, batches both smaller and larger than
the persistent worker set, actual exact-once item execution and worker-slot
ownership, functional parity, concurrent calls sharing a context, and
an application-owned nested thread group, plus deterministic lowest-index
failures.  Sorted, reverse-sorted, alternating, single-heavy and mixed-parity
batches verify that heterogeneous work selects the dynamic exact-once fallback.
The batch-aliasing test warms a multi-thread context and verifies
that repeated scalable encode/decode execution makes no allocations.

## Promotion status

This series remains isolated and is not selected by the integrated production
branch.  A directional 30-thread screen on the development host rejected a
blanket uniform-work promotion.  For GF8 high `(K,R)=(240,16)`, 1 KiB shards
and a batch of 64, static encode took about 19 percent longer while decode was
neutral.  At 64-byte shards and a batch of 256, static encode again took about
24 percent longer, while decode improved by about 11 percent.  A small
low/balanced/GF16 matrix likewise showed path-dependent wins and regressions.

Those measurements are diagnostic, not promotion-grade evidence.  They show
that equal byte counts do not imply equal effective CPU cost on an unpinned,
SMT-capable host.  Before integration, either a narrower deterministic
codec/path/cost gate must pass isolated neighboring ABBA cells, or the dynamic
queue must remain the default and static assignment must be recorded as a
rejected experiment.

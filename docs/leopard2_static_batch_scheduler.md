# Leopard2 static batch scheduler

The production context pool assigns regular encode and decode batch items with
deterministic static ranges instead of a shared atomic fetch-add queue.  The
calling thread is participant zero.  Persistent workers keep the increasing
participant slots assigned when they are created.  For `Q` items and `W`
current participants, participant `i` receives

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
- account for heterogeneous item sizes, for which equal item counts can be
  imbalanced even though ordinary fixed-size stripe batches are regular; and
- define safe fork, context-destruction, and affinity-restoration behavior for
  any future placement layer.

No NUMA, affinity, huge-page, or topology-discovery behavior is enabled by this
milestone.  Such behavior must be separately selectable and validated on each
supported platform before production use.

Static ranges remove shared queue contention but trade away dynamic balancing.
No throughput improvement is claimed by this correctness milestone.  In
particular, batches with heterogeneous `shard_bytes` can assign very different
byte totals to equal-item ranges and may regress.  Performance promotion needs
an isolated heterogeneous-size neighbor matrix or a deterministic policy guard
that retains dynamic scheduling for that regime.

## Correctness evidence

The public API contract test exhaustively checks the range invariants across
representative task counts, configured scheduling scales through 128
participants, `UINT32_MAX`, and `SIZE_MAX`.  It also checks deterministic replay,
invalid range queries, lazy worker growth, batches both smaller and larger than
the persistent worker set, actual exact-once item execution and worker-slot
ownership, functional parity, concurrent calls sharing a context, and
an application-owned nested thread group, plus deterministic lowest-index
failures.  The batch-aliasing test warms a multi-thread context and verifies
that repeated scalable encode/decode execution makes no allocations.

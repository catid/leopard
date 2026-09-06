# GFNI source-staging timing filter

Bead: `leopard-79h.38.5.4.11`. Date: 2026-09-06.
Status at preregistration: **validated front end, not yet timed**.

The correctness milestone `0d0ac30` proves removal of the explicit input-copy
pass at AUTO/GFNI K1000/R200/64 KiB, not a speedup. This new filter asks whether
the existing out-of-place GFNI first stage is faster than the current copy-first
policy. It is not a retry of the exhausted current-versus-Leopard1 `.10` plan.

One new executable handles OFF and ON, with identical code layout and archive.
The only codec call intercepted is the existing source-policy function. There
are no libc-copy wrappers. A bounded 64-record trace remains in both modes;
this is a diagnostic, not a zero-instrumentation production benchmark. The
front end delegates the unchanged existing `current_route_screen.cpp` loop:
one initial encode, four additional warmups and 21 timed public encodes.
It checks all pass arguments and counts after the loop, then prints one trace.
The maximum is 52 records; overflow refuses before delegating to the codec.

The public workload, arithmetic, tiling, output binding, archive and production
sources remain unchanged. Only the previously validated narrow source-policy
predicate can alter a pass. Same-current controls use the identical executable
path, inode and arguments. `--exercise` separately runs 26 clock-free check
workloads to test full trace capacity; it does not model the timing loop's
allocation or warmup behavior.

## Untimed checks

- All 24 Release/sanitizer × OFF/ON × six-cell records match the original
  workload. Twelve Release parity files match every byte of exact Leopard1,
  totaling 122,028,032 bytes.
- Four 26-encode target exercises pass, including both sanitizer modes, with
  exactly 52 matching passes and either zero or 52 changed policies.
- Both 16- and 64-record capacities pass overflow-before-delegation, reset and
  neighbor checks in Release and ASan+UBSan+LSan.
- Six rebuilt default copy probes match all old workload/copy/policy records
  and parity files. Sixteen malformed or invalid-timing requests refuse before
  entering a workload. An inherited `OMP_NUM_THREADS=64`/dynamic setting is
  overridden locally by the driver, without changing another process.
- Seven pure collector tests pass with normal Python and `-O`.

The initial build's final binary-equality assertion failed: changing the guard
from `calls == 16` to `calls >= capacity` changes its branch encoding. Retained
disassembly shows exactly that comparison/jump difference in the wrapper.
All compilation commands succeeded; the failed assertion and original recipe
are retained. The corrected recipe removes that invalid equality expectation;
the six functional comparisons above replace it. The default 16-record limit
and all codec/resource checks remain enforced. Build peak was 88,977,408 bytes under 512 MiB. All successful native
checks stayed under 256 MiB, maximum 141,807,616 bytes, with all six memory-event
counters and swap zero.

## Fixed decision rule

The committed JSON fixes foureyes CPU22/sibling86, a ten-second passive gate,
one attempt, six cells, three OFF/ON ABBA rounds and three same-current ABBA
rounds per cell: 144 timed processes. The target needs at least 5% aggregate
gain and a gain in all three rounds. All six same-current and five mechanically
unchanged OFF/ON controls must lie within the 2% equivalence band. Otherwise
the filter rejects the candidate or is inconclusive, as specified in the plan.
No partial inference, retry, CPU substitution or movement of other workloads.

The filter never claims confidence intervals, a production promotion, an exact
Leopard1 speed comparison or completion of v19 qualification. The experiment
and parent performance-gap tasks remain open. Review is Codex self-review and
deterministic/sanitizer validation under the user's Claude opt-out, not an
independent-model `CONVERGED` result.

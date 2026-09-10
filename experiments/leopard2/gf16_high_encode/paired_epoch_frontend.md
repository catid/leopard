# Three-epoch frontend: completed clock-free qualification

Tracker `leopard-79h.38.5.4.19.1.4.3`, 2026-09-10.
**Clock-free qualification is complete. No real timing, production change,
new performance result, or AUTO R199/32KiB enablement.**

This implements the separately [reviewed diagnostic design](paired_epoch_method.md).
It does not replace the qualified two-endpoint frontend or reuse any consumed
timing screen. The codec archives and their original native-Leopard1 comparator
are unchanged. The three-epoch method remains diagnostic-only even if a future
timed version passes every inherited performance gate.

## Implemented boundary

One codec/buffer lifetime contains three complete epochs. Each epoch records a
before snapshot, performs four actual-route preflights, four warmup schedule
passes and21 exercise passes, then checks accounting/input/probe state and an
after snapshot. All six snapshots must agree with endpoint0. Reference parity
is initialized only by epoch0/slot0, and the external public-call witness never
resets. Exactly252 sample slots are reserved before the first epoch.

The public-encode lambdas and grouped call/timing region remain source-identical
to the pinned original; no claim of machine-code identity is made. GF8 retains
the AUTO request and groups1/256. All nine cells, four L2 schedules, native
comparison, and ordinary one-item batch semantics remain in scope.

Only abort/synthetic clock executables are built. Actual ELF inspection refuses
real steady-clock imports; `--measure` refuses before allocation. The verifier
reconstructs cumulative public-call/state/API/order witnesses, six boundary
marks, per-epoch counts and the intervening preflight/warmup calls at synthetic
endpoints168/336. Arithmetic faults are checked separately in each epoch.

A scalar-only global destructor emits selection/snapshot progress on ordinary
return and `std::exit`, independently of success records. It does not inspect
freed buffers. A captured snapshot whose equality check failed can be included
in the progress count; that count does not assert successful validation. Missing
progress after a signal is incomplete evidence, never an inferred zero count.

## Bugs found during qualification

1. The existing codec diagnostic mode getter maps both normalized raw1 and
   armed raw3 to1. Thus the initial epoch guard could miss a rearmed probe whose
   counter happened to match the last preflight, notably a zero-probe slot.
   The new shared `RequireQuiescentProbe` also requires that the finish-probe
   operation return false. On true raw1 this does not write anything. On leaked
   raw3 it returns true and the guard immediately fails; it never normalizes
   and continues. No codec implementation or diagnostic API was changed.
2. Import-discovery alone could accept a missing retained local module if its
   file and digest entry were removed together. Exact required builder/runner
   tool inventories now supplement transitive import and artifact checks.
3. The first verifier assumed the pinned older build manifest contained a
   `root` field. It does not. The verifier now uses that manifest's independently
   pinned original command root. This error stopped verification before any
   codec execution; its source/log remain in the first run's `precheck-failure`.
4. Checking retained tool digests and the executing main verifier alone did
   not bind its executing dependencies. Final replay now checks every loaded
   local module against the frozen source hashes, including duplicate/aliased
   main modules. Isolated subprocess regressions first pass with matching
   tools, then reject drifted unit/prior/overlay dependencies and an aliased
   main while the retained files remain unchanged.

The user-authorized local read-only reviewer independently confirmed the first
two findings, identified the fourth, and accepted the corrections with no
remaining material finding in its bounded static review. This is Codex review and
deterministic evidence under the user's Claude opt-out, not an independent
Claude fixed-point `CONVERGED` claim.

## Final validation

- All366 matrix processes completed with the expected outcomes:270 positive,
  18 deliberate clock aborts,36 arithmetic-clock refusals and42 CLI refusals.
- Positive records contain1,434,240 public calls,57,240 selections,1,620
  snapshots and22,680 synthetic spans. All261 full native-reference parity
  comparisons pass, covering1,582,128,320 bytes. These are not timed samples.
- Twenty-two overlay/oracle/provenance tests pass, including four real-record
  mutation methods, plus twelve retention tests. Both Python modes pass with
  no skips in the worktree and sealed copy (34 methods per sealed mode).
- 42 supplemental native processes pass: original-native, Release and full
  ASan/UBSan with leak detection. Geometry/probe units have33/45/45 cases.
- The actual shared probe guard is tested against the old predicate's false
  pass, rearmed probes at each of three boundary iterations, repeated normalized
  checks, disabled modes and counter mismatch.
- Native units accept312 selections and six snapshots, reject excess capacity
  and changed later endpoints, and exercise ordered/cumulative mark refusals.
- Actual synthetic-wrapper units accept all504 endpoints, refuse the505th
  before writing, and reject unknown faults and malformed fault epochs.
- Every matrix record has its own validated resource envelope. The highest
  child peak is182,046,720 bytes; controller peak45,821,952 bytes, elapsed
  30:16.55 and exit0. The supplemental unit maximum is17,457,152 bytes.

| Final scope | Peak bytes | Cap |
| --- | ---: | ---: |
| Frontend build / unit build |219262976 /168275968 |512MiB |
| Source tests, normal / optimized |38699008 /46440448 |256MiB |
| Source retention tests, normal / optimized |14860288 /16924672 |256MiB |
| Raw replay, normal / optimized |82419712 /83603456 |256MiB |
| Private retention |72028160 |256MiB |
| Sealed replay, normal / optimized |82890752 /83910656 |256MiB |
| Sealed34-test suite, normal / optimized |37920768 /45608960 |256MiB |

All listed final scopes and native children have zero memory-event counters
and zero swap, with no raised caps. Normal/optimized full replays independently
check the original build/run inventories, ELF/recipe identities, all raw
records/parity bytes, supplemental units and final24-file tool closure.

The first full matrix was deliberately stopped after the probe-guard discovery.
It retains241 records:240 validated records, including178 positive processes,
then the explicitly terminated owned sanitizer child with exit143. This is not
a complete qualification. Nothing from it is labeled timing evidence.

The corrected366-process matrix is terminal at
`/tmp/leopard-paired-epoch-quiescent.qGx924`; the earlier root
`/tmp/leopard-paired-epoch-final.LVkuUB` is preserved unchanged. Each native
child has the canonical lock,256MiB/no swap,60-second CPU and120-second wall
bounds. Builds are serial512MiB/no-swap with compiler GC10/4096. No SSH worker
or unrelated process/host setting is used or changed.

## Retained evidence and limits

The read-only bundle `.research/leopard-79h/paired-epoch-qualified.unzse9`
contains1984 files and2,967,079,659 bytes, including the separately labeled
`stopped_attempt`. Outer manifest SHA256:
`a4749f44c89c6528c1280699f8fce4a22b3b1e917e78be3a72a5a8f4ed09d772`.
All four actual raw/sealed normal/optimized replay results have SHA256
`9bc766dc89c1ac37f1f7257bc7e4028907b958f4c6ac8086abe1e8b683438973`.

A separate stdlib-only audit checks every raw/sealed byte digest, exact sealed
file/directory namespace and read-only modes, distinct source/copy inodes,
both complete and stopped inventories, replay equality and final test/resource
logs. No source or earlier evidence was overwritten. Exit143 alone does not
establish deliberate termination; the recorded owned-process stop supplies
that history, and the incomplete run is never accepted as qualification.

Delivery logs and independent audit source are at
`/tmp/leopard-paired-epoch-delivery.cYar1M`. The
[structured result](results/paired_epoch_frontend_20260910.json) records the
full totals and delivery identities. Final replay tools are a separate pinned
inventory:14 Python files, eight baseline assets and two unit C++ sources.
The original captured build/collection tools are not rewritten.

The implementation stays on `codex/claude-fable-5-1-audit`; no mainline merge
or public release is asserted. Next qualify the separate steady-clock frontend
and collector/replayer, then independently review and commit **and push** the
fresh one-attempt diagnostic preregistration before any real clocks. Even a
future all-pass diagnostic cannot promote R19932 AUTO. Historical control-shift
cause, valid AUTO integration and the user-facing release remain open.

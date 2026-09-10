# Three-epoch frontend: implementation and boundary-check checkpoint

Tracker `leopard-79h.38.5.4.19.1.4.3`, 2026-09-10.
**Full qualification is still in progress. No real timing, production change,
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

The user-authorized local read-only reviewer independently confirmed the first
two findings and reviewed the corrections. This is bounded Codex review and
deterministic evidence under the user's Claude opt-out, not an independent
Claude fixed-point `CONVERGED` claim.

## Validation completed at this checkpoint

- Eleven pure overlay/oracle/provenance tests pass. The four real-record test
  methods are explicitly deferred until the complete new matrix exists.
- 42 supplemental native processes pass: original-native, Release and full
  ASan/UBSan with leak detection. Geometry/probe units have33/45/45 cases.
- The actual shared probe guard is tested against the old predicate's false
  pass, rearmed probes at each of three boundary iterations, repeated normalized
  checks, disabled modes and counter mismatch.
- Native units accept312 selections and six snapshots, reject excess capacity
  and changed later endpoints, and exercise ordered/cumulative mark refusals.
- Actual synthetic-wrapper units accept all504 endpoints, refuse the505th
  before writing, and reject unknown faults and malformed fault epochs.
- Frontend build peak219,262,976 bytes and unit build peak168,275,968 bytes,
  each under512MiB. Pure-test peak14,204,928 bytes and the unit-controller
  peak28,016,640 bytes are under256MiB. These scopes have all memory events and
  swap zero. Individual native units also pass their256MiB/no-swap scopes.

The first full matrix was deliberately stopped after the probe-guard discovery.
It retains241 records:240 validated records, including178 positive processes,
then the explicitly terminated owned sanitizer child with exit143. This is not
a complete qualification. Nothing from it is labeled timing evidence.

The corrected366-process matrix is separately running from
`/tmp/leopard-paired-epoch-quiescent.qGx924`; the earlier root
`/tmp/leopard-paired-epoch-final.LVkuUB` is preserved unchanged. Each native
child has the canonical lock,256MiB/no swap,60-second CPU and120-second wall
bounds. Builds are serial512MiB/no-swap with compiler GC10/4096. No SSH worker
or unrelated process/host setting is used or changed.

The Bead stays open for the full corrected matrix, real-record adversarial
tests, independent replay, evidence sealing and final delivery. A subsequent
steady-clock frontend/collector and committed-and-pushed timing preregistration
remain separate gates; this checkpoint authorizes no real benchmark clocks.

## Pending evidence-retention checkpoint

`retain_paired_epoch.py` and its twelve-test companion are committed as pending
work, not validated delivery. They preserve the complete corrected run and the
separately stopped earlier run in private read-only copies. Static review added
refusals for an existing final-tools manifest and overlapping input histories.
The overlap regression now constructs an otherwise valid nested stopped history,
so the new guard must reject before any destination file is created.

These retention tests have not yet run: the corrected native matrix still owns
the serial qualification workflow. Final-tool inventory/provenance verification,
normal/optimized tests and replays, and sealed-copy delivery remain required.
The implementation and this pending work share the single integration branch
`codex/claude-fable-5-1-audit`; no master merge or release is asserted.

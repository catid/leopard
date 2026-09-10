# Three-epoch plain frontend and analysis qualification

Tracker `leopard-79h.38.5.4.19.1.4.4`, 2026-09-10.
This is an intermediate qualification checkpoint, **not a timing campaign,
performance result, AUTO enablement, or completed collector qualification**.
The parent task and GitHub release goal remain open.

## Native qualification

The separate frontend reuses the [qualified epoch codec archives](paired_epoch_frontend.md).
Only driver/clock/mark objects and their links are new:24 commands produce12
executables across original-native, Release and full ASan/UBSan/LSan profiles.
Every profile has real-steady, plain-abort, witnessed-synthetic and witnessed
fake-steady variants. The real-steady executable runs only `--check` and
`--exercise` during qualification. No real benchmark clock was read.

Fake-steady deliberately advertises `steady` and exercises `--measure` using
the qualified synthetic clock implementation. Its stdout is not real timing
evidence. Executable/source/recipe hashes and ELF bindings distinguish it;
the real measurement paths are separate from all qualification variants.

The common driver captures `__real_*` metadata addresses. Plain links use
linker aliases to the actual public API symbols, not forwarding wrappers.
ELF checks require identical alias/public addresses, and runtime snapshots
are checked against the actual linked image, including PIE load bias.

The terminal189-process run passed with the expected outcomes:

| Outcome | Count |
| --- | ---: |
| Successful native records |126 |
| Deliberate clock aborts |3 |
| Synthetic measured-branch arithmetic faults across all three epochs |36 |
| CLI/mode refusals |24 |
| Full retained native-reference parity comparisons |90 |
| Full parity bytes compared |455,944,704 |

The other36 successful records exercise fake-steady measurement, whose CLI
forbids parity dumps. Those use expected input/output hashes, internal parity,
guards and public-call witnesses; they are not additional full file comparisons.
The old full codec and three-epoch matrices were not rerun.

## Analysis and review

`paired_epoch_analysis.py` and `replay_paired_epoch_analysis.py` independently
project full process records and apply the inherited arithmetic to each epoch.
They preserve318 ordered processes,252 spans per process,300 epoch-specific
controls and264 homogeneous-process epoch ratios. No epoch is dropped,
averaged, selected or pooled. Even all-pass maps only to
`diagnostic_all_gates_pass`; promotion/default/cause flags remain false.

Tests cover a complete synthetic318-file campaign, filename/hash/sibling and
missing/duplicate/reordered-entry errors, symlinks, and an updated-digest
corruption in the last process's last epoch. Separate actual native-record
mutations cover all12 profile/variant combinations. Synthetic fixture output
does not authenticate a campaign or establish measured throughput.

Final replay/test sources are frozen separately from the unchanged captured
build/collection tools. Final verification binds the executing verifier and
every loaded local dependency, including aliases and either analysis module,
to those final source hashes. It checks the exact flat file inventory too.

The user-authorized read-only reviewer found missing file-backed tests, the
executing-dependency binding gap, an omitted extra-directory check, an
unavailable original guard asset in the initial freezer, and an ineffective
overlap-test fixture. Those findings are addressed with focused regressions.
This is static Codex review plus deterministic tests under the Claude opt-out,
not an independent-model `CONVERGED` claim.

## Preserved failures and scope

The first frontend build stopped on its fifth command with an undefined
`__real_leo_encode` link reference, before codec execution. Its original
captured tools, inputs and failure log remain in
`/tmp/leopard-paired-epoch-timing.tYABFS`. The corrected run is separate at
`/tmp/leopard-paired-epoch-timing-alias.jpvPUs`. Initial provenance tests also
retained nine fixture-setup errors for the unavailable guard asset; the final
freezer keeps only the two original source assets actually read by its tests.
Generated header/guard inputs remain bound in the build inventory.

Build peak181,379,072 bytes/512MiB; largest native child181,202,944 bytes/256MiB;
native controller143,273,984 bytes/256MiB. These final scopes exited0 with all
six memory-event counters0 and swap0. Compilation and checks stayed serial
under the canonical lock and existing CPU/wall limits on the local machine.

The final36-test worktree suite passes in normal and optimized Python without
skips (125,165,568 and125,173,760-byte peaks/256MiB, all events/swap0).
This includes11 frontend tests,9 analysis tests,9 final-provenance tests and7
retention tests. The final29-file source inventory is separate from the
original11-file build and15-file collection closures.

The private readonly bundle is
`.research/leopard-79h/paired-epoch-timing-qualified.8oRYqG`:705 files,
671,957,991 bytes, including the original failed build. Its outer manifest is
`b2970bbcc95e0bc9eaf63a96b53bd20039c809187a5930c7426d593871ef4abd`.
Retention peak68,100,096 bytes/256MiB, all events/swap0. Empty source
directories are not copied; the retained regular-file inventory is exact.

Raw and sealed normal/optimized frontend replays all agree, SHA256
`03e8feb5755acf8ef2cc844d9ca003991ff25f33d5f8c3c79d9502cb29d4ba14`.
Final-tool manifest SHA256 is
`f310055a8bdad13604a370b01602133167b36896681cd9db1e16ca16854184a0`.
Delivery logs are at `/tmp/leopard-paired-epoch-timing-delivery.6fvOKU`.
The sealed36-test suite also passes in both modes without skips. A separate
stdlib-only audit checked all705 retained files, raw/sealed byte equality,
private-copy inodes, exact readonly namespace and final resource/test logs;
both modes agree. Audit peaks37,085,184/40,792,064 bytes/256MiB, all
events/swap0. The source and sealed replay/test scopes likewise all pass.
The [structured checkpoint](results/paired_epoch_timing_20260910.json)
records their exact resource peaks and identities.

## Remaining work

The analysis helpers are **not full campaign authenticators**. Implement and
qualify the collector, immutable-input freeze and independent full replay
that bind these helpers to actual plain-steady executable identities,
27 preflights, raw exit/stderr/resource evidence, stopped-service conditions,
isolation and exclusive one-attempt consumption. These must be tested without
real clocks first. A separate reviewed, committed **and pushed** exact
preregistration is still required before the diagnostic can run.

The diagnostic remains non-promotional even if it passes all inherited gates.
R199/32KiB AUTO remains OFF pending separate valid integration evidence.

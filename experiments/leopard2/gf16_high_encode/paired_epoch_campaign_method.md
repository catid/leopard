# Three-epoch diagnostic collector qualification

Tracker `leopard-79h.38.5.4.19.1.4.4`. This is a readiness-false implementation,
not a timing preregistration, performance result, or permission to enable AUTO.
The R199/32-KiB route remains off. No codec archives are rebuilt here.

The campaign reuses the completed `b64efb7` plain-frontend qualification.
Its manifest, build inventory, replay and final tool closure are pinned; only
the actual measurement/native/steady and measurement/release/steady images
are accepted. Synthetic and witnessed executables cannot be substituted.
New sources, transitive old dependencies, plan, condition script, proof files
and private executable copies form an exact immutable flat inventory.

Qualification uses a separate immutable `qualification` bundle with readiness
false and runs only 27 `--check` preflights. It never calls the timing entry
point, passive observation, or `--measure`. Full synthetic campaign fixtures
are tests, not timing evidence. Both Python modes and a read-only code review
must pass before a later, separately reviewed and pushed readiness-true plan.
The timing freezer and collector also independently replay the actual readonly
27-check qualification and its exact resource log. Every input pin except the
readiness-different plan must match. Missing or writable qualification evidence
refuses timing before publication checks or attempt creation.

Both immutable stages share exactly one fixed `attempt1` path. The timing
entry point refuses false readiness before any subprocess, passive clock,
lock acquisition or attempt creation. A consumed attempt cannot be resumed,
replaced or made reusable by another freeze. Each child launch has a flushed,
fsynced intent before execution and a separate durable outcome afterward;
failed/timeout output remains evidence, never a successful measurement.

The campaign runner is the sole owner of the canonical lock and CPU lease.
The scope wrapper must invoke the runner directly; callers must not add an
outer `flock` around it. This prevents a self-deadlock from being mistaken for
an experiment result.

The controller command is also part of the evidence identity. Launchers must
use the absolute `/usr/bin/timeout`, `/usr/bin/prlimit`, and `/usr/bin/python3`
paths emitted by `controller_command`, wrapped by `/usr/bin/time -v` so the
resource footer includes its exit-status record. Equivalent bare command names
or an omitted timing footer are not accepted by the independent verifier.

The unchanged diagnostic consists of nine cells, 27 preflights, 318 measured
processes, three complete epochs per process, 21 sample passes, four warmups,
252 spans, 100 controls per epoch and 264 homogeneous process trajectories.
GF8 uses AUTO and groups of 256. No epoch is selected, pooled or averaged.
Target/native gain gates remain 5%; equivalence bounds remain symmetric 2%.
Even all-pass is diagnostic only and cannot authorize production promotion.

Local host work, CPU 26/sibling 90/controller 0; one canonical lock and CPU
lease, a 256-MiB/no-swap scope, child CPU 30 seconds and wall 60 seconds.
Read-only service/container checks must match before and after. The independent
replayer reconstructs every launch, environment, outcome, preflight, full
metadata record, passive/sibling condition, resource envelope and file name
before deriving per-epoch results. It does not import the campaign collector.
The scope wrapper records the exact cgroup and controller argv, so an unrelated
successful job's resource log cannot satisfy this gate. Recorded logical output
paths allow byte-identical evidence to be replayed after retention elsewhere.

On September 13 the previous test handle and `/tmp` preparation/log directories
were missing following the September 11 reboot. That test's outcome is unknown,
not passing. The original immutable `.3` and `b64efb7` frontend evidence survives.
Before any campaign freeze or timing preregistration, the unarmed plan was moved
to a durable `.research` parent and updated from kernel `6.8.0-137-generic` to
the observed `6.8.0-139-generic`. No consumed timing attempt is reused, and the
earlier performance results are not relabeled measurements of this environment.

The first actual clock-free collector qualification stopped in its first service
check (`rg` absent from the restricted PATH), before any codec launch. Its
immutable `qualification` bundle, incomplete `qualification-checks` journal and
resource log remain untouched under the same campaign parent. The three exact
line checks now use system `grep -Fxq` with the same expected service output.
The corrected qualification uses distinct `qualification-v2` inputs and
`qualification-checks-v2` output; the sole future timing path remains `attempt1`.

No Claude review is used under the explicit opt-out. Deterministic tests and
the authorized read-only Codex reviewer are not an independent-model
`CONVERGED` gate. Collector acceptance and the public release goal stay open
until their own required evidence is complete.

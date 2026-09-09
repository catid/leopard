# Current routes after the user-requested Slipgate shutdown

Bead: `leopard-79h.38.5.4.10`. Date: 2026-09-09. Status: preregistration.

The user explicitly stopped Slipgate Forge and OBS, disabled their boot/login
startup, and stopped their dedicated PostgreSQL/Redis containers. The new host
condition, tracked by closed `leopard-0ge` and commit `cb837e9`, is verified
before and after this one local diagnostic. This does not identify which task
caused the prior interference, prove kernel isolation, or authorize moving other
threads. The servers remain out of scope.

This successor uses `current_route_screen_post_slipgate_plan.json`: same local
9980X `work` CPU26/sibling90/controller0, unchanged codec archives and C++ driver,
six workloads, three ABBA rounds with identical-path controls,21 samples per
child,144 timed invocations,10-second passive zero-sibling gate, and one attempt.
The later passive readiness survey at00:52:13–00:52:23 UTC showed both CPUs idle;
it was not a codec benchmark. No CPU substitution or further host changes occur.

All control aggregates must pass the unchanged2-percent equivalence gate before
any directional inference. Any sibling work stops the attempt, with no partial
analysis or retry. Neither the24 samples from exhausted local `b8b9b1a` nor
server results are reused. The collector distinguishes the post-shutdown plan
from the old local plan and retains all existing protocol/type/identity checks.
The plan and collector must be committed and pushed before launch.

Artifact and correctness provenance remain those documented in
`current_route_screen.md` and the original read-only158-entry bundle
`.research/leopard-79h/gf16-current-route-failed.STc10h`, outer manifest
`e83a217e680a23088e53208d9ba7747742e7918a1377e74efe9c35399b4ba342`.
Recheck all original file hashes and full parity evidence, then make fresh
lane-owned read-only copies. The12 fresh untimed checks must match exact route,
source, scratch and output identities before encode clocks are permitted.

No production code, routing threshold or kernel changes are included. This
diagnostic is separate from unfinished v19/production qualification. Review is
Codex self-review plus deterministic checks under the user's Claude opt-out,
not independent-model `CONVERGED`.

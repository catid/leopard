# Current routes after the user-requested Slipgate shutdown

Bead: `leopard-79h.38.5.4.10`. Date: 2026-09-09. Status: complete diagnostic.

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

## Results

Preregistration `23bcdd7` was pushed before launch. All12 untimed checks and
144 timed invocations passed. Sibling90 stayed at570054 non-idle jiffies during
the10.000241142-second passive window, and every timed invocation had delta0.
The shutdown snapshots before and after matched. All six aggregate controls
passed the unchanged2-percent equivalence gate.

Ratios below are Leopard2 throughput divided by Leopard1 throughput, computed
from the preregistered within-cell ABBA comparisons; greater than1 favors L2.
They are diagnostic results on this9980X, not confidence intervals or claims
about other CPUs, layouts, APIs or shard sizes.

| K/R/shard bytes | Requested / actual L2 route | L2/L1 throughput | Same-L2 control | Classification |
| --- | --- | ---: | ---: | --- |
| 1000/200/65536 | AUTO / GFNI | 1.415585 | 1.006356 | L2 advantage |
| 1000/200/65536 | AVX2 / AVX2 | 0.969768 | 1.000609 | Investigate L2 deficit |
| 1000/200/65536 | AVX512 / AVX512 | 1.047702 | 0.995393 | L2 advantage |
| 1000/200/32768 | AUTO / AVX2 | 0.971964 | 0.996557 | Investigate L2 deficit |
| 1000/199/65536 | AUTO / AVX2 | 0.956348 | 0.996111 | Investigate L2 deficit |
| 4096/512/4096 | AUTO / AVX2 | 1.009573 | 0.999681 | Near parity / uncertain |

Leopard1 uses its pinned native compiler policy, whereas an explicitly requested
L2 AVX2 route is ISA-restricted. That comparison is not ISA-matched and cannot
by itself attribute the deficit to API or dispatch overhead.

The standalone `replay_current_route_post_slipgate.py` rehashed all16 frozen
inputs, checked every raw stdout/stderr record against the journal and exact
workload identities, and rebuilt all36 round ratios and12 aggregates without
importing the collector or executing any codec. Both normal and optimized
Python runs agreed with the collector; each rejected15 mutations and verified
that a synthetic bad control suppresses every directional classification.

The complete native scope exited0 after73.52 seconds, with memory peak
133,763,072 bytes under268,435,456; all six memory event counters and swap were
zero. No archive, codec kernel, source-routing predicate or production flag was
changed. Earlier failed samples remain excluded. The diagnostic's attempt1/1
is consumed; this is not an invitation to repeat it.

## Evidence-driven next work

The old broad assumption that the current AUTO64-KiB target trails Leopard1 is
not supported: it is41.6% ahead in this diagnostic. The remaining measured
deficits are AVX2 cases, including two default-route boundaries.

Source and history inspection explain why these neighbors stay on AVX2:
`CodecMayUseAutoGF16GFNIEncode` restricts K/R to1000/200 and
`UseAutoGF16GFNIEncode` requires exactly64KiB. The original `8da10dd` selector
promotion qualified one exact cell and inactive neighbors; it did not establish
that selecting GFNI for those neighbors would be slower. This is a routing
hypothesis to test, not permission to widen production immediately.

- `leopard-79h.38.5.4.17` (P0): validate and measure GFNI at the32-KiB andR199
  default-route deficit boundaries, then consider a tightly scoped AUTO change
  only if the full correctness and performance gates pass.
- `leopard-79h.38.5.4.18` (P1): attribute the explicit-AVX2 deficit separately;
  an explicitly requested backend must not be silently changed to GFNI.
- Combined first-inverse/terminal GFNI fusion remains an unpromoted later
  experiment, now lower priority than closing observed default-route deficits.

The full objective remains open, as do broader production/v19 qualification
requirements. This milestone identifies where to improve; it does not claim a
new production optimization has shipped.

Raw evidence and validation are retained locally in
`.research/leopard-79h/gf16-current-post-slipgate.f4j4of`, referencing the
unchanged original correctness bundle above. Structured results are in
`results/current_route_post_slipgate_20260909.json`.

The sealed bundle has351 manifest entries; its outer `SHA256SUMS` hash is
`429564bd1eef3d66b219f991b78f6b02951240bacb71e94761d1285926f435e8`.
The frozen-copy replay passed locally, and the original evidence is retained
without requiring either SSH host. Retention also stayed under256MiB with all
memory events and swap zero.

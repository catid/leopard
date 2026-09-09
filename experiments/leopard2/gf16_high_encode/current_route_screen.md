# Current-route GF16 encoder diagnostic

Bead: `leopard-79h.38.5.4.10`. Date: 2026-09-06.
Status: **historical attempts exhausted; a distinct post-Slipgate successor now has valid diagnostic results**.

See [the 2026-09-09 post-shutdown result](current_route_post_slipgate.md) for
the complete144-invocation local comparison. The failures documented below
remain unchanged and are excluded from that result.

## Local-only successor (2026-09-08)

The user redirected all work to this machine because the SSH server is busy.
No further inspection, computation or evidence transfer on either SSH server is
authorized. The earlier proposal to move unrelated server threads is superseded,
not approved. No unrelated local process will be stopped or have its affinity
changed, and no persistent host settings will change.

`current_route_screen_work_plan.json` separately preregisters one local attempt
on `work`, Threadripper 9980X, kernel6.8.0-137-generic, CPU26/sibling90 with
controllerCPU0. A read-only ten-second CPU survey showed both idle; no codec
timings were used for CPU selection. The original plan and failure remain
unchanged. The collector accepts only the two explicit host/CPU profiles, not
arbitrary CPU substitution. Pinning applies only to our controller and children.

Reuse the exact locally built binaries and archives listed below, in fresh
lane-owned read-only copies. Replaying the retained correctness evidence again
verified all13 file pins, six sanitizer records and61,014,016 full parity bytes.
All six workloads,144 timed children,3 ABBA rounds,21 samples, identical-path
controls,2-percent thresholds,10-second passive and zero-sibling gates remain
unchanged. No pooling with server observations, retries, production promotion
or v19 closure is permitted. Commit and push this successor before launch.

Preregistration `b8b9b1a` was pushed before the sole local launch. All12 fresh
untimed checks passed, and sibling90 remained at569301 non-idle jiffies during
the10.000076453-second passive window. The run stopped on its24th timed
invocation (cell0, round2, final same-current control): sibling90 accumulated
4 non-idle jiffies. The preceding23 invocations had zero sibling deltas.
The journal has `complete:false` and no analysis. No partial target ratio,
cross-host comparison or performance conclusion is valid. Budget1/1 is
exhausted; there was no retry or CPU substitution.

The local scope exited1 after24.47 seconds, peaked at132,247,552 bytes under
256MiB, and recorded all six memory event counters and swap as zero. No
unrelated workload was moved or stopped. Neither SSH server was contacted.
Raw evidence is `/tmp/leopard-current-local.BUcGF2`. The standalone
`replay_current_route_work_failure.py` verifies all14 frozen inputs,12 raw
checks,24 raw timing records, ordering, route/source identities, exact failure
and resource evidence without computing a ratio. Normal and optimized Python
replays also reject11 semantic mutations each. The original archive's complete
158-entry manifest and full parity replay were rechecked locally before launch.

Next measurement requires a controlled **local** CPU window and a separately
preregistered successor. This user instruction changes the host, not the
correctness, resource or isolation gates. Combined-fusion timing is still
pending; the busy servers are not a fallback.

The read-only102-entry local bundle is
`.research/leopard-79h/gf16-current-route-work-failed.xhiwvW`, with outer
`SHA256SUMS` hash `847a0721b8939b1c03a88cf3afc753e60bf1846643d0a7692aa1479fe5cdfbc6`.
It contains the frozen artifacts, raw records, tests/replays and logs. The
unchanged original correctness archive below is a required separate reference,
not duplicated. The sealed-copy replay passed locally; no second-host copy was
made. Structured details are in `results/current_route_work_20260908.json`.

## Question and scope

The historical 0.976778x AVX-512 result on a 9950X3D does not establish a
current AUTO deficit on a Threadripper model08h: the exact K1000/R200/64-KiB
AUTO operation now selects GFNI. The prior intra-butterfly cache-block
candidate did not affect that path and was rejected in a clean same-source
screen. This separate experiment compares current operation routes with
standalone Leopard1 to identify useful optimization targets.

`current_route_screen_plan.json` fixes six workloads, three ABBA rounds per
workload and three interleaved same-binary control rounds, 21 samples per
process, foureyes CPU22/sibling86, a 10-second passive gate, and one attempt.
It must be committed and pushed before timing. No old samples are pooled.
The two control labels invoke the identical current executable path with
identical arguments. Controls must all fall inside `[1/1.02,1.02]` before
any directional classification. No confidence intervals or promotion claims
are produced. The existing `.8`/v19 qualification and runtime/build-handoff
requirements remain separate, incomplete, and unchanged.

## Source and workload identity

Current codec source is `36dc0c8f66604b8d974468e687c51e6183ecb61d`; the
fresh production Release build includes both fields and default backends,
with tests/hooks/benchmarks/CUDA disabled. Leopard1 is separately built from
clean detached `6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198`, using the existing
standalone adapter's native policy: `-march=native -Wall -Wextra -fopenmp
-g -O0 -O3`, without `-DNDEBUG`. Its archive contains only four original
Leopard1 translation units. It cannot share current Leopard2 kernels.

Both archives were built locally with GCC13.3. The local 9980X and remote
9985WX compiler `-march=native -###` expansions compare byte-for-byte,
including target features, tuning and cache parameters (SHA
`35ca2ea962b79f83eca0ffe19383864ee98d841fa385d6323109b57457e0afce`).
The screen is not an ISA-matched AVX2 comparison: native Leopard1 may use the
compiler's native instruction choices, while current routing is recorded
explicitly.

The new driver is compiled independently for the two APIs. Both use the
same aligned XorShift32 source with seed20260906 and full recovery output.
Every sample includes one public encode, but excludes initialization, codec
creation, allocations, source generation and hashing. Leopard1's result
remains in its first R work buffers; no extra copy is imposed on it.
Leopard2 uses its separate public output buffers and scratch. There are no
decode, batching, or plan-setup amortization measurements.

## Untimed validation

All six local pairs produce byte-identical full parity files, not merely
matching digests. Check mode reads no clock. Initial measured-route
expectations are:

| Cell | K/R/bytes | Current requested/selected route | Current scratch bytes | Leopard1 work bytes |
| --- | --- | --- | ---: | ---: |
| 0 | 1000/200/65536 | AUTO / GFNI | 16,808,512 | 33,554,432 |
| 1 | 1000/200/65536 | AVX2 / AVX2 | 16,808,512 | 33,554,432 |
| 2 | 1000/200/65536 | AVX512 / AVX512 | 33,585,728 | 33,554,432 |
| 3 | 1000/200/32768 | AUTO / AVX2 | 16,808,512 | 16,777,216 |
| 4 | 1000/199/65536 | AUTO / AVX2 | 16,808,512 | 33,554,432 |
| 5 | 4096/512/4096 | AUTO / AVX2 | 4,308,992 | 4,194,304 |

The Release builds peaked at 383,315,968 bytes (current) and 177,393,664
bytes (main) under512MiB. Full parity checks, including output file cache,
peaked at229,519,360 bytes under256MiB. All six memory event counters and
swap were zero in these jobs. The first five pure collector tests passed.

Two full dual-field GCC sanitizer builds, at `-O1` and then `-O0`, exceeded
the512MiB scope while compiling the heavily force-inlined GF8 T16/Q2 unit.
Scopes `run-u283966` and `run-u283994` report `oom-kill`; final memory peaks
and counters did not survive and are not inferred. No field or sanitizer
check was disabled. A bounded single-file test at the original `-O1` with
compiler-only garbage-collection parameters `ggc-min-expand=10` and
`ggc-min-heapsize=4096` passed at521,043,968 bytes, zero events/swap.
The full dual-field `-O1` ASan+UBSan build with that GC policy then passed at
536,637,440 bytes under536,870,912, with zero memory events and zero swap.
This is very little headroom, not a general guarantee for other compilations.
All six sanitizer check-mode records exactly match Release; leak detection
and fail-on-error were enabled, with8MiB quarantine and64KiB thread-local
quarantine. Test peak was140,664,832 bytes under256MiB, zero events/swap.
No field, route or sanitizer check was removed; Release binaries are unchanged.
The five new pure collector tests pass normally and with Python `-O`.

Scratch evidence is `/tmp/leopard-gf16-current-routes.BMj72w`. Review uses
Codex self-review and deterministic checks under the user's Claude opt-out;
there is no independent-model `CONVERGED` claim.

## Terminal attempt and retained evidence

Preregistration `7887d56` was pushed before the one allowed server attempt.
All 12 untimed server checks matched the fixed local expectations. During
the immediately following 10.000058885-second passive window, sibling86's
non-idle counter rose from191252 to191260. The collector exited with
`ValueError: passive sibling activity; attempt stopped`, `complete:false`,
zero timed invocations and no analysis. No ratio, performance direction or
gap-closure conclusion can be inferred. The budget is exhausted: no retry,
CPU substitution, pooled partial data or relaxed threshold under this plan.

Scope `run-rda0aad05282d4aa89af2dfd3b7e3cf2e.scope` exited1 after12.87s,
peaked at132,370,432 bytes under256MiB, and recorded all six memory event
counters and swap as zero. No other workload was stopped or moved. A request
to temporarily exclude unrelated user-space threads from CPUs22/86, then
restore their affinities, remains unapproved; no such intervention occurred.

The separately written `replay_current_route_failure.py` imports no collector
and executes no codec. Normal and optimized Python runs both verified all13
frozen file pins, the raw12 server checks, the six sanitizer records, and
all61,014,016 parity bytes across six local old/new pairs. Peaks were8,445,952
and11,898,880 bytes under256MiB with all events and swap zero. An initial
launcher permission error is retained; it executed no replay or codec.

The158-file read-only bundle is retained locally and on ripper at
`.research/leopard-79h/gf16-current-route-failed.STc10h`. It includes full
source dependencies, build recipes and metadata, frozen Release artifacts,
sanitizer artifacts, raw parity/check files, both failed sanitizer builds'
logs, the terminal journal, host evidence, and the separate replay/logs.

- Outer `SHA256SUMS`: `e83a217e680a23088e53208d9ba7747742e7918a1377e74efe9c35399b4ba342`
- Attempt journal: `f2089e56a1e787bc0dd41ce7fd3c6ec6b54bc89d0b1ef759d403f060fbf4af51`
- Scope log: `78510b95f0b065afc906813e65d3f2e89dcc9b1cd00a14ede5dcc0b69977c3f2`

Existing v19 and K65 campaign gates remain unchanged. The current-route
comparison and the overall Leopard1 performance-gap objective remain open.

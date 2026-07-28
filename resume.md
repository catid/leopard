# Session resume state — 2026-07-28, branch `codex/leopard2`

Beads is the durable task source. Use the Beads 1.x binary explicitly:

    /home/catid/.nvm/versions/node/v24.1.0/bin/bd

The legacy `~/.local/bin/bd` is 0.47 and must not touch this checkout.

## P0 pre-resume bug gate completed

`leopard-79h.49` and all fifteen children are closed. The reviewable code
milestones are:

- `a9e0e3e` — defined `NextPow2` edge behavior, overflow-safe GF16 live-set
  policy, warning-clean reduced-field builds, and backend-correct scalar
  validation;
- `a8881a2` — fail-closed lab process/TID/affinity containment and real runtime
  observation evidence;
- `cf96061` — build-time benchmark source attestation, exact current/historical
  compile provenance, and shadow/nested/superproject defenses.

Final validation at those commits:

- GNU 13 Release: 104/104 CTests passed twice;
- Clang 18 ASan+UBSan: 36/36 focused production tests passed;
- portable strict dual-field build: 97/97;
- strict GF8-only and GF16-only: 7/7 each, warning-free;
- forced-scalar high-pruned regression: 1/1;
- run_abba: 56/56; Visual Studio/CMake graph: 134/134;
- lab, benchmark-matrix, balanced-promotion, attestation-refresh, CUDA-optional,
  ISA-isolation, and field-option self-tests passed;
- three independent review partitions and a final repeat found no unresolved
  actionable defect.

Disposable validation builds are under:

- `/tmp/leopard-bughunt-final-release-1785263016`
- `/tmp/leopard2-bughunt-asan-ubsan-1871043`
- `/tmp/leopard2-bughunt-matrix-1876619`

## Resume production optimization

Continue `leopard-79h.46` first: derive the GF16 live-set target and tiling
threshold from the cache available to the context. Existing measurements show
that the fixed 32-MiB-CCD policy loses roughly 3–9% on this host's 96-MiB
V-cache CCD for 34–67 MiB working sets, while tiling remains beneficial once
the set exceeds the larger cache.

Design constraint: do not scale the target by thread count; that variant was
measured and rejected. A safe implementation should use the minimum relevant
last-level cache across the context's allowed CPUs, retaining the calibrated
32-MiB fallback when topology is unavailable or spans asymmetric cache
domains. A process pinned wholly to the 96-MiB CCD should be able to select the
larger policy. Add deterministic detection/policy tests before timing.

After `.46`, the current priority order remains:

1. `leopard-79h.43` exact-main campaign — requires a bare terminal with no
   Codex/Claude session alive because the fail-closed affinity supervisor sees
   the session's process churn;
2. `leopard-79h.44` explicit-GFNI campaign and possible AUTO promotion;
3. `leopard-79h.45` remaining GF16 table/footprint disposition if still open;
4. `leopard-79h.47` second-microarchitecture validation;
5. `leopard-79h.48` final verified Leopard1 comparison documentation.

Never use screen measurements in the published Leopard1 comparison. Only a
verified immutable `run_abba.py` evidence bundle is authoritative.

## Host and coordination

The allowed CPU set currently contains 30 logical CPUs:
`0-14,16-30` (15 SMT pairs), one NUMA node. Cap parallel build/test work at
30. Serialize builds, correctness runs, and timings with:

    /tmp/leopard-gf8-authoritative.lock

Authoritative timing additionally requires frozen lane-owned binaries and the
provenance checks in `AGENTS.md`. The exact-main campaign remains a bare-
terminal-only task on this host.

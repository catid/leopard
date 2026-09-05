# V19 build relocation feasibility — 2026-09-05

Fresh, unprivileged scratch builds can reproduce all four pinned v19
executables/archives without overwriting the historical build tree. This is a
build-feasibility result for the frozen source identities below, **not a
performance result or permission to acquire v19 measurements**.

The candidate is `cf7a7056e0bd7f54b8da436a39cae857beab10c1`; Leopard1 is
`6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198`. The original preflight and its
preregistration were not changed. The compact machine-readable evidence is
[the relocation checkpoint](results/v19_build_relocation_feasibility_20260905.json).

## Why the original path mattered

The exact-main adapter intentionally reproduces Leopard1's Release flags
`-g -O0 -O3`, including debug information. Native inspection found the original
source and build paths in the pinned ELF's DWARF sections. Changing a scratch
directory therefore changes the full-file identity even if execution code
would be equivalent. Code-section-only equality is not the pinned contract.

Private namespace shadowing was unavailable: the unprivileged `unshare`
capability probe failed while writing `/proc/self/uid_map`. No privilege,
security-policy change, mount, or historical-tree replacement was attempted.

The successful alternative adds one explicit flag to the baseline build:

```text
-ffile-prefix-map=<physical-workspace>=/home/catid/leopard-v19-build-preflight.QSq8egne
```

Both workspaces retain the child names `candidate-source`, `leopard1-source`,
`candidate-build`, and `baseline-build`. The flag is supplied through
`CMAKE_CXX_FLAGS`; the original Release/native flags remain in effect. The
candidate build uses the original preflight recipe with no additional flag.
Actual physical paths remain in the diagnostic's commands and build metadata;
the mapped debug paths must not be substituted for source-ownership evidence.

## Contrasts and resource result

| Experiment | Result |
| --- | --- |
| Mapped baseline, first scratch path | Both full-file hashes equal the pinned baseline hashes |
| Unmapped baseline control, different scratch path | Both hashes differ, as expected |
| First full fresh-source build | Candidate build completed; post-build memory-event gate rejected the run before baseline build |
| Full fresh-source build with owned-file cache eviction | All four full-file hashes match; host, lock, retained-preflight and resource checks pass |

The full builds use fresh depth-one detached candidate and baseline clones and
the pinned `sse2neon` commit. The earlier baseline-only cases copy only two
adapter files, checking their hashes against the retained candidate source
manifest; they are not full candidate-source builds.

The successful full probe ran on ripper under the actual canonical lock,
one-job builds, 512 MiB cgroup limit, no swap, and a 65,536 soft open-file limit.
It took 32.61 seconds for the whole build experiment. Its cgroup peak was
521,965,568 bytes (497.8 MiB), with all six memory-event counters zero. Maximum
child RSS was 387,120 KiB; that is **not** the total cgroup peak.

Before compilation, the probe flushed and requested cache eviction for 1,341
independently cloned regular files totaling 114,767,874 bytes, using `fsync`
and `POSIX_FADV_DONTNEED`. It rejected aliased files and never requested
eviction for original source/preflight or shared Git objects, or used global
cache controls.
Observed cgroup usage around that operation fell from 175,935,488 to
52,125,696 bytes. This is an advisory operation with measured endpoints, not a
general guarantee about cache residency.

The rejected first full run remains excluded. Its helper did not preserve the
exact rejected cgroup counters before scope teardown; no particular counter
value is reconstructed from its error message. The final helper records raw
cgroup endpoints on exit. All benchmark executables remained unexecuted.

## Evidence and integration boundary

The retained diagnostic bundle on ripper is:

```text
/home/catid/leopard/.research/leopard-79h/v19-build-relocation-feasibility.NS1qbI
```

It contains the four experiment records, their scripts and command logs, the
owner modules, final build caches/compile commands/link metadata, and copies
of all four matching outputs. All 80 manifest-listed files passed checksum
verification. Files are mode 0444 and directories 0500. The outer
`SHA256SUMS` digest is
`16edea6a9fcc5c66c5afdc401c1c38cf961186ae4ff69cf2426b05c051ab85b1`.
The checkpoint records individual result, final-probe and output hashes.

For reproduction, restore the bundle's `owners` directory as `repo` beneath a
fresh private diagnostic root, alongside the chosen experiment's `probe.py`
and `preregistration.json`. The final probe accepts `full`, `mapped`, or
`unmapped`; it requires ripper's native preflight conditions and acquires the
canonical lock itself. Its immutable retained preflight inputs must still
exist at their original location. Do not wrap it in a second canonical lock,
run it on occupied scratch output, or treat its result as campaign authority.

The current exact-main compile-command comparison requires exact argument
lists. Byte equality does not silently authorize the additional mapping flag.
A v19-only recipe/closure change must explicitly bind the actual workspace,
original debug-path root, exact mapping argument, source identities, and four
output hashes without relaxing v1–v18 replay. Full source/runtime lifetime
ownership, physical v18 archive verification, failure sealing and acquisition
wiring are still required under Beads `leopard-79h.38.5.4.8.2.2.2` and its
parent. The observed 14.2 MiB memory headroom also means that a future source
owner must avoid retaining whole source-file contents during compilation.

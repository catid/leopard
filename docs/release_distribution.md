# Release distribution

The public source release is produced with `git archive`.  It contains the
library sources and headers, CMake configuration, tests, examples, benchmark
programs, portability support, license, and user-facing documentation.

Repository-only material is retained in Git history but marked
`export-ignore`: `.research/`, experimental campaigns other than the small
`experiments/leopard2/performance_atlas/` reproducibility tooling, Beads/session
metadata, agent instructions, and working research notes.  The excluded files
are not required to build or use Leopard2 and may contain machine-local paths
or large evidence bundles.

The archive must be made from the release commit, not from a build directory:

```sh
git archive --format=tar.gz --prefix=leopard-2/ HEAD > leopard-2.tar.gz
```

After extraction, a clean release build uses only the documented CMake
workflow.  No `/home/catid`, `.research`, temporary build, or private toolchain
path is a release dependency.  Performance evidence and reproducibility notes
that are useful to users live under `docs/performance/` and remain included.

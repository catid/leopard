# Low encoder coefficient-copy A/B evidence

`run_abba.py` is the authoritative evidence collector for the LOW_V1 encoder
change between control commit `4070e4e527935026fb87593567587558f0a08d51`
and candidate commit `6d3afee213b94d486cf5f8145ac18078883ebc20`.
It does not build either side.  It accepts only the clean, tests-disabled,
Release production builds described in the reproduction guide, verifies their
source/object/archive/link/runtime closure, cleanly recompiles every retained
source/object pair byte-for-byte, and then runs the fixed matrix. New evidence
uses raw/manifest schema v4 and failure schema v4, binding the canonical CMake
target `leopard`, archive `libleopard.a`, and `leopard.dir` object closure.
V4 retains the exact bounded UTF-8 archive `link.txt` content, binds its byte
size and SHA-256 to the existing recipe-file identity, and parses those bytes
to require the declared archive and ordinary-object target directory (backend
object-library directories remain distinct allowed targets). Retained
schema-v3 evidence replays against its exact historical
`libleopard`/`liblibleopard.a` identity and record shape. A schema/path relabel
that retains the old recipe bytes is rejected.

The evidence digests are unkeyed integrity checks, not an external
authenticity anchor. They detect edits relative to the retained evidence and
the v4 semantic binding prevents re-signing only the schema and path labels
around old recipe bytes. They cannot prevent an attacker from replacing every
evidence byte and recomputing a completely internally consistent bundle.
Preventing that stronger evidence re-authoring requires a trusted external
signature, transparency log, or equivalent authenticity anchor.

The frozen control and candidate commits both predate canonical-target commit
`0a2a666`. They cannot produce a fresh v4 bundle while also satisfying the
collector's exact-clean-commit rule: both source trees define only the
historical target. Do not work around this with an unrecorded generated-build
edit, copied archive, or dirty CMake overlay. Fresh v4 timing requires either
new behavior-equivalent paired commits carrying the canonical rename or a
separately versioned, fully authenticated external build wrapper. Until then,
v3 replay is the reproducible evidence for this frozen comparison.

The evidence exit statuses are deliberately distinct:

- `0`: valid evidence and every encode promotion/no-regression gate passed;
- `2`: valid, reviewable evidence whose performance policy failed;
- `1`: malformed input, invalid evidence, unsafe isolation, child failure, or
  verified failed-run diagnostics (never a performance result). Portable
  `--no-current-input-check` replay also exits 1 because it is structural
  inspection, not authoritative validation of the current inputs.

The `self-test` command uses only `mock_benchmark.py`. It executes all 24 fixed
cells through five mock A1/B1/B2/A2 rounds, then exercises raw-sample
validation, pass and policy-failure replay, bounded inputs and outputs,
mutation rejection, immutable executable staging, reservation replacement,
stable pair-lease interoperability (including path replacement), shared stable
anchor interoperability with the current butterfly collector, exact-affinity
checking, independently fixed C++ loss vectors, Linux subreaper/pidfd teardown
of a `setsid()` double-fork daemon, and transactional no-replace publication
preflight. The sibling `../test_process_containment.py` gate injects post-spawn
attachment, procfs, pidfd, and teardown failures through all three evidence
runners (18 real runner/fault combinations) and verifies the independent
emergency cleanup path. Neither test
executes either real benchmark binary, and they must not be cited as timing
evidence.

`production-smoke` closes that mock/production boundary without creating
performance evidence. It executes one tiny correctness invocation per requested
backend against both exact real binaries, validates the complete benchmark-v2
schema and checked round trip, and requires the original-data, transmitted-parity,
and repaired-original digests to match between control and candidate. The source
and repaired digests intentionally cover different domains and are not compared
to one another.

See `docs/reproduction/leopard2_low_encode_copy.md` for the exact build,
reservation, collection, and replay commands.

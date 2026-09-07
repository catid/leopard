# Correction: sanitizer instructions are not a Release ISA gate

Bead: `leopard-79h.38.5.4.15`. Date: 2026-09-07.
Resolution: the extra whole-sanitizer-archive scan in `726a202` was applied
outside the production audit's existing scope. No codec or gate fix is needed.

The observed ZMM stores and nonzero scanner exits are real and remain retained.
They do **not** establish a production ISA failure: the project has explicitly
separated sanitizer correctness from unsanitized archive ISA auditing since
commit `e0cfa2fa` on 2026-07-18. `CMakeLists.txt:1650` explains the reason:
instrumentation can insert its own target-specific instructions. Sanitized
configurations register checker self-tests, while unsanitized configurations
register the actual archive audit. A sanitized Release build with
`LEO2_PORTABLE_ISA_RELEASE_AUDIT=ON` is rejected at configure time.

My initial report overlooked that policy and incorrectly made the additional
sanitizer scan a prerequisite for terminal timing. The task is resolved by
verifying the applicable policy, not by forcing a sanitized archive through an
inapplicable gate. Its raw failure has not been rewritten as a pass. The full
Release archive scan and ASan+UBSan correctness checks from the terminal
milestone remain valid and separate. Full Release metadata/v19 qualification,
performance controls and independently linked Leopard1 comparisons remain
outstanding; this correction does not waive any of those requirements.

## Evidence

Frozen source from `726a202` was configured without building or executing codec
targets. Both fields and ordinary correctness-test registration remain enabled.

- Fresh strict Release and ASan+UBSan RelWithDebInfo configurations each pass
  `leopard2_sanitizer_classification` and `leopard2_portable_isa_registration`:
  four CTest passes. The classifier's nine cases include mixed configurations,
  quoted flags, empty and false-like configuration names. Registration checks
  require ordinary field tests to remain present.
- A separately configured sanitized **Release** strict audit fails with the
  explicit unsanitized-Release requirement, proving that renaming the build
  type does not make sanitized code production ISA evidence.
- The unchanged ISA checker passes its complete self-test/negative-control
  mode. No mnemonic, width check or metadata rule was changed.
- CMake and checker files match the working checkout after verification.
  CMake SHA-256 is `465ae5088e4cc6e0b7ddb7719b573ed8de14d81e44116b3865c6b22a45ab5be6`;
  checker SHA-256 is `43ac8cb56b76ad1dc114634e0725babbf51b934c65a7306207eac5a43f26b158`.
- The initial configuration scope peaked at 65,998,848 bytes, final policy
  verification at 40,931,328, and checker controls at 28,794,880, all under
  512 MiB with zero memory events and swap. Compiler probes/fixtures were
  serialized through the canonical lock; no full codec build was needed.

The initial negative fixture used RelWithDebInfo, so it hit the earlier
Release-build-type guard instead of the intended sanitizer guard. Its script
and failing harness log are retained; the corrected sanitized-Release fixture
then reached the intended guard. No project code was changed for either check.

The read-only `sanitizer-isa-policy.DVv0Pp` bundle under
`.research/leopard-79h/` is retained locally and on ripper, including frozen
policy/test sources, generated inventories, configure output, both fixture
attempts and resource logs. The preceding `gfni-terminal.nEDGNa` bundle and
its original failure/replay records remain untouched.

Next: prepare and qualify a wrapper-free terminal-only timing front end.
First-stage forwarding stays OFF, and the original 5% target, 2% control,
immutable-artifact and zero-sibling timing gates remain unchanged. No timing
attempt is authorized by this report alone; preregistration is still required.
Review is Codex self-review and deterministic native checks, with no Claude
or independent-model `CONVERGED` claim.

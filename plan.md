BEGIN CODEX PROMPT

COPY-SAFE FORMAT NOTE

This instruction is intentionally plain text. It contains no Markdown code fences and no inline backtick quoting. Copy everything from BEGIN CODEX PROMPT through END CODEX PROMPT into Codex as one instruction.

MISSION

Implement Leopard2, a production-grade systematic Reed-Solomon erasure codec, directly on top of the existing Leopard-RS codebase:

https://github.com/catid/leopard

Do not stop at a design document or patch sketch. Inspect the code, initialize project tracking, establish independent correctness oracles, implement the high-rate and low-rate encoders and decoders, run the full correctness and performance program, iterate on failures and regressions, make reviewable local commits, and leave a clean repository plus an evidence-based final report.

Assume the current working directory is a checkout of catid/leopard. If it is not, clone the repository and enter it before continuing. Do not overwrite unrelated local changes. Inspect git status before any modification.

The target machine has 128 available CPU cores. Use all 128 cores for compilation, exhaustive and differential testing, fuzzing, variant search, batch scaling, NUMA experiments, and independent experiment sweeps whenever those workloads are genuinely parallel. Do not corrupt benchmark validity merely to display full utilization. Authoritative single-core, cache-sensitive, and memory-bandwidth measurements must be pinned and isolated; immediately return to full-machine parallel work after each isolated measurement phase. Record why any phase intentionally used fewer than 128 cores.

OPERATING DEFINITIONS

Use these terms consistently in code, tests, documentation, and reports.

K = number of original or systematic data shards visible to the application.
R = number of transmitted recovery or parity shards visible to the application.
K + R = transmitted code length.
N = internal power-of-two parent-code length unless an explicitly versioned exact-size profile is used.
P = ceil_pow2(K), used by the low-rate profile.
T = ceil_pow2(R), used by the high-rate profile.
GF8 = GF(2^8), with the legacy Leopard polynomial, Cantor basis, symbol representation, and coordinate order when compatibility is claimed.
GF16 = GF(2^16), with the legacy Leopard representation when compatibility is claimed.
LCH-FFT = additive FFT in the Lin-Chung-Han polynomial basis used by Leopard and the cited papers.
Shortening = fixing parent-code systematic coordinates to known zero values and not transmitting them.
Puncturing = omitting selected parent-code parity coordinates from transmission.
Erasure = a missing symbol whose location is known.
Unknown error = a corrupt symbol whose location is not known; this is outside the production erasure hot path.
Wire profile = the complete versioned definition of field, coordinate order, parent construction, shortening, puncturing, and parity ordering. A CPU-dependent optimization may change kernels, but it may not silently change the wire profile.
Plan setup = work depending on K, R, the profile, and an erasure pattern but not on shard byte values.
Execution = byte-heavy encoding or decoding work applied to one or more stripes after reusable setup.

SUCCESS PRIORITY

Correctness, MDS behavior, deterministic wire format, old-data compatibility where claimed, memory safety, and reproducible evidence outrank an attractive benchmark number. Experimental code must never become the default merely because one microbenchmark is favorable.

1. NON-NEGOTIABLE OBJECTIVES

Build a new leopard2 implementation on top of the existing Leopard codebase with:

1. A legacy-compatible high-rate profile using the existing Leopard coordinate layout and parity format wherever bit-for-bit compatibility can be proven.
2. A new low-rate profile optimized when the number of original shards is the smaller side.
3. Fast systematic encoding and erasure decoding:
   - high-rate work proportional to the redundancy-side transform depth;
   - low-rate work proportional to the message-side transform depth;
   - separation of erasure-pattern setup from the byte-heavy execution step.
4. Arbitrary public original_count = K and recovery_count = R, including R > K, initially through shortening and puncturing around power-of-two LCH parent codes.
5. A safe scalar oracle, mature SIMD backends, runtime dispatch, reusable erasure plans, direct small-loss paths, batched operation, and scalable multicore execution.
6. Preservation of the old public API and old encoded data. Add a new API rather than silently changing old contracts.
7. Experimental evaluation of exact-size transforms and newer implementation ideas. Promote only those that are correct, measurably useful, maintainable, and compatible with the declared profile.
8. No unknown-error-correction work in the normal erasure hot path. Reserve a clean extension point and run a separate experiment only after the erasure codec is complete.

Correctness, MDS behavior, wire-format stability, and reproducible evidence are more important than an attractive benchmark number.

2. REQUIRED SOURCE MATERIAL AND COMPLETE URL BIBLIOGRAPHY

Read the local repository first. At minimum inspect:

README.md
leopard.h
leopard.cpp
LeopardCommon.h
LeopardFF8.h
LeopardFF8.cpp
LeopardFF16.h
LeopardFF16.cpp
CMakeLists.txt and other build files
tests/
Benchmarks.md
docs/ if present
open issues and recent commits relevant to fields, ARM, recovery counts, progressive encoding, and errors versus erasures

Create an ignored research cache such as .research/leopard2/. Download papers and reference implementations there when useful. Do not commit third-party PDFs, generated archives, or copied repositories unless repository policy and licensing clearly justify it. Record exact URLs, retrieval dates, and commit SHAs in docs/leopard2_math_and_sources.txt or docs/leopard2_math_and_sources.md.

The following bibliography is part of this task. Every external source mentioned elsewhere in this prompt is identified here with a URL. Read the primary papers before transcribing algorithms. Use implementations as differential references, not as substitutes for understanding the mathematics.

PROJECT, WORKFLOW, AND COMPARISON IMPLEMENTATIONS

REFERENCE R01
Name: Leopard-RS
Author/owner: catid
Repository: https://github.com/catid/leopard
Role: Production base, legacy API, field representation, coordinate layout, mature scalar and SIMD kernels, current high-rate parity format, compatibility oracle, and baseline.

REFERENCE R02
Name: XDRS
Author/owner: fastecc; repository README identifies Chao Chen as author
Repository: https://github.com/fastecc/xdrs
Role: Executable GF(256), N=256 research reference for switchable low-rate and high-rate LCH-FFT encoding and decoding. Treat as research code. Explicitly select and test the intended algorithms rather than assuming the default benchmark macro does so.

REFERENCE R03
Name: Beads, bd
Repository: https://github.com/gastownhall/beads
Documentation: https://beads.gascity.com/
Role: Dependency-aware issue and experiment tracking. This task requires bd init and the Codex integration.

REFERENCE R04
Name: Intel ISA-L
Repository: https://github.com/intel/isa-l
Role: External erasure-coding implementation and benchmark comparison where its supported code families and parameters are comparable.

REFERENCE R05
Name: Jerasure
Repository: https://github.com/ceph/jerasure
Role: External erasure-coding comparison and source of historical SIMD/table design context.

REFERENCE R06
Name: FastECC
Repository: https://github.com/Bulat-Ziganshin/FastECC
Role: Related fast erasure-coding implementation and comparison target.

REFERENCE R07
Name: ECC-Benchmark
Repository: https://github.com/Bulat-Ziganshin/ECC-Benchmark
Role: Cross-library benchmark methodology and optional comparison harness. Verify fairness and parameter equivalence before quoting results.

REFERENCE R08
Name: SSE2NEON
Repository: https://github.com/DLTcollab/sse2neon
Role: Existing portability aid used by Leopard. Native NEON kernels should still be evaluated rather than assuming translated x86 intrinsics are optimal.

PRIMARY LCH-FFT AND REED-SOLOMON PAPERS

REFERENCE R10
Title: Two Fast Erasure Decoding Algorithms for Reed-Solomon Codes Based on LCH-FFT
Authors: Chao Chen, Sian-Jheng Lin, Nianqi Tang, Yunghsiang S. Han, Suihua Cai, Leilei Yu, Zhongwei Li, Baoming Bai, and Bo Bai
PDF: https://i4ai.org/hanyunghsiang/IT2026.pdf
Role: Primary source for the low-rate and high-rate erasure decoders. Reconstruct the exact low-rate Algorithm 4 and high-rate Algorithm 5 behavior, including the message-only high-rate optimization, locator factors, interpolation identities, and complexity decomposition. The paper uses full-length notation; active-parent generalization is a proof and testing obligation.

REFERENCE R11
Title: Parallel Welch-Berlekamp Algorithm
Authors: Chao Chen, Yunghsiang S. Han, Nianqi Tang, Xiao Ma, and Baoming Bai
PDF: https://i4ai.org/hanyunghsiang/IT2025.pdf
Role: Unknown-error key-equation research, including parallel, early-terminating, and frequency-domain variants. Do not place this machinery in erasure-only encoding or decoding. Use only for a later isolated errors-plus-erasures experiment and hardware-oriented analysis.

REFERENCE R12
Title: Fast Error and Erasure Decoding Algorithm for Reed-Solomon Codes
Authors: Nianqi Tang, Chao Chen, and Yunghsiang S. Han
PDF: https://i4ai.org/hanyunghsiang/CL2024.pdf
Role: Arbitrary-parameter error-and-erasure context, partial FFT use, practical RS(255,223) evidence, and future optional unknown-error work. It does not replace a specialized pure-erasure path.

REFERENCE R13
Title: New Decoding of Reed-Solomon Codes Based on FFT and Modular Approach
Authors: Nianqi Tang and Yunghsiang S. Han
Abstract: https://arxiv.org/abs/2207.11079
PDF: https://arxiv.org/pdf/2207.11079
Role: Detailed arbitrary-epsilon FFT and IFFT algorithms, arbitrary-parameter construction, fast modular approach, and implementation details. Read Appendix A. Its arbitrary-(n,k) code definition is not automatically identical to classical RS or a legacy Leopard profile; treat it as a new profile unless generator-matrix equivalence is proved.

REFERENCE R14
Title: Fast Encoding and Decoding Algorithms for Arbitrary (n,k) Reed-Solomon Codes Over GF(2^m)
Authors: Nianqi Tang and Yun Lin
IEEE page: https://ieeexplore.ieee.org/document/8955804/
Role: Partial-FFT construction for arbitrary public parameters. Use with the same wire-equivalence caution as R13.

REFERENCE R15
Title: Fast Transforms over Finite Fields of Characteristic Two
Author: Nicholas Coxon
Abstract: https://arxiv.org/abs/1807.07785
PDF: https://arxiv.org/pdf/1807.07785
Role: Truncated additive transforms and conversions among LCH, Lagrange, Newton, and monomial bases. This is the preferred starting point for pruning a dyadic parent while preserving its wire profile.

REFERENCE R16
Title: Novel Polynomial Basis and Its Application to Reed-Solomon Erasure Codes
Authors: Sian-Jheng Lin, Wei-Ho Chung, and Yunghsiang S. Han
Abstract: https://arxiv.org/abs/1404.3458
PDF: https://arxiv.org/pdf/1404.3458
Role: Foundational LCH basis and additive-transform work underlying Leopard. Also reconcile this preprint with the journal citation and bundled documents named in the Leopard repository.

REFERENCE R17
Title: FFT Algorithm for Binary Extension Finite Fields and its Application to Reed-Solomon Codes
Authors: Sian-Jheng Lin, Tareq Y. Al-Naffouri, and Yunghsiang S. Han
Abstract: https://arxiv.org/abs/1503.05761
PDF: https://arxiv.org/pdf/1503.05761
Role: Reformulated binary-extension-field FFT and fast RS decoding context. Useful for active subspaces, syndrome transforms, and the legacy mathematical lineage.

REFERENCE R18
Title: Reed-Solomon Coding Algorithms Based on Reed-Muller Transform for Any Number of Parities
Authors: Leilei Yu, Sian-Jheng Lin, Hanxu Hou, and Zhengrui Li
IEEE page: https://ieeexplore.ieee.org/document/10086680/
Role: Alternative arbitrary-parity-count encoder and erasure decoder. Use as a correctness and performance comparison for non-power-of-two redundancy and as a source of potentially useful transform ideas.

REFERENCE R19
Title: A Fast Decoding Algorithm for Generalized Reed-Solomon Codes and Alternant Codes
Authors: Nianqi Tang, Yunghsiang S. Han, Danyang Pei, and Chao Chen
Abstract: https://arxiv.org/abs/2502.02356
PDF: https://arxiv.org/pdf/2502.02356
Role: Generalized RS and arbitrary-evaluation-set context. Relevant if a future profile cannot be expressed as a simple shortened or punctured LCH prefix code.

REFERENCE R20
Title: Efficient Erasure Decoding of Reed-Solomon Codes
Author: Frederic Didier
Abstract: https://arxiv.org/abs/0901.1886
PDF: https://arxiv.org/pdf/0901.1886
Role: Generic erasure-decoding algorithms and complexity comparison outside the LCH-specific construction.

REFERENCE R21
Title: Polynomial Evaluation and Interpolation on Special Sets of Points
Authors: Alin Bostan and Eric Schost
PDF: https://mathexp.eu/bostan/publications/BoSc05.pdf
Role: Fast generic multipoint evaluation and interpolation. Use as an independent exact-size fallback and low-rate comparison where K is much smaller than the padded parent.

REFERENCE R22
Title: Notes on the Truncated Fourier Transform
Author: Joris van der Hoeven
PDF: https://www.texmacs.org/joris/tft/tft.pdf
Role: General truncated-transform design principles, jump removal near powers of two, and practical implementation ideas. Adapt only after proving applicability to additive finite-field transforms.

CURRENT-LITERATURE REFRESH ENDPOINTS

REFERENCE R23
ArXiv search: https://arxiv.org/search/?query=Reed-Solomon+LCH-FFT&searchtype=all
IEEE Xplore search: https://ieeexplore.ieee.org/search/searchresult.jsp?newsearch=true&queryText=Reed-Solomon%20LCH-FFT
DBLP search: https://dblp.org/search?q=Reed-Solomon%20LCH-FFT
Role: Before design freeze and again before final reporting, search for later papers that cite or improve R10-R22, especially work published after the dates of those papers. Record queries, date, candidates reviewed, and why each does or does not change the implementation plan. Prefer papers, official author copies, and publisher pages over summaries.

PROCESSOR, PROFILING, AND TESTING REFERENCES

REFERENCE R30
Name: Intel Intrinsics Guide
URL: https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html
Role: Authoritative intrinsic semantics, throughput and latency notes, and ISA requirements for SSSE3, AVX2, AVX-512, GFNI, and carryless multiplication experiments.

REFERENCE R31
Name: Intel 64 and IA-32 Architectures Software Developer Manuals
URL: https://www.intel.com/content/www/us/en/developer/articles/technical/intel-sdm.html
Role: Authoritative instruction semantics and feature detection details.

REFERENCE R32
Name: Arm Neon Intrinsics Reference
URL: https://arm-software.github.io/acle/neon_intrinsics/advsimd.html
Role: Native Neon implementation reference.

REFERENCE R33
Name: Arm C Language Extensions
URL: https://arm-software.github.io/acle/main/acle.html
Role: Arm feature detection, Neon, SVE, SVE2, and related intrinsic contracts.

REFERENCE R34
Name: Linux perf documentation
URL: https://perfwiki.github.io/
Role: Reproducible performance-counter collection when available.

REFERENCE R35
Name: Agner Fog optimization manuals and instruction tables
URL: https://www.agner.org/optimize/
Role: Supplemental microarchitecture guidance. Treat measured end-to-end results on the target machines as authoritative.

REFERENCE R36
Name: RISC-V Vector Extension specification
URL: https://github.com/riscv/riscv-v-spec
Role: Optional RVV backend experiment if suitable hardware and toolchain are present. No release dependency may be introduced merely for this experiment.

REFERENCE R37
Name: RISC-V vector intrinsic documentation
URL: https://github.com/riscv-non-isa/rvv-intrinsic-doc
Role: Optional RVV intrinsic implementation details.

REFERENCE R38
Name: NVIDIA CUDA documentation
URL: https://docs.nvidia.com/cuda/
Role: Optional GPU batch experiment only when supported hardware and toolchain are already available.

REFERENCE R39
Name: AMD HIP documentation
URL: https://rocm.docs.amd.com/projects/HIP/
Role: Optional GPU batch experiment only when supported hardware and toolchain are already available.

REFERENCE R40
Name: LLVM libFuzzer documentation
URL: https://llvm.org/docs/LibFuzzer.html
Role: Public-API and parser/state-machine fuzzing.

REFERENCE R41
Name: Clang AddressSanitizer documentation
URL: https://clang.llvm.org/docs/AddressSanitizer.html
Role: Memory-safety verification.

REFERENCE R42
Name: Clang UndefinedBehaviorSanitizer documentation
URL: https://clang.llvm.org/docs/UndefinedBehaviorSanitizer.html
Role: Undefined-behavior verification.

REFERENCE R43
Name: Clang ThreadSanitizer documentation
URL: https://clang.llvm.org/docs/ThreadSanitizer.html
Role: Concurrent context, codec, plan, batch, and thread-pool verification.

SOURCE-USAGE RULES

1. Read R10 in full before implementing either specialized decoder.
2. Read the relevant local Leopard files and R16-R17 before changing field representation, basis normalization, transform shifts, or coordinate order.
3. Use R02 only as a differential oracle after identifying its selected compile-time algorithm, field representation, coordinate convention, and license obligations.
4. Use R13-R15 and R18-R22 for arbitrary-size and fallback experiments. Do not silently import a different RS code definition under a legacy profile.
5. Keep R11, the Welch-Berlekamp work, out of the normal erasure path.
6. Map every implemented formula, normalization constant, and coordinate map to a source or a clearly labeled new derivation in docs/leopard2_math_and_sources.md.
7. For every new derivation, include an independent direct-algebra test and, where feasible, an exhaustive small-field test.
8. If a URL is inaccessible, record the failure in Beads and use an official mirror, author page, DOI landing page, or locally bundled copy. Record the alternate URL.
9. Preserve all applicable licenses and notices. Prefer clean reimplementation from papers plus differential testing over copying research code.
10. Do not claim that a method is state of the art merely because it appears in this bibliography. Complete the R23 literature refresh and report scope precisely.

Write and commit docs/leopard2_math_and_sources.md. It must explain every change needed for an active parent subspace of size 2^n embedded in GF(2^m). Do not blindly replace m with n. Derive and test normalization constants, subspace polynomials, shifts, basis indices, coordinate maps, shortening and puncturing semantics, locator evaluations, locator derivatives, and Forney-like output factors.

3. INITIALIZE BEADS BEFORE SUBSTANTIVE WORK

This project must use Beads for all task tracking.

At the start, execute the following commands from the repository root:

    set -euo pipefail
    git status --short

    if ! command -v bd >/dev/null 2>&1; then
        echo "bd is required. Install Beads from https://github.com/gastownhall/beads or https://beads.gascity.com/ and then resume." >&2
        exit 1
    fi

    test -d .beads || bd init
    bd setup codex
    bd prime

If bd is absent, use an official installation method documented at https://github.com/gastownhall/beads and verify the installed binary before continuing. Common official methods include brew install beads and npm install -g @beads/bd. Do not clone the Beads repository into Leopard merely to use the CLI.

Read the resulting or updated AGENTS.md and obey it. Do not create Markdown TODO lists as a substitute for Beads.

Create one P0 epic for Leopard2 and a dependency graph containing at least:

- repository baseline and reproducibility;
- mathematical/coordinate specification;
- scalar finite-field and direct-RS oracle;
- active-subspace LCH transform;
- legacy high-rate encoder;
- legacy-compatible high-rate decoder;
- low-rate encoder/profile;
- low-rate decoder;
- arbitrary counts through shortening/puncturing;
- reusable erasure-plan API;
- direct small-loss recovery;
- SIMD/backend refactor;
- 128-core batch/NUMA work;
- correctness, sanitizers, fuzzing, and compatibility;
- benchmark and experiment harness;
- each speculative experiment listed later;
- documentation and release readiness.

Use bd ready, bd show, bd update <id> --claim, bd dep add, bd remember, and bd close throughout. A bead is closed only with concrete evidence: commit, test command, result artifact, benchmark result, or a documented negative result.

If multiple Codex workers or worktrees are used, serialize Beads mutations through one coordinator or a lock because the default embedded store is a single-writer system. Worker processes may write result files independently but must not race the issue database.

4. GIT AND CHANGE-MANAGEMENT RULES

1. Create a working branch such as codex/leopard2.
2. Keep the existing implementation available as the baseline and compatibility oracle.
3. Commit in reviewable milestones; do not make one giant final commit.
4. Do not push, publish a release, or open a pull request unless the environment explicitly authorizes it.
5. Use separate worktrees or compile-time experimental modules for risky variants. Only the coordinator merges successful work.
6. Keep the default build free of experimental behavior.
7. Preserve BSD notices from Leopard. If source is adapted from Apache-2.0 xdrs, retain the required notices and clearly identify modified files. Prefer reimplementation from the papers and differential testing over wholesale copying.
8. Finish with a clean git status.

5. BASELINE BEFORE MODIFICATION

Before changing algorithms:

1. Build all existing configurations that are supported on the host.
2. Run current tests.
3. Record compiler, flags, CPU topology, NUMA layout, OS, git SHA, and runtime ISA support.
4. Run representative current Leopard encode/decode benchmarks in GF8 and GF16.
5. Save machine-readable baseline results under an ignored or generated-results directory, plus a concise committed summary.
6. Generate legacy golden vectors for boundary and random (K,R) pairs, including parity bytes and recovered outputs.
7. Check whether old FF8 and FF16 outputs are deterministic across scalar, SSSE3, AVX2, and NEON-capable builds.
8. Instrument or count field additions, fixed-constant multiplications, memory passes, and transform butterflies so algorithmic changes can be distinguished from ISA changes.

Do not quote old README throughput as the new baseline.

6. PUBLIC API AND ARCHITECTURE

Add a new public header, preferably leopard2.h, while preserving leopard.h and all existing behavior.

Use opaque, immutable or caller-owned objects along these lines:

    typedef struct leo2_context leo2_context;
    typedef struct leo2_codec leo2_codec;
    typedef struct leo2_decode_plan leo2_decode_plan;

    typedef enum leo2_profile {
        LEO2_PROFILE_AUTO = 0,
        LEO2_PROFILE_LEGACY_HIGH_V1,
        LEO2_PROFILE_LOW_V1,
        LEO2_PROFILE_EXACT_EXPERIMENTAL_V1
    } leo2_profile;

    typedef enum leo2_field {
        LEO2_FIELD_AUTO = 0,
        LEO2_FIELD_GF8,
        LEO2_FIELD_GF16
    } leo2_field;

Provide, with names adapted to repository style:

- context/backend creation and destruction;
- codec creation from (K,R,profile,field,options);
- profile/field/parent-size introspection;
- encode scratch-size query and encode;
- decode-plan creation from missing-coordinate information;
- decode scratch-size query and plan execution;
- one-shot decode wrapper;
- optional parity rebuild;
- batch encode/decode;
- explicit thread-count/executor options;
- result strings and robust error reporting.

Required semantics:

- No allocation in encode/decode hot paths after context/codec/plan creation.
- Immutable codec and plan objects are safe for concurrent execution.
- Caller-provided scratch is size-queryable and alignment-queryable.
- The new API handles arbitrary positive byte lengths, including vector tails. Do not weaken the old API's documented behavior.
- Define input/output aliasing precisely and test it.
- Default decode restores missing originals only. Parity rebuilding is explicit.
- A no-loss decode is a true no-op in the new API.
- Old leo_* APIs continue to work. Route them through Leopard2 only after bit-for-bit and behavior compatibility is proven.
- Old API restrictions, including any R <= K or buffer-size rules, must not silently change. The new API may remove those restrictions.
- Profiles and fields are part of the persistent code identity. Never silently decode a codeword under a different profile or field.

Suggested internal separation:

    src/leopard2/
      context
      codec
      coordinate_map
      plan
      encode_high
      encode_low
      decode_high
      decode_low
      decode_direct
      lch_scalar
      locator
      scheduler
      backends/
        scalar
        ssse3
        avx2
        avx512
        gfni
        neon
        sve2
    tests/leopard2/
    bench/leopard2/
    experiments/leopard2/
    docs/

Adapt this to the repository rather than forcing gratuitous reorganization.

7. COORDINATE PROFILES

7.1 LEGACY-COMPATIBLE HIGH PROFILE

For public counts K and R, define:

    T = ceil_pow2(R)
    N = ceil_pow2(K + T)
    D = N - T
    S = D - K

Use parent code [N,D] with coordinate order:

    [ parity: 0 .. T-1 ]
    [ actual message: T .. T+K-1 ]
    [ shortened known-zero message: T+K .. N-1 ]

Transmit parity coordinates 0 .. R-1 and puncture R .. T-1.

This is expected to match the current Leopard construction, but compatibility is a test result, not an assumption. Prove it for all practical boundary values and a large random sample in GF8 and GF16. If any case differs:

- preserve the old encoder for LEO2_PROFILE_LEGACY_HIGH_V1;
- put the mathematically cleaner version in a distinct profile;
- document the exact difference.

Implement the high-rate encoder by accumulating shifted T-point inverse LCH transforms of message blocks and applying a truncated T-point forward transform to the requested parity outputs. Retain and improve the existing zero-skip, partial-final-block, XOR-accumulation, and fused-layer behavior.

7.2 NEW LOW PROFILE

Define:

    P = ceil_pow2(K)
    N = ceil_pow2(P + R)

Use parent code [N,P] with coordinate order:

    [ actual message: 0 .. K-1 ]
    [ shortened known-zero message: K .. P-1 ]
    [ transmitted parity: P .. P+R-1 ]
    [ punctured parity: P+R .. N-1 ]

Implement systematic low-rate encoding by:

1. Inverse-transforming the padded P message values to LCH coefficients.
2. Evaluating those coefficients on parity blocks with shifted P-point transforms.
3. Skipping completely punctured blocks and truncating the final partial block.
4. Avoiding repeated whole-block copies where an out-of-place or tiled evaluation can be used.
5. Supporting R > K.

The effective transmitted code is [K+R,K]. Verify the MDS property rather than relying only on the parent-code argument.

7.3 PROFILE SELECTION

For backward compatibility:

- the old API uses the legacy high profile;
- the new API's conservative default uses legacy high when appropriate for compatibility and low profile when R > K;
- a later cost model may select the lower measured cost for new code identifiers.

Never make AUTO nondeterministic. If machine-dependent autotuning affects parity format, it must select only kernels within a fixed profile, not a different code.

7.4 FIELD SELECTION

Initially:

- GF8 when the selected parent fits in 256 coordinates;
- GF16 when it fits in 65,536;
- explicit override permitted;
- reject impossible combinations safely.

The selected field is part of the code identifier. Add an experiment for cases where dyadic padding unnecessarily pushes a code from GF8 to GF16 even though K+R <= 256.

8. TRANSFORM AND FINITE-FIELD FOUNDATION

First implement a simple, obviously correct scalar layer:

1. GF8 and GF16 add, multiply, inverse, exponent/log conversion, and basis conversion.
2. LCH forward and inverse transforms parameterized by:
   - active transform length 2^n;
   - field size 2^m, with n <= m;
   - shift/coset offset;
   - known nonzero input count;
   - requested output mask or count.
3. Formal derivative in the normalized LCH basis for the active parent.
4. Direct multipoint evaluation/interpolation for use as an oracle.
5. Generation of skew factors, subspace-polynomial constants, normalization factors, coefL, and coefH for every supported active parent.
6. A test-only small field such as GF(2^4) for exhaustive proofs.

Then optimize without changing results:

- retain radix-2 as the oracle;
- port the existing fused two-layer/radix-4 kernels;
- generalize dependency pruning beyond the current final FFT;
- precompute branch-free transform schedules in codecs/plans;
- skip known-zero leaves and unused output subtrees;
- specialize multipliers 0 and 1 to memset/XOR/copy;
- fuse copy, zero-fill, multiplication, and the first transform layer where profitable;
- fuse a final transform layer with scatter and inverse-locator multiplication;
- fuse block IFFT output into XOR accumulation for encoding/syndrome construction;
- tile over bytes to keep working sets cache-resident.

Do not duplicate nearly identical FF8/FF16 algorithms if a template or generated-kernel layer can share structure without slowing the hot path.

8A. ENCODER-SPECIFIC PRODUCTION REQUIREMENTS

Encoding is a first-class deliverable, not a side effect of decoder work. Establish a direct systematic generator-matrix oracle before optimizing either profile.

High-rate legacy-compatible encoder requirements:

1. Preserve the legacy field representation, coordinate order, parity ordering, and exact output bytes wherever LEO2_PROFILE_LEGACY_HIGH_V1 claims compatibility.
2. Refactor the existing block-IFFT accumulation and truncated final FFT into planned kernels without changing arithmetic order unless bit identity is proved.
3. Skip shortened zero inputs, skip completely unused parity subtrees, handle a partial final message block, and produce only requested transmitted parity coordinates.
4. Fuse copy, zero-fill, fixed multiplication, IFFT butterflies, and XOR accumulation where measurement supports it.
5. Add direct paths for R=1 and very small R.
6. Make the encoder usable through a reusable codec plan and through a progressive accumulator experiment.
7. Test encoding from unaligned buffers, arbitrary positive byte counts, batches, and multiple threads.

Low-rate encoder requirements:

1. Pad the K systematic values to P=ceil_pow2(K) with known shortened zeros.
2. Perform one P-point inverse transform to obtain coefficients in the defined active-parent basis.
3. Evaluate only transmitted parity blocks and only requested outputs, with shifted P-point transforms or a proven exact-size alternative.
4. Avoid copying the coefficient block once per parity block; use tiled out-of-place evaluation, reusable scratch, or generated schedules.
5. Add direct paths for tiny K, tiny R, and small total byte counts.
6. Support R greater than K and test maximum legal redundancy for GF8 and GF16 profiles.
7. Re-encoding after message recovery must reproduce parity exactly.

Common encoder requirements:

1. Precompute code-dependent coefficients, skew factors, coordinate maps, and pruned schedules in the immutable codec object.
2. Provide a scratch-size query, an alignment query, and no hot-path allocation.
3. Define exact input/output aliasing rules. Reject unsupported overlap rather than invoking undefined behavior.
4. Make output deterministic across scalar and SIMD backends.
5. Count memory passes, bytes read and written, XORs, and fixed multiplications in test instrumentation.
6. Expose batch encoding so short shards and many stripes can fill SIMD units and all CPU cores.
7. Permit generating any requested parity subset without first materializing all R parity shards.
8. Add an experimental incremental-parity-update API, but do not expose it as stable until it is proven to reproduce full encoding exactly.

9. DECODER IMPLEMENTATIONS

9.1 RETAIN A GENERIC FALLBACK

Keep the current full O(N log N) LCH erasure decoder as:

- a compatibility oracle;
- a fallback for unsupported patterns/configurations;
- a benchmark baseline.

Correctness may never depend solely on the new specialized algorithms.

9.2 LOW-RATE DECODER

Implement IT2026 Algorithm 4 for the active low-profile parent:

- separate locator/normalization setup from byte-heavy execution;
- perform blockwise P-point inverse transforms;
- compute only the required derivative/weighted combination;
- perform the final P-point transform;
- recover erased message coordinates only;
- rebuild parity by re-encoding only when requested.

Target main-step work proportional to N log P, excluding reusable plan setup.

9.3 HIGH-RATE DECODER

Implement IT2026 Algorithm 5, and the message-only optimization, for the active high-profile parent:

- form the T-coefficient syndrome-like polynomial from blockwise inverse transforms;
- apply locator factors;
- inverse-transform to the evaluator polynomial;
- evaluate only blocks/subtrees containing requested missing originals;
- omit parity-only derivative/recovery work unless explicitly requested.

Target main-step work proportional to N log T, excluding reusable plan setup.

9.4 ACTIVE-SUBSPACE PROOF OBLIGATION

The paper presents full-field notation. Existing Leopard uses smaller power-of-two parent subspaces inside GF16. Before enabling either specialized decoder:

- derive the active-N form of every use of s_m, basis index, derivative, and Forney-like normalization;
- validate it symbolically or exhaustively in GF(2^4);
- validate it against direct interpolation in GF8/GF16;
- validate GF8 full-length results against xdrs;
- validate the legacy high profile against existing Leopard.

A failed derivation must leave the generic decoder selected rather than being papered over with test exceptions.

10. ERASURE PLAN AND LOCATOR IMPROVEMENTS

Create a reusable leo2_decode_plan containing:

- fixed shortened/punctured positions;
- dynamic missing positions;
- deterministic selection of enough received coordinates;
- virtual erasures for surplus received symbols;
- locator evaluations at all input positions actually used;
- inverse locator derivatives or equivalent output factors;
- active block list and pruned transform schedule;
- direct-repair coefficients when the dispatcher selects that path;
- output scatter map;
- backend-formatted fixed multipliers when useful.

Precompute permanent contributions from punctured positions in the codec.

Implement and compare locator construction strategies:

1. Existing FWHT-based evaluation, restricted to the active parent rather than the entire field order.
2. Sparse direct construction for few erasures.
3. Product-tree or epsilon-transform construction where appropriate.
4. Incremental update for a sequence of patterns with small Hamming distance.

For incremental plans, exploit the additive indexing relation omega_i + omega_j = omega_(i xor j) in characteristic two. Investigate updating log-locator evaluations in O(N) per added/removed erasure, while handling the zero/self term and derivative-at-erasure values correctly. Do not promote this until exhaustive tests cover add, remove, swap, puncture, and virtual-erasure cases.

When more than K symbols are available, experiment with selecting the subset that minimizes pruned-transform work. Compare:

- lowest-index deterministic choice;
- prefer surviving systematic shards;
- block-aligned greedy selection;
- exact/dynamic-programming selection for small N.

The chosen subset must not affect decoded data, only cost.

11. DIRECT AND HYBRID PATHS

Transform algorithms are not always fastest.

Implement:

1. R=1 XOR parity and one-loss decode.
2. Tiny K, tiny R, and tiny loss-count direct matrix paths.
3. A small-loss repair solver that uses surviving originals plus a selected set of parity equations, solves only an L x L system for L missing originals, and precomputes byte-execution coefficients in the plan.
4. Direct systematic evaluation for parameter cases where dyadic padding is excessive.
5. A dispatcher driven by measured cost, not asymptotics alone.

The dispatcher should consider:

    K, R, parent N, P or T, number of missing originals,
    number of available parity shards, shard bytes, batch size,
    plan reuse count, field, SIMD backend, and thread count.

Provide deterministic built-in thresholds and an optional offline-generated calibration table. Never benchmark in the first user call unless the API explicitly requests autotuning.

12. SIMD AND BACKEND WORK

Refactor fixed-constant multiplication, multiply-add, XOR, copy, and butterfly operations behind compile-time-specialized backends with runtime dispatch.

Required production backends where supported by the repository/host:

- scalar;
- existing SSSE3;
- existing AVX2;
- native ARM NEON rather than depending solely on translation wrappers.

Build ISA-specific translation units or target functions. Do not compile the entire library with an ISA flag that can produce illegal instructions before runtime dispatch. Add backend self-tests during initialization.

Evaluate:

12.1 AVX-512 BW/VBMI

- 512-bit nibble-table multiplication;
- wider XOR/copy;
- radix-4 and radix-8/fused-three-layer butterflies;
- register pressure and spill behavior;
- AVX-512 frequency effects.

Promote only when end-to-end throughput improves on relevant CPUs, not merely when the inner instruction count falls.

12.2 GFNI AFFINE FIXED-MULTIPLIER BACKEND

Multiplication by a fixed element of GF(256) is an 8-by-8 GF(2)-linear map. Precompute the matrix for every multiplier in Leopard's actual Cantor/polynomial representation and test whether VGF2P8AFFINEQB can apply that fixed map directly to each byte.

Required validation:

- determine the instruction's matrix bit order empirically;
- exhaust all 256 multipliers and all 256 byte inputs;
- test XMM/YMM/ZMM forms that the host supports;
- compare against scalar and current nibble-table paths;
- include coefficient 0/1 specialization;
- benchmark complete transforms and codecs, not only mul_mem.

Also test VGF2P8MULB with explicit basis conversion as an alternative, but do not convert the stored wire representation merely for the kernel.

12.3 GF16 ALTERNATIVES

Experiment with:

- improved existing 16-bit table layout;
- PCLMULQDQ/VPCLMULQDQ constant multiplication and reduction;
- a quadratic composite-field GF((2^8)^2) representation using GFNI;
- linear isomorphism to the legacy field for wire compatibility.

Any representation change must either map exactly to legacy symbols at the API boundary or use a new field/profile identifier.

12.4 ARM

Implement or prototype:

- native NEON table lookups and fused butterflies;
- PMULL constant multiplication;
- SVE/SVE2 scalable-vector kernels when toolchain and hardware are available.

Cross-compile-only success is not a performance claim.

13. MULTICORE, BATCH, AND NUMA DESIGN

Single-stripe transform stages can require barriers and may become memory-bandwidth bound. Prefer parallelism in this order:

1. independent stripes/codewords in a batch;
2. independent codecs/plans;
3. byte tiles of large shards;
4. independent transform blocks;
5. within-stage butterflies only when granularity justifies synchronization.

Add a batch API so small shards can fill vectors and cores.

Use:

- a persistent optional thread pool or caller executor;
- no per-call thread creation in the hot path;
- static scheduling for regular transform work;
- per-worker scratch;
- first-touch allocation;
- NUMA-node partitioning;
- core affinity in benchmarks;
- separate per-node accumulators where cross-node writes would hurt.

Measure thread counts 1,2,4,8,16,32,64,128. Report both aggregate input throughput and generated/repaired-output throughput. Identify the memory-bandwidth saturation point rather than claiming linear scaling.

14. USE ALL 128 CORES PRODUCTIVELY

Detect the allowed CPU set instead of assuming CPUs are numbered 0..127.

Create a reproducible experiment runner, for example:

    tools/leopard2_lab.py

It must:

- detect process affinity, logical/physical cores, sockets, and NUMA nodes;
- generate a deterministic JSON manifest of jobs;
- assign stable random seeds;
- pin jobs to allowed CPUs or CPU groups;
- cap per-job memory;
- write one result JSON per job;
- resume incomplete runs;
- merge results deterministically;
- retain stdout/stderr for failures;
- use timeouts;
- print live utilization/progress;
- default to 128 workers when 128 CPUs are allowed.

Use all 128 cores for:

- cmake --build ... -j128;
- ctest -j128;
- exhaustive small-field cases sharded by (N,K,R,pattern);
- GF8 differential cases;
- sanitizer build matrices;
- 128 independent fuzz seeds;
- compiler/ISA/tile/fusion variant search;
- all-core batch and NUMA scaling;
- exact-size-transform sweeps;
- direct-vs-transform crossover sweeps.

Typical bootstrap:

    CPU_COUNT="$(nproc)"
    JOBS="$CPU_COUNT"
    if [ "$JOBS" -gt 128 ]; then JOBS=128; fi
    export CMAKE_BUILD_PARALLEL_LEVEL="$JOBS"

    cmake -S . -B build/release \
      -DCMAKE_BUILD_TYPE=Release \
      -DLEO2_BUILD_TESTS=ON \
      -DLEO2_BUILD_BENCHMARKS=ON
    cmake --build build/release -j "$JOBS"
    ctest --test-dir build/release -j "$JOBS" --output-on-failure

Do not hardcode this if the repository uses another supported build flow.

For authoritative microbenchmarks:

- pin to one physical core or one NUMA node;
- warm caches in a controlled way;
- avoid simultaneous memory-intensive jobs that invalidate results;
- record governor/turbo/frequency state when readable;
- use multiple repetitions and confidence intervals.

Scientific validity overrides the desire to display 128 busy cores for that short phase. Resume full saturation immediately afterward.

15. CORRECTNESS AND VERIFICATION PROGRAM

15.1 SCALAR/DIRECT ORACLE

Implement a mathematically independent direct RS encoder/interpolator for tests. Do not use the optimized transform code to verify itself.

15.2 EXHAUSTIVE SMALL-FIELD TESTS

Using GF(2^4) or another test-only field:

- enumerate every valid small parent/profile;
- encode every basis message and many/all messages as feasible;
- verify systematic form;
- verify every subset of K transmitted coordinates reconstructs;
- verify every erasure pattern up to R;
- compare generic, low, high, direct, pruned, and incremental-plan paths.

15.3 GF8/GF16 DIFFERENTIAL TESTS

Cover counts around every boundary:

    1,2,3,4,7,8,9,15,16,17,31,32,33,63,64,65,
    127,128,129,191,192,193,223,224,225,239,240,
    241,247,248,249,255,256

Use valid combinations and corresponding larger GF16 boundaries.

Cover:

- R=1;
- K=1;
- R>K;
- rate just below/equal/above 1/2;
- K or R just above a power of two;
- K+R <= 256 but dyadic parent >256;
- zero shortened tail;
- large shortened tail;
- partial final transform block;
- all originals present;
- exactly one missing original;
- maximum correctable erasures;
- extra available parity symbols;
- missing parity only;
- mixed original/parity loss;
- repeated plan execution;
- incremental pattern changes.

15.4 BUFFER TESTS

Test byte sizes:

    1,2,3,7,8,15,16,17,31,32,33,63,64,65,
    127,128,129,255,256,257,1023,1024,1025,
    64 KiB, 1 MiB, 16 MiB

Test unaligned buffers, page boundaries, aliasing allowed by the API, null missing pointers, guard pages where practical, and large-count overflow checks.

15.5 LEGACY COMPATIBILITY

For the legacy profile:

- compare parity byte-for-byte with old Leopard;
- compare recovered originals;
- compare behavior/error codes through the old wrappers;
- test GF8/GF16 field choice and endian behavior;
- retain committed compact golden vectors.

Do not route old APIs to the new engine until this gate passes.

15.6 REFERENCE COMPARISON

For full GF8 N=256 power-of-two cases:

- build xdrs locally;
- explicitly select its low Algorithm 2 and high Algorithm 3;
- compare parity and decode results after accounting for coordinate conventions;
- document any expected representation difference.

15.7 TOOLING

Run, in separate builds as appropriate:

- ASan;
- UBSan;
- TSan for concurrent plan/context use;
- MSan if the toolchain supports it without excessive setup;
- Clang and GCC warnings at high levels;
- libFuzzer or AFL-style public-API fuzz targets;
- deterministic stress tests;
- valgrind only where it adds information.

Fuzz workers must use different logged seeds and all 128 cores.

16. BENCHMARK MATRIX AND REPORTING

Benchmark old Leopard and every candidate path across:

RATES/COUNTS

    low:
      (K,R) = (8,248), (16,240), (32,224), (64,192),
              (100,156), (127,129), and larger GF16 analogues

    balanced:
      (128,128), (256,256), larger analogues

    high:
      (192,64), (224,32), (240,16), (248,8),
      (1000,200), (4096,512), and arbitrary non-power cases

    padding stress:
      (129,1), (129,100), (200,50), (225,30),
      and cases where parent inflation changes the field

Adjust invalid pairs while preserving the intent.

SHARD SIZES

    64 B through 16 MiB, logarithmically spaced

LOSS COUNTS

    0,1,2,4,8,R/4,R/2,maximum valid

REUSE/BATCH

    one stripe per plan, 8, 64, 1024, and larger batches

THREADS

    1,2,4,8,16,32,64,128

Record:

- setup time;
- encode/decode execution time;
- input GB/s;
- parity/repaired-output GB/s;
- cycles per byte;
- instructions per byte;
- cache/TLB misses when available;
- memory bandwidth when available;
- scratch bytes;
- code/table footprint;
- operation counts;
- confidence interval or median/MAD;
- compiler and ISA.

Use perf stat or equivalent when available without requiring privileged installation.

A candidate is promoted by evidence, not by one favorable cell. Default promotion rule:

- complete correctness gates;
- statistically credible improvement of at least 5% in its target regime;
- no unexplained regression greater than 2% in key neighboring regimes;
- acceptable code size and maintenance burden;
- dispatcher can avoid losing cases.

Stricter thresholds are appropriate for complicated ISA-specific code.

17. REQUIRED PRODUCTION IMPROVEMENTS

The release path should include, subject to correctness:

1. Legacy-compatible high-rate encode.
2. IT2026 high-rate original-recovery decoder.
3. New low-rate encode/profile.
4. IT2026 low-rate original-recovery decoder.
5. Arbitrary K,R through shortening/puncturing.
6. Parent-subspace rather than full-field locator/transform setup.
7. Reusable immutable erasure plans.
8. Fixed contribution for permanent punctures.
9. Direct paths for no loss, one parity, and few losses.
10. Pruned input/output transform schedules.
11. Tiled low-memory scratch, targeting O(min(P,T)) shard slots rather than O(N) full shard slots where practical.
12. Existing scalar/SSSE3/AVX2 behavior retained and cleaned up.
13. Native NEON path where testable.
14. Runtime field/backend dispatch.
15. Arbitrary byte tails in the new API.
16. Batch API and 128-core scaling.
17. Machine-readable benchmark/test artifacts.
18. Full documentation of wire profiles and compatibility.

18. EXPERIMENTAL RESEARCH TRACKS

Each experiment gets its own Bead with:

- hypothesis;
- implementation plan;
- correctness oracle;
- benchmark cells;
- promotion threshold;
- kill criterion;
- final result, including negative results.

Do not let speculative work block the production codec.

EXPERIMENT A: GFNI AFFINE MULTIPLICATION

Test the direct 8x8 affine-map implementation of every fixed GF8 multiplier. Compare against AVX2 nibble tables, AVX-512 nibble tables, and GF2P8MULB with basis conversion.

Kill if instruction throughput, matrix setup, downclocking, or surrounding loads make end-to-end codecs no faster.

EXPERIMENT B: THREE/FOUR-LAYER BUTTERFLY FUSION AND GENERATED KERNELS

Search radix-2, fused-2, fused-3, fused-4, and mixed schedules by transform size and ISA. Inspect generated assembly and spills. Consider offline code generation or constexpr tables, not runtime executable memory in the default library.

Promote per-size schedules only where whole-codec gains survive.

EXPERIMENT C: NON-POWER-OF-TWO TRANSFORMS AND EXACT-LENGTH CODES

This is a staged research program, not one vague experiment. Create a Beads epic for non-power-of-two work and child beads C0 through C10 below. Run the stages in order unless a stage is clearly independent. Preserve negative results.

Distinguish two separate goals at all times:

1. Exact computation of the existing punctured or shortened dyadic parent code. This must preserve parity bytes and wire compatibility. It removes unnecessary computation but does not change the code.
2. A new exact-(K,R) code profile. This may change coordinate sets, generator matrices, or parity bytes. It must have a new profile identifier and independent MDS proof and compatibility documentation.

Never describe an exact-length implementation as legacy compatible merely because it decodes correctly. Legacy compatibility requires byte-for-byte parity equivalence for the same K, R, field, and input data.

C0. SYMBOLIC COST AND DEPENDENCY SIMULATOR

Before writing irregular SIMD kernels, implement a symbolic simulator for:

- current padded Leopard transforms;
- prefix-pruned dyadic transforms;
- general input/output dependency pruning;
- recursive truncated transforms;
- binary decomposition into aligned dyadic blocks;
- full-block plus exact or dense tail hybrids;
- direct matrix interpolation/evaluation;
- generic subproduct-tree interpolation/evaluation.

For every valid 1 <= K,R <= 256, and representative GF16 ranges, count or estimate:

- radix-2 butterfly equivalents;
- XOR bytes moved;
- fixed-constant multiplications;
- loads and stores;
- temporary shard slots;
- irregular boundary operations;
- transform passes over shard memory;
- parent inflation;
- whether padding forces GF16 despite K+R <= 256.

Produce machine-readable matrices and heat maps for estimated gain:

G(K,R) = padded_cost(K,R) / candidate_cost(K,R).

Concentrate on counts just above and just below powers of two:

3,5,7,9,15,17,31,33,63,65,127,129,255.

Use the simulator to prioritize implementation. Do not assume that an algorithm with fewer field operations has lower memory traffic or lower measured runtime.

C1. WIRE-COMPATIBLE RECURSIVE DEPENDENCY PRUNING

Keep the exact existing dyadic parent code and coordinate map. Build a scalar transform evaluator that recursively tracks:

- whether a subtree has any nonzero dynamic input;
- whether a subtree can affect any requested output.

Skip a node when its result cannot affect a requested output. Specialize fixed multipliers 0 and 1. Eliminate copies, zero fills, and butterflies made unnecessary by permanent punctures or shortened zeros.

Implement and compare four execution forms:

1. recursive branching;
2. a precompiled flat operation list;
3. a list of complete radix-4 or fused subtransforms plus boundary operations;
4. generated kernels for common K,R pairs.

The correctness oracle is the full padded transform. Every output must match bit-for-bit. Measure branch overhead, schedule size, instruction-cache pressure, and memory-traffic reduction over shard sizes from tiny to large.

C2. TRUNCATED DYADIC LCH TRANSFORM PRESERVING THE PARENT CODE

Implement true truncated forward and inverse LCH transforms with:

- an enclosing power-of-two parent;
- an exact active input prefix or sparse active input mask;
- an exact requested output prefix or output mask;
- no materialization of the unused suffix;
- identical output to the padded parent transform.

Represent arbitrary active lengths by complete dyadic subtrees plus a recursive boundary. Run full optimized kernels for complete blocks and irregular logic only on the ragged boundary. Use Coxon's truncated additive-transform and basis-conversion ideas as the main reference.

Required APIs may be internal, conceptually:

lch_truncated_forward(active_inputs, requested_outputs, parent_size, shift)

lch_truncated_inverse(active_values, requested_coefficients, parent_size, shift)

Do not expose these names publicly unless they fit the repository design.

C3. TRUNCATED INVERSE AND BASIS-CONVERSION STUDY

Treat the inverse transform as a separate research problem. Compare:

- padded LCH inverse transform;
- Coxon-style Lagrange-to-LCH conversion;
- Newton-to-LCH conversion;
- Tang-Han epsilon-point inverse transform;
- direct Gaussian elimination or precomputed inverse matrices for small sizes;
- generic subproduct-tree interpolation.

Test arbitrary prefix lengths and shifted cosets. Record setup time separately from shard-byte execution time. Determine crossovers by K, shard bytes, batch size, and plan reuse. Expect dense methods to win for some small K values.

C4. FULL DYADIC BLOCK PLUS EXACT OR DENSE TAIL

For q = 2^a + d with 0 < d < 2^a, test a hybrid:

- execute a normal optimized 2^a-point Leopard transform for the large block;
- handle the d-element tail by one of: direct matrix, smaller padded transform, recursive truncation, or Newton-basis conversion;
- combine with precomputed cross-block factors.

Sweep every d for practical a. Generate an offline decision table selecting the best tail method. This hybrid is a leading production candidate because it preserves regular SIMD execution over most of the problem.

C5. BINARY DYADIC-BLOCK DECOMPOSITION

Decompose an arbitrary q as:

q = 2^a1 + 2^a2 + ... + 2^as, with a1 > a2 > ... > as.

Use aligned additive cosets for each block. Process each block with existing optimized LCH kernels. Investigate cross-block conversion or joining in LCH, Newton, and Lagrange bases.

First build a scalar algebra and operation-count model. Measure:

- cross-block matrix sparsity;
- fixed multiplications;
- scratch;
- whether the overall map becomes dense;
- SIMD regularity;
- MDS behavior of any new coordinate profile.

Abandon this route if cross-block joining becomes dense enough to erase the savings.

C6. EXACT GF256 FIELD-BOUNDARY RESCUE

Prioritize cases where K+R <= 256 but dyadic parent inflation forces GF16. Compare:

1. padded legacy GF16;
2. direct systematic GF256 matrix encode/decode;
3. wire-compatible truncated GF256 where possible;
4. Tang-Han epsilon GF256;
5. binary-block GF256 profile;
6. generic GF256 interpolation/evaluation.

Test counts immediately above powers of two on either side. A new exact GF256 profile is allowed and may be valuable even if it performs somewhat more field operations, because GF256 kernels and tables may be substantially cheaper than GF65536.

C7. EXACT LOW-RATE PROFILE

Develop a separately versioned exact low-rate profile only after C0-C6 establish that it is worthwhile. For arbitrary K:

- interpolate from exactly K systematic evaluations;
- evaluate exactly R parity coordinates;
- avoid padded systematic coordinates in the hot algebra;
- compare coordinate layouts: prefix points, a union of aligned dyadic cosets, and an offline-searched transform-friendly set.

Implement the encoder before adapting the decoder. Prove systematic behavior and MDS behavior exhaustively in a small field and by large random GF8/GF16 differential tests.

C8. EXACT HIGH-RATE PARITY SOLVE

For arbitrary R, formulate high-rate systematic encoding as solving exactly R parity constraints. Compare:

- partial or transposed LCH transforms;
- Tang-Han epsilon transforms;
- Newton interpolation;
- precomputed R-by-R inverse matrices;
- dyadic-block Schur-complement methods;
- the existing T = ceil_pow2(R) regular encoder.

Do not assume exact wins. The padded T-point encoder may remain faster because T < 2R and its kernels are highly regular.

C9. EXACT SPECIALIZED ERASURE DECODERS

Only after exact encoders and coordinate profiles are frozen, adapt the IT2026 specialized decoders.

For low rate, investigate replacing P-point operations with exact-K interpolation and evaluation while preserving the weighted block reduction.

For high rate, investigate replacing the T-coefficient syndrome and evaluator objects with exact-R representations.

For a general coordinate set S, derive the required vanishing polynomial:

Z_S(x) = product over alpha in S of (x - alpha),

and its derivative and quotient evaluations. Determine how these replace the simple subspace-polynomial factors available for additive subspaces.

Do not promote an identity based on numerical success alone. Document the derivation and verify it independently. The extra normalization multiplications may make exact decoding slower even when the transform is smaller.

C10. END-TO-END CROSSOVER MAP AND DISPATCH

Build an offline crossover map covering:

- K and R;
- loss count and pattern shape;
- shard bytes;
- batch size;
- plan reuse count;
- field;
- ISA backend;
- thread count;
- legacy versus new profile.

Compare at least:

- direct matrix;
- padded generic decoder;
- padded IT2026 low/high decoder;
- dependency-pruned parent transform;
- truncated parent transform;
- full-block plus tail hybrid;
- exact profile where implemented.

The expected result is a region-based dispatcher, not one universal winner. Generate a compact deterministic decision table. Do not benchmark online in normal user calls.

NON-POWER-OF-TWO PROMOTION RULES

A candidate may be promoted only when:

1. It passes exhaustive small-field MDS and differential tests.
2. A wire-compatible candidate matches padded Leopard parity byte-for-byte for every claimed profile.
3. Setup and execution are reported separately and at multiple reuse counts.
4. It reduces measured memory traffic or runtime, not only symbolic field-operation count.
5. It provides at least a 10 percent credible improvement in a meaningful parameter region, unless the code is exceptionally simple and broadly useful.
6. A dispatcher avoids regressions outside that region.
7. Schedule, table, generated-code, and instruction-cache costs remain acceptable.
8. Exact-code variants have a new profile identifier, serialized code identity, and documented coordinate map.

Treat Tang-Han's arbitrary-parameter code definition as a new profile unless generator-matrix equivalence to an existing Leopard profile is proved.

EXPERIMENT D: GF8 FIELD-BOUNDARY RESCUE

For cases with K+R <= 256 but ceil_pow2 parent construction exceeds 256, compare:

- dyadic GF16;
- exact/direct systematic GF8;
- epsilon/truncated GF8;
- low versus high profile.

A new exact GF8 profile is acceptable if versioned and demonstrably useful.

EXPERIMENT E: INCREMENTAL LOCATOR PLANS

Update a cached locator/derivative plan when one or a few erasure positions change. Test stream-like patterns, RAID rebuilds, and packet-loss windows. Compare O(N) delta update with FWHT rebuild and sparse direct rebuild.

EXPERIMENT F: RECEIVED-SUBSET OPTIMIZATION

Choose which surplus received symbols to use to minimize active transform blocks and memory traffic. Compare deterministic baseline, greedy block packing, and exact search for small codes.

EXPERIMENT G: COMPOSITE GF16 AND CARRYLESS MULTIPLICATION

Prototype quadratic GF((2^8)^2), GFNI-assisted arithmetic, and VPCLMUL reduction. Include basis-isomorphism cost and wire-compatibility conversion. Reject a fast microkernel that loses end-to-end.

EXPERIMENT H: GF12

The old code mentions 12-bit fields. Prototype GF(4096) in 16-bit containers:

- table/cache footprint;
- maximum coordinate count;
- vector multiplication;
- wasted high bits;
- wire packing cost;
- parent-size sweet spots.

This is a new field/profile only. Do not complicate the default release unless it clearly wins a useful range.

EXPERIMENT I: TRANSFORM-AWARE STREAMING/INCREMENTAL ENCODE

Expose an optional encoder accumulator:

- high profile naturally accumulates block IFFTs;
- investigate a linear low-profile accumulator;
- reduce peak memory;
- support progressive input;
- finalize requested parity.

Output must match ordinary encode exactly.

EXPERIMENT J: PLAN CACHE AND FIXED-PATTERN BATCHING

Implement an optional bounded LRU keyed by codec ID plus erasure bitmap, or expose enough hooks for callers to manage plans. Measure hit/miss overhead and thread contention.

EXPERIMENT K: NUMA TASK GRAPH

Compare:

- one global pool;
- per-NUMA pools;
- stripe partition by node;
- byte partition by node;
- duplicated read-only tables;
- first-touch scratch;
- huge pages where available without requiring privileges.

EXPERIMENT L: OFFLINE AUTOTUNING

Build a reproducible tuner that chooses:

- direct versus transform;
- low versus high for new profiles;
- GF8 versus GF16 where both are legal;
- fusion depth;
- tile size;
- thread count;
- backend.

Generate a compact deterministic decision table. Do not add opaque online benchmarking to normal startup.

EXPERIMENT M: CHECKSUM-ASSISTED CORRUPTION-TO-ERASURE

Many applications can identify a bad shard with a checksum and then use the much faster erasure path. Prototype an optional validation callback or batch prevalidation layer. Keep checksumming outside RS algebra and do not silently trust unauthenticated data.

EXPERIMENT N: UNKNOWN-ERROR CORRECTION EXTENSION

Only after all production erasure milestones are complete:

- reserve an API/module boundary for errors plus erasures;
- use the 2022 fast modular approach as the primary software candidate;
- test IT2025 EPWB/FEPWB as a research implementation or hardware-oriented design;
- ensure none of this changes the erasure fast path or ABI;
- stop if scope threatens release readiness.

EXPERIMENT O: LINEAR-CIRCUIT OR TABLE SYNTHESIS FOR TINY CODES

For fixed tiny (K,R) profiles, search for generated XOR/affine networks and fused multi-source kernels. Compare against direct mul-add and transform paths. Cap generated code size.

EXPERIMENT P: INCREMENTAL PARITY UPDATES AFTER SOURCE CHANGES

Exploit code linearity when one or a few source shards change. Compare:

- direct delta multiply-add into every transmitted parity shard;
- a transform-domain delta update;
- batching several source deltas before finalization;
- full re-encoding.

Require bit-for-bit equality with a fresh encode. Measure update count, changed byte range, K, R, shard size, cache behavior, and break-even points. Consider a range-update API only if semantics remain simple and race-free.

EXPERIMENT Q: SIMD ACROSS INDEPENDENT STRIPES

For tiny shards, vectorizing within one shard may underfill vectors and cores. Prototype array-of-structures, structure-of-arrays, and AoSoA batch layouts that place the same byte offset from multiple independent stripes in SIMD lanes. Include transpose and packing costs, latency, batch thresholds, and API ergonomics. Preserve the ordinary layout at the external API boundary unless a separately versioned packed-batch API is justified.

EXPERIMENT R: CACHE-OBLIVIOUS, SIX-STEP, STOCKHAM, AND TRANSPOSED TRANSFORM LAYOUTS

Compare the existing iterative schedule against recursive cache-oblivious traversal, six-step-style decomposition, Stockham/autosort layouts, and explicitly transposed block evaluation. Measure copies, cache and TLB misses, barrier count, prefetch behavior, and NUMA traffic. The additive LCH transform is not a multiplicative complex FFT; derive each schedule rather than mechanically copying a Cooley-Tukey implementation.

EXPERIMENT S: MEMORY-TRAFFIC AND STORAGE-POLICY OPTIMIZATION

Evaluate software prefetch, non-temporal stores, write-combining, zero-page or shared-zero handling, huge pages, page coloring where practical, scatter versus compact output, in-place versus out-of-place stages, and first/last-layer fusion with application buffers. Promote only end-to-end wins. Never use non-temporal stores when they slow data that is immediately consumed.

EXPERIMENT T: BIT-SLICED AND TOWER-BASIS ARITHMETIC

Prototype table-free bit-sliced GF8, vectorized Boolean linear maps, tower-field representations, and alternate internal bases with explicit conversion at the profile boundary. Compare against nibble tables and GFNI. Account for conversion, transpose, register pressure, code size, and tiny versus large batch behavior. A legacy profile must emit exactly the legacy bytes.

EXPERIMENT U: GENERATED PER-SIZE AND PER-PROFILE KERNELS

Use offline code generation, templates, or constexpr schedules for common transform sizes and profiles. Specialize constant multipliers, zero branches, output masks, and fusion depth. Do not use runtime executable-memory JIT in the default library. Set code-size budgets and verify every generated kernel against the scalar oracle.

EXPERIMENT V: OPTIONAL RVV AND GPU BATCH BACKENDS

If suitable hardware and toolchains are already available, prototype RISC-V RVV and CUDA or HIP batch backends using R36-R39. Keep them in isolated build options, require the scalar wire-format oracle, report transfer overhead, and do not make release completion depend on unavailable hardware. GPU tests must include host-device transfer and realistic batch sizes.

EXPERIMENT W: SERIALIZED CODE IDENTIFIERS AND LONG-TERM WIRE COMPATIBILITY

Design and test a compact, endian-stable serialization for profile version, field, K, R, parent length, padded side, coordinate-map version, and optional algorithm-independent metadata. Include malformed-input tests, forward-compatibility rules, golden vectors, and migration documentation. The serialized identifier must describe mathematics and wire layout, not the current CPU backend or autotuning decision.

19. PERFORMANCE AND CODE-QUALITY TRAPS TO AVOID

- Do not benchmark only 1024-byte GF8 shards and generalize.
- Do not count plan construction as free; report it separately and amortized at several reuse counts.
- Do not use all 128 cores during a single-core latency benchmark.
- Do not select a different wire profile based on the current CPU.
- Do not assume power-of-two padding is harmless near a field boundary.
- Do not materialize an N-shard workspace when a P- or T-wide tiled accumulator suffices.
- Do not run the full-field 65,536-point FWHT for a small active GF16 parent without proving it is best.
- Do not copy all surviving originals on a no-loss decode.
- Do not add branches inside the byte inner loop when a precomputed schedule can remove them.
- Do not promote AVX-512 based on instruction count while ignoring frequency and memory bandwidth.
- Do not use the optimized transform as its own only oracle.
- Do not claim legacy compatibility from matching one (K,R) pair.
- Do not integrate IT2025's unknown-error solver into erasure-only execution.
- Do not leave global mutable codec parameters like the research xdrs code.
- Do not require root-only performance tooling or mandatory third-party runtime dependencies.
- Do not hide failed experiments; record them in Beads and the final report.

20. MILESTONE ORDER

Follow this order unless a discovered dependency requires a documented adjustment.

MILESTONE 0: BOOTSTRAP AND BASELINE

- Beads initialized.
- Current code builds/tests.
- Baseline and golden vectors recorded.
- Research notes and coordinate maps drafted.

MILESTONE 1: INDEPENDENT SCALAR ORACLE

- Generic small-field/GF8/GF16 direct RS oracle.
- Active-subspace scalar LCH transforms.
- Exhaustive small-field suite.

MILESTONE 2: LEGACY HIGH ENCODER

- Refactor/preserve existing encoder.
- Establish parity compatibility.
- Add arbitrary tails in new API.
- Add streaming experiment scaffold.

MILESTONE 3: HIGH-RATE SPECIALIZED DECODER

- Scalar active-parent Algorithm 5.
- Message-only recovery.
- Legacy differential tests.
- SIMD port after scalar correctness.

MILESTONE 4: LOW PROFILE

- Low encoder.
- Algorithm 4 decoder.
- R > K.
- MDS and parity rebuild tests.

MILESTONE 5: ARBITRARY COUNTS AND PLANS

- Shortening/puncturing.
- Deterministic received selection.
- Reusable plan.
- parent-size locator setup.
- small-loss direct path.

MILESTONE 6: PRODUCTION SIMD AND MEMORY

- shared backend layer;
- SSSE3/AVX2 cleanup;
- NEON;
- pruned schedules;
- tiled scratch;
- batch API.

MILESTONE 7: 128-CORE SCALING

- scheduler;
- NUMA;
- batch scaling;
- race tests;
- all-core benchmarks.

MILESTONE 8: SPECULATIVE EXPERIMENTS

Run A-W in parallel where code dependencies allow. Promote winners behind tests and dispatch. Keep losers as concise reports, not dead production code.

MILESTONE 9: RELEASE READINESS

- all tests and sanitizers;
- compatibility gate;
- benchmark report;
- API docs;
- wire-format docs;
- migration guide;
- clean Beads graph and git tree.

21. DEFINITION OF DONE

The work is complete only when:

1. Old Leopard tests pass.
2. New API tests pass under scalar and every available SIMD backend.
3. Legacy high parity is bit-identical wherever claimed.
4. Low and high profiles recover all valid tested erasure patterns.
5. Exhaustive small-field MDS tests pass.
6. Random GF8/GF16 differential tests pass at large scale.
7. ASan/UBSan and concurrency tests are clean.
8. Fuzzers complete their allocated campaigns with no unresolved crash.
9. The 128-core correctness/experiment matrix completes and is resumable.
10. Benchmarks compare old/new/direct/high/low and report setup separately.
11. Dispatcher never knowingly chooses a slower path in its calibrated target cells.
12. No default experimental ABI or wire-format instability remains.
13. Documentation explains profiles, fields, parent inflation, preprocessing, scratch, threading, compatibility, and limitations.
14. All promoted changes have Beads evidence and local commits.
15. Rejected experiments have a short explanation and retained result data.
16. git status is clean.
17. The complete URL bibliography and literature-refresh log are committed and every implemented mathematical identity has a source or labeled derivation.
18. Encoder-specific compatibility, partial-output, streaming, and incremental-update experiments have explicit results.

22. FINAL RESPONSE FORMAT

At the end, report:

1. Branch and commit list.
2. Beads epic and closed/open bead summary.
3. Architecture and public API added.
4. Exact high/low profile definitions.
5. Compatibility result with old Leopard.
6. Correctness evidence and test counts.
7. Sanitizer/fuzzer results.
8. Benchmark tables for representative low, balanced, and high rates.
9. 1-to-128-core scaling results and NUMA findings.
10. Memory/scratch reductions.
11. Experiments promoted, rejected, and still inconclusive.
12. Known limitations and next highest-value work.
13. Exact commands to reproduce the build, tests, fuzzing, and benchmark report.
14. A source appendix containing the exact paper, repository, documentation, and literature-search URLs actually used.

Do not stop after planning. Continue implementing, testing, measuring, fixing, and documenting until the definition of done is met or a genuine external blocker is encountered. For any blocker, record it in Beads with the exact command, error, and the best available workaround, then continue all independent work.

END CODEX PROMPT

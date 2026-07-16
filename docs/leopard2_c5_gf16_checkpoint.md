# Leopard2 C5 GF16 and memory checkpoint

## Outcome

No C5 implementation is promoted.  The sampled parent-compatible dyadic-block
identities are correct, but the coarse materialized legacy-GF16 stage model
reduces logical memory traffic by only 2.78%, below the 10% promotion gate.  An ideal
ceiling that makes factor, tail-input, and accumulation traffic free reaches
11.11% at `q=257`; it is not an implementation or measurement and describes
the general C1/C2 pruned-parent fusion problem, not a distinct C5 wire route.
The naive exact-prefix route remains coupled across every sampled input/output
block pair and adds an explicit local-to-global join.  These are bounded kills
of the studied routes, not a proof that no factored or pruned exact algorithm
can exist.

The later bounded C++ executor in `docs/leopard2_c5_cpp_execution.md` completes
the remaining SIMD, sanitizer, tail, setup/execution, scratch, batch/reuse, and
pinned timing gates.  Its strongest fused `q=2^a+1` implementation also fails
the region-wide 10% promotion rule, so the disposition here is unchanged and
the previously open measurement obligation is resolved.

The deterministic program and artifact are:

    experiments/leopard2/c5_gf16_checkpoint.py
    experiments/leopard2/c5_gf16_checkpoint_result.json

This checkpoint extends the scalar C5 result in
`experiments/leopard2/non_power_of_two/c5/`; it does not modify the codec,
public API, dispatcher, ABI, or wire format.

## Two constructions, two compatibility claims

`parent_wire_block_accumulation` splits an active LCH coefficient prefix into
canonical aligned blocks and accumulates their contributions to the existing
power-of-two parent.  It computes the same mathematical map and does not define
a new code identity.  The simulator verifies representative instances of

    X_(o+t)(x) = X_o(x) X_t(x)

for every selected GF16 prefix, where block offset `o` is aligned to block size
and `t` is a local coefficient index.

`exact_prefix_join` interpolates exactly `K=q` source points and evaluates 17
following parity points.  That coordinate set is a new experimental exact
profile unless its full generator matrix is proved equal to an existing
profile.  Nothing in this checkpoint claims parent-wire or legacy parity
compatibility for it.

## Executable legacy-GF16 algebra

The simulator uses Leopard's polynomial `0x1002d`, 16-element Cantor basis,
coordinate order, and symbol values.  Multiplication converts additive
coordinates through independently generated polynomial-basis log/exp tables.
All 16 recurrence-generated normalization values `s_j(v_j)` are compared with
direct products over `V_j`, including the 32,768-factor final product.

The retained validation covers:

| Check | Count |
| --- | ---: |
| Direct GF16 normalizer-product comparisons | 16 |
| Deterministic multiply/inverse checks | 36 |
| Log multiply versus carryless shift/reduce | 36 |
| Non-basis subspace recurrence versus direct products | 32 |
| Parent block-identity samples | 930 |
| Parent coset factors evaluated | 391,871 |
| Exact GF16 join matrices | 8 |
| Dense global parity coefficients | 17,170 / 17,170 |
| Effective local-to-parity nonzeros | 14,604 / 17,170 |
| Coupled input/output dyadic block pairs | 92 / 92 |

For each exact case, the LCH local-to-global join followed by LCH parity
evaluation is compared entry-by-entry with an independently constructed direct
Lagrange map applied to the local block evaluations.

Representative join density is:

| q | Blocks | Join nonzeros / entries | Cross-block nonzeros / entries | Effective parity nonzeros / entries | Coupled pairs |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 17 | 2 | 33 / 289 | 16 / 32 | 257 / 289 | 10 / 10 |
| 33 | 2 | 65 / 1,089 | 32 / 64 | 495 / 561 | 10 / 10 |
| 63 | 6 | 321 / 3,969 | 232 / 2,604 | 1,013 / 1,071 | 12 / 12 |
| 65 | 2 | 129 / 4,225 | 64 / 128 | 973 / 1,105 | 10 / 10 |
| 129 | 2 | 257 / 16,641 | 128 / 256 | 1,913 / 2,193 | 10 / 10 |
| 191 | 7 | 959 / 36,481 | 711 / 18,732 | 2,056 / 3,247 | 14 / 14 |
| 255 | 8 | 1,793 / 65,025 | 1,418 / 43,180 | 4,088 / 4,335 | 16 / 16 |
| 257 | 2 | 513 / 66,049 | 256 / 512 | 3,809 / 4,369 | 10 / 10 |

The joins immediately above powers of two are sparse and triangular, but the
effective parity maps remain materially populated and every dyadic input block
affects every requested output block.  Multi-block prefixes add substantially more
cross-block structure.  Density and coupling do not rule out a fast
factorization; they do reject treating independent block transforms as if they
were already an exact encoder.

## Parent-compatible operation and scratch evidence

The parent sweep includes 18 GF16 prefixes from 257 through 65,535.  A logical
butterfly reads and writes two shards.  A fixed multiplication reads and writes
one shard, and an accumulation XOR reads source and destination then writes the
destination.  Both materialized candidate and padded baseline are optimistic
transform-stage models: they exclude repeated input staging, final output
placement, packing, scatter, allocator effects, and cache-line effects.  They
are not executable end-to-end memory totals.  A second, deliberately
unattainable fusion bound charges only butterfly traffic and makes all other
work free.

| q | Blocks | Candidate / padded butterflies | Factor muls | Join XORs | Materialized improvement | Free-fusion ceiling | Scratch slots candidate / padded |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 257 | 2 | 2,048 / 2,304 | 0 | 256 | 2.78% | 11.11% | 256 / 512 |
| 513 | 2 | 4,608 / 5,120 | 0 | 512 | 2.50% | 10.00% | 512 / 1,024 |
| 1,000 | 6 | 12,032 / 5,120 | 1,792 | 2,560 | -190.00% | -135.00% | 512 / 1,024 |
| 4,096 | 1 | 24,576 / 24,576 | 0 | 0 | 0.00% | 0.00% | 4,096 / 4,096 |
| 4,097 | 2 | 49,152 / 53,248 | 0 | 4,096 | 1.92% | 7.69% | 4,096 / 8,192 |
| 16,385 | 2 | 229,376 / 245,760 | 0 | 16,384 | 1.67% | 6.67% | 16,384 / 32,768 |
| 32,769 | 2 | 491,520 / 524,288 | 0 | 32,768 | 1.56% | 6.25% | 32,768 / 65,536 |
| 65,535 | 16 | 2,211,840 / 524,288 | 442,368 | 491,520 | -434.38% | -321.88% | 32,768 / 65,536 |

The `2^a+1` family halves local scratch, but its materialized traffic advantage
shrinks with transform depth.  Only the free-fusion ceiling for `q=257` exceeds
10%; `q=513` merely equals it.  That ceiling omits the tail read and every
accumulation cost, so it is neither a promotion result nor a performance claim.
Any implementable fusion or cache benefit computes the existing parent map and
belongs in the more general C1/C2 pruning and scheduling work rather than a
distinct C5 wire path.  Prefixes with many set bits regress sharply even under
the free-fusion bound.

## Byte, batch, and plan-reuse matrix

The artifact contains 1,440 cells over:

    shard bytes: 64, 1 KiB, 64 KiB, 1 MiB, 16 MiB
    batch:       1, 8, 64, 1024
    plan reuse:  1, 8, 64, 1024

Modeled transform-stage traffic scales exactly by `shard_bytes * batch * reuse`.  Candidate
setup is modeled once as a 32-byte header plus 16 bytes per block and active
coset, then amortized over `batch * reuse` stripes.  This is an explicit record
layout estimate, not measured allocator or cache traffic.

Each cell also records scratch bytes per stripe and the conservative case where
all batch scratch is simultaneously resident.  For `q=257` and 1 MiB shards,
this is 256 MiB per candidate stripe versus 512 MiB for the padded parent; an
eight-stripe resident batch is 2 GiB versus 4 GiB.  This capacity benefit does
not turn into a 10% logical-traffic win.

Selected cells are:

| q | Bytes | Batch | Reuse | Materialized improvement | Free-fusion ceiling | Setup bytes per stripe |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 257 | 64 | 1 | 1 | 2.07% | 10.40% | 4,192 |
| 257 | 1 KiB | 8 | 8 | 2.78% | 11.11% | 65.5 |
| 257 | 16 MiB | 1,024 | 1,024 | 2.78% | 11.11% | 0.004 |
| 513 | 64 | 1 | 1 | 1.87% | 9.37% | 8,288 |
| 4,097 | 64 KiB | 64 | 64 | 1.92% | 7.69% | 16.02 |
| 65,535 | 64 | 1 | 1 | -435.16% | -322.66% | 1,048,864 |
| 65,535 | 16 MiB | 1,024 | 1,024 | -434.38% | -321.88% | 1.0003 |

The best materialized cell remains 2.78%; batch and reuse remove setup overhead
but cannot change the execution ratio.  The 11.11% ceiling is explicitly
unimplemented and unmeasured.  No runtime timing is quoted because Python
matrix time is not codec evidence and the only potentially interesting fused
parent schedule is assigned to C1/C2.

## Disposition

- Kill the coarse materialized parent-wire block schedule as a distinct
  implementation candidate.  Preserve the
  algebraic identity and free-fusion ceiling as evidence for C1/C2 pruned or
  tiled parent schedules; neither is a distinct C5 production route.
- Kill explicit same-basis exact join followed by parent evaluation.  It adds
  join storage and execution before an already dense parity map.
- Do not infer that all exact codes are slow.  A different factored exact
  construction remains inconclusive and belongs to C6/C7/C8, with a new profile
  identifier and independent MDS proof.
- Promote no code, profile, dispatcher cell, or wire-format change from C5.

This scalar checkpoint did not itself satisfy the full C5 C++/SIMD/timing
gate.  The only modeled ceiling above the promotion threshold required fused
parent pruning.  The later C++ checkpoint implements and rejects the bounded
`2^a+1` fusion; broader general-purpose pruning remains a separate C1/C2
scheduler question.

## Reproduction

From the repository root, use all available work up to the 26 independent jobs:

    PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0 python3 -X dev \
        experiments/leopard2/c5_gf16_checkpoint.py \
        --workers 26 \
        --output /tmp/c5-gf16-parallel.json

    PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0 python3 -X dev \
        experiments/leopard2/c5_gf16_checkpoint.py \
        --workers 1 \
        --output /tmp/c5-gf16-serial.json

    cmp /tmp/c5-gf16-parallel.json /tmp/c5-gf16-serial.json
    cmp /tmp/c5-gf16-parallel.json \
        experiments/leopard2/c5_gf16_checkpoint_result.json

The worker count is deliberately excluded from the JSON identity.  Parallel
and serial generation must be byte-identical.

The following independent seeded check evaluates the direct `q`-term LCH
polynomial and the dyadic-block accumulation separately for ten irregular
prefixes.  It covers 1,092 enclosing-parent output points and must print
`independent direct-vs-block output points 1092 PASS`:

    PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0 python3 -X dev - <<'PY'
    import importlib.util
    import pathlib
    import random

    path = pathlib.Path('experiments/leopard2/c5_gf16_checkpoint.py')
    spec = importlib.util.spec_from_file_location('leopard2_c5_checkpoint', path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    field = module.make_gf16()
    factors = module.normalizer_factors(field)
    rng = random.Random(0xC5)
    checks = 0

    for q in (3, 5, 7, 9, 17, 33, 63, 65, 129, 257):
        parent = module.ceil_power_of_two(q)
        coefficients = [rng.randrange(65536) for _ in range(q)]
        direct = []
        for point in range(parent):
            values = module.lch_values(field, factors, point, q)
            value = 0
            for coefficient, basis_value in zip(coefficients, values):
                value ^= field.multiply(coefficient, basis_value)
            direct.append(value)

        candidate = [0] * parent
        for offset, size in module.binary_prefix_blocks(q):
            for coset in range(0, parent, size):
                factor = module.lch_value(field, factors, coset, offset)
                if not factor:
                    continue
                for local_point in range(size):
                    point = coset + local_point
                    values = module.lch_values(field, factors, point, size)
                    value = 0
                    for index, basis_value in enumerate(values):
                        value ^= field.multiply(
                            coefficients[offset + index], basis_value
                        )
                    candidate[point] ^= field.multiply(factor, value)

        if direct != candidate:
            raise SystemExit(f'direct-vs-block mismatch q={q}')
        checks += parent

    print('independent direct-vs-block output points', checks, 'PASS')
    PY

Retained SHA-256 values for this checkpoint are:

    781b8b3cc5d0b5a8e2843eabb2e3c949a9e07ed901d6dd50a872401371826d70  c5_gf16_checkpoint.py
    d577538fc7962fb2e66ac2eb01a9027d5a124c87b0393a32c1e145a9f8837cfb  c5_gf16_checkpoint_result.json

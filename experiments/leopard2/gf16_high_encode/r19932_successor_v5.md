# AUTO R199/32 KiB successor v5 (inconclusive, 2026-09-15)

Beads: `leopard-79h.38.5.4.19.1.6`. This is a distinct successor to the
consumed v3 attempt and the method-invalid v4 attempt. The readiness-true
preregistration was pushed at commit `cf67a257fbe42b844bf29f8b588c191ae9d18040`
before the one timing attempt. No previous attempt was retried or overwritten.

The run used the existing qualified default-off candidate at commit
`45e2effd869859c9b3aa48190eff6f4738817c61`, exact native Leopard1 at
`6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198`, CPU 26/sibling 90, and a
256 MiB/no-swap systemd scope. It completed 27 preflights and 318 timed
processes with zero sibling jiffies. The captured resource footer reports a
184,995,840-byte peak, `memory.max=268435456`, zero memory events, and zero
swap. The raw attempt is retained in the local read-only
`.research/leopard-79h/paired-epoch-campaign-v5-attempt.20260915` bundle.

## Fixed decision

The preregistered thresholds remain unchanged: every target must exceed 1.05×
versus OFF and native Leopard1 in every round, while all controls and neighbors
must stay within ±2%. Both target APIs pass their gain requirements in every
epoch, but fixed controls fail in all three epochs, so this is
`inconclusive_controls`; the route remains default-off and no production
promotion is authorized.

| Epoch | Ordinary OFF | Batch OFF | Ordinary native | Batch native |
| ---: | ---: | ---: | ---: | ---: |
| 0 | 1.518× | 1.525× | 1.424× | 1.414× |
| 1 | 1.521× | 1.524× | 1.425× | 1.417× |
| 2 | 1.524× | 1.530× | 1.417× | 1.411× |

The aggregate GF8 K17/R7/B64 control is outside the fixed interval in every
epoch (for example 1.040×/0.953× in epoch 0 and 1.021×/0.980× in epoch 2),
and the per-round native control also crosses the bound. These are retained,
not trimmed or pooled. The result is evidence of a stable large target gain
under this candidate, but not valid permission to enable AUTO R19932.

Raw evidence hashes:

| File | SHA-256 |
| --- | --- |
| `resource.log` | `06e9f980a89986894908d5763f4af2385e9b53bb9dd9daf4deb7537aca1f4812` |
| `attempt1/attempt.json` | `b6a58f5b2e0bfebd56eefeb4f707b6b2ce58e42dd8059911d5792c76113c6fc4` |
| `attempt1/launches.jsonl` | `0d596dd3847948b6e34af30ec00666f13cff9d7eb7eae5cf1ebea89116658638` |

The original v3 inconclusive attempt and v4 missing-footer attempt remain
separately retained. A release can continue with this route disabled and the
limitation documented; no threshold relaxation or timing retry is justified.

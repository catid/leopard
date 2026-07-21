# R=1 exact Leopard-main gate

This checkpoint compares production Leopard2 commit
`e6211091784df3e0c746790c7650450d8528119f` with the independently linked
Leopard default-branch commit
`6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198`.  The rebuilt Leopard2 executable
and archive are byte-identical to the earlier e621109 production artifacts.
Test hooks, fuzzers, and CUDA were disabled.

Ratios are exact-main time divided by Leopard2 time, so values above one favor
Leopard2.  Each row is a separate three-round ABBA campaign with nine retained
samples per process, two warmups, reuse eight, one thread, and exact parity and
recovery digest matching.  The reserved SMT sibling accumulated zero non-idle
jiffies in every accepted campaign.

| K, bytes | Encode, 95% CI | Decode first use, 95% CI | Decode reuse 8, 95% CI |
| --- | ---: | ---: | ---: |
| 128, 4 KiB | 1.051 [1.005, 1.100] | 1.009 [0.982, 1.037] | 1.047 [1.020, 1.076] |
| 129, 4 KiB | 1.104 [1.024, 1.189] | 1.013 [1.008, 1.019] | 1.052 [1.046, 1.059] |
| 130, 4 KiB | 1.081 [1.054, 1.110] | 1.017 [0.995, 1.040] | 1.055 [1.032, 1.077] |
| 129, 64 KiB | 1.329 [1.176, 1.503] | 1.183 [1.076, 1.299] | 1.185 [1.078, 1.302] |
| 129, 1 MiB | 1.098 [1.085, 1.112] | 1.103 [1.095, 1.111] | 1.103 [1.095, 1.111] |

The clean K=130 boundary result disproves the apparent decode regression in
the earlier contaminated multi-cell run.  No selector adjustment is needed.
The coarse R=1 XOR path has credible gains over exact main at the larger shard
sizes, and the 4 KiB neighbors show no credible regression over two percent.

The signed manifests are committed beside this file.  They retain the complete
source/build/runtime identities, statistical summaries, raw-bundle hashes, and
isolation records.  The raw streams remain in the ignored experiment cache.
`summary.json` records both rejected attempts: one pre-timing stale-lease
failure and one 1 MiB run rejected because its reserved sibling performed two
user jiffies.  Neither rejected timing contributes to this table.

To verify a complete retained bundle while the ignored raw directory is still
available, run:

    python3 -I -S experiments/leopard2/main_compare/run_abba.py verify \
        --manifest /tmp/leopard-r1-authoritative-e621109/k129-64k/manifest.json

# AUTO R199/32 KiB: disabled after inconclusive controls

The `K=1000,R=199`, 32-KiB AUTO GFNI extension remains **default-off**.
The 2026-09-15 successor-v5 diagnostic observed large target gains, but failed
its preregistered stability controls. These are diagnostic ratios, not
qualified release speedups or permission to enable the route.

The fixed gates required every target to exceed 1.05× versus both the disabled
route and native Leopard1 in every round, with controls and neighbors within
±2%. Both ordinary and one-item batch targets passed their gain requirements;
controls failed in all three epochs.

| Epoch | Ordinary vs OFF | Batch vs OFF | Ordinary vs native | Batch vs native |
| ---: | ---: | ---: | ---: | ---: |
| 0 | 1.518× | 1.525× | 1.424× | 1.414× |
| 1 | 1.521× | 1.524× | 1.425× | 1.417× |
| 2 | 1.524× | 1.530× | 1.417× | 1.411× |

The GF8 `K=17,R=7,B=64` aggregate control fell outside its interval in every
epoch, and a per-round native control also crossed its bound. No attempts
were pooled, trimmed, or retried, and no threshold was relaxed.

The diagnostic used default-off candidate
`45e2effd869859c9b3aa48190eff6f4738817c61` and native Leopard1
`6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198`, after preregistration
`cf67a257fbe42b844bf29f8b588c191ae9d18040` was pushed. Its 27 preflights and
318 timed processes ran on CPU 26 with zero activity on SMT sibling 90.
Peak memory was 184,995,840 bytes under a 256-MiB/no-swap limit, with zero
memory events and swap use.

The [full research report](https://github.com/catid/leopard/blob/986c922/experiments/leopard2/gf16_high_encode/r19932_successor_v5.md)
records retained evidence identities. Raw local research bundles are not part
of the source distribution. Earlier failed attempts remain separately
retained; this result does not justify another timing retry. The already
qualified R200/32-KiB, R200/64-KiB, and R199/64-KiB routes are unchanged.

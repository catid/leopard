# 01 — GFNI affine multiplication

**Disposition: LANDED** (behind `LEO2_GFNI_VARIANT`, not yet promoted to a
production backend member)

## Idea
A fixed GF multiply is a GF(2)-linear map, so each byte-wide component is an
8x8 bit matrix and `VGF2P8AFFINEQB` evaluates it in one instruction. GF8: one
instruction replaces six. GF16: the 16x16 matrix splits into four 8x8 blocks
over Leopard's split low/high byte layout, so four affine ops plus two XORs
replace four ANDs, two shifts, eight shuffles and eight XORs.

## Versus leopard1 (exact main), encode
| Cell | stock Leopard2 | with GFNI |
| --- | ---: | ---: |
| GF8 K=128 R=128, 4 KiB | 0.972 | **1.549** |
| GF8 K=128 R=128, 64 KiB | 1.070 | **1.598** |
| GF8 K=240 R=16, 4 KiB | 1.168 | **1.540** |
| GF16 K=255 R=129, 64 KiB | 0.946 | **1.240** |
| GF16 K=1000 R=200, 4 KiB | 1.091 | **1.555** |
| GF16 K=4096 R=512, 4 KiB | 1.119 | **1.462** |
| GF8 K=2 R=2, 1 KiB | 0.574 | 0.577 |

Every cell that lost to leopard1 flips except `K=2 R=2` and `K=4 R=4` at 1 KiB,
which are fixed per-call cost, not transform cost.

## Why it worked
Kernel screen: GF8 butterfly 1.32-2.08x, GF16 multiply-add 1.33-2.47x against
the nibble form. It reduces operation count on an ALU-bound kernel.

## Caveats
GFNI is Intel Ice Lake+/Tremont+ and AMD Zen 4/5 only. Haswell..Skylake and
Zen 1-3 get nothing. The evaluation build fails `leopard2_portable_isa` by
design — a production member needs its own qualified translation unit.

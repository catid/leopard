# GF16 Walsh locator AVX2 result

The same-process successor screen preregistered in
[`gf16_walsh_locator_avx2_preregistration_v2.md`](gf16_walsh_locator_avx2_preregistration_v2.md)
passed all six cells. The qualified AVX2 setup path was 3.902×–18.596× faster
than the scalar active-parent setup, with zero reserved-sibling non-idle ticks.

| Parent | Erasures | Scalar median | AVX2 median | Speedup |
|---:|---:|---:|---:|---:|
| 32 | 32 | 0.908875 µs | 0.048875 µs | 18.596× |
| 64 | 40 | 1.045250 µs | 0.082625 µs | 12.651× |
| 256 | 64 | 1.862750 µs | 0.240375 µs | 7.749× |
| 1,024 | 64 | 5.092625 µs | 1.034000 µs | 4.925× |
| 4,096 | 256 | 19.104750 µs | 4.629375 µs | 4.127× |
| 65,536 | 8,192 | 336.918750 µs | 86.337500 µs | 3.902× |

The complete provenance and resource record is in
[`gf16_walsh_locator_avx2_v2.json`](gf16_walsh_locator_avx2_v2.json). This is
setup-only evidence; it does not claim an end-to-end Leopard1 throughput gain,
and no AUTO selector was changed.

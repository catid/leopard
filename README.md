# Leopard-RS

Windows users should use the CMake workflow described in
[`docs/leopard2_windows_build.md`](docs/leopard2_windows_build.md).  The
checked-in Visual Studio 2015-format solution remains available for legacy
consumers and is structurally checked against the production Leopard2 source
graph.  Native Visual Studio 2015 load/build validation is still outstanding.

## Portable builds

The default CMake build is portable for the target platform; it does not add
`-march=native`.  On x86-64, baseline code stays at the platform's SSE2 floor,
while SSSE3, AVX2, and AVX-512VL kernels are compiled in
separate translation units and selected only after runtime CPU and
operating-system checks.  An `AUTO` context reports the established AVX2
baseline when it is available.  On the offline-calibrated AMD family
1Ah/model 44h host class, a legacy-high GF16 codec may use the qualified
AVX-512VL table for a full-output encode only when `K >= 8`, `N >= 16`,
`2 <= R <= 4096`, and the shard length is a multiple of 64 bytes from 64 bytes
through 4 MiB inclusive.  All other AUTO operations retain the context
baseline; explicit backend requests retain their exact table.  A normal
distributable build therefore needs no host-CPU option:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target leopard
```

Do not add `-march=native` when producing binaries for other machines: it can
allow the compiler to emit unsupported instructions anywhere in the program.
For cross-compilation, select the destination architecture in a CMake toolchain
file and keep its baseline compatible with every deployment CPU.  The
Apple build specifically treats a single `CMAKE_OSX_ARCHITECTURES` entry as the
target architecture rather than trusting the build host's processor.  Universal
Apple archives deliberately omit the optional x86 SSSE3/AVX object libraries:
one object-library compile option cannot safely describe both its x86 and Arm
slices, so those builds retain the portable scalar/native-Arm graph.

The `LEO2_BACKEND_VARIANT` setting is a diagnostic control for forcing a
qualified backend during testing; accepted values are `auto`, `scalar`,
`ssse3`, `avx2`, and `avx512`.  Lower forced variants cap explicit contexts at
the selected backend, while `avx512` retains access to qualified lower tables.
This setting is not a portability target or a wire-format choice.  `auto`
continues to report AVX2 as its context baseline even when one bounded encode
operation is eligible for the model-scoped AVX-512VL policy above.
When its POSIX shell and disassembly tools are available, the normal test build
registers `leopard2_portable_isa` to audit the x86-64 archive and any available
build metadata.  Missing audit tools remain non-fatal for ordinary developer,
cross-platform, and dependency-light builds.

Sanitizer instrumentation is not production ISA evidence: compiler-generated
ASan/UBSan helpers may use instructions outside the codec object's declared
floor.  A sanitizer-instrumented build therefore registers
`leopard2_portable_isa_checker_self_test` to retain all adversarial classifier
controls, while the archive audit remains exclusive to a clean build.  The
classification follows the active custom build type and its
`CMAKE_CXX_FLAGS_<CONFIG>` variable.  Multi-configuration generators classify
each configuration independently: sanitized configurations run the checker
self-test, while clean configurations retain the archive audit (select one with
CTest's `-C` option).  An empty single-config build type remains valid and is
classified from the general C++ flags.

The strict release-audit option below rejects sanitizer flags instead of
silently substituting the weaker mode.

```sh
ctest --test-dir build -R '^leopard2_portable_isa$' --output-on-failure
```

Release CI can make every part of that evidence mandatory with the default-off
strict mode:

```sh
cmake -S . -B build/release-audit -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
  -DLEO2_PORTABLE_ISA_RELEASE_AUDIT=ON
cmake --build build/release-audit --target leopard
ctest --test-dir build/release-audit \
  -R '^leopard2_portable_isa$' --output-on-failure
```

This mode fails closed if tests, the static-archive tool, `objdump` or
`llvm-objdump`, POSIX `sh`, `compile_commands.json`, or the archive audit itself
is unavailable.  It currently supports native, single-configuration x86-64
release builds on non-Apple Unix with GCC or Clang.  MSVC/Windows uses the
portable runtime-dispatch build plus the separately documented Visual Studio
structural audit; enabling this strict ELF/GNU-archive audit there is an
intentional configuration error.

## Leopard2 field selection and reduced builds

The new `leopard2.h` API lets an application select `LEO2_FIELD_GF8` or
`LEO2_FIELD_GF16` explicitly when it creates a codec. `LEO2_FIELD_AUTO` is a
wire-stable convenience, not a performance tuner: it uses GF8 when the selected
profile's power-of-two parent has at most 256 coordinates and GF16 otherwise.
That mapping does not depend on the CPU or on which fields were compiled. In
particular, explicitly choosing GF16 for a small code changes the persisted
field/code identity, and an AUTO codec whose canonical field was omitted
returns `LEO2_UNSUPPORTED` instead of silently emitting different parity.

Both fields are included by default. Applications which only need one may
reduce code and table footprint with:

```sh
cmake -S . -B build-gf8 -DLEOPARD_ENABLE_GF16=OFF  # GF8 only
cmake -S . -B build-gf16 -DLEOPARD_ENABLE_GF8=OFF  # GF16 only
```

At least one field must remain enabled. Non-CMake builds can use the equivalent
`NO_LEO_HAS_FF8` and `NO_LEO_HAS_FF16` preprocessor controls. The installed
CMake package reports `leopard_GF8_ENABLED` and `leopard_GF16_ENABLED`, while
`leo2_context_field_mask()` reports the fields usable by a successfully created
context. The legacy `leo_*` API retains its historical canonical field choice;
if that field is absent, a nontrivial operation returns `Leopard_Platform`
rather than switching wire formats.

GF8 has much smaller initialization tables and accepts arbitrary positive shard
byte lengths. GF16 supports larger parent codes and its native layout requires
complete two-byte symbols. Passing an odd physical `shard_bytes` value to a
native explicit-GF16 codec deterministically returns `LEO2_UNSUPPORTED`; it is
never silently padded or projected to one byte. For an odd application payload
of `B` bytes, explicitly select `LEO2_SHARD_LAYOUT_GF16_PADDED_ODD_V1`, query the
physical `B+1` size, and retain that extra byte in every parity shard. Relative
throughput depends on the code shape, shard size, backend, and target CPU;
neither field is universally faster.

## MDS Reed-Solomon Erasure Correction Codes for Large Data in C

Leopard-RS is a fast library for Erasure Correction Coding.
From a block of equally sized original data pieces, it generates recovery
symbols that can be used to recover lost original data.


#### Motivation:

It gets slower as O(N Log N) in the input data size, and its inner loops are
vectorized using the best approaches available on modern processors, using the
fastest finite fields (8-bit or 16-bit Galois fields with Cantor basis {2}).

It sets new speed records for MDS encoding and decoding of large data,
achieving over 1.2 GB/s to encode with the AVX2 instruction set on a single core.

Example applications are data recovery software and data center replication.


#### Encoder API:

Preconditions:

* The original and recovery data must not exceed 65536 pieces.
* The recovery_count <= original_count.
* The buffer_bytes must be a multiple of 64.
* Each buffer should have the same number of bytes.
* Even the last piece must be rounded up to the block size.

```
#include "leopard.h"
```

For full documentation please read `leopard.h`.

+ `leo_init()` : Initialize library.
+ `leo_encode_work_count()` : Calculate the number of work_data buffers to provide to leo_encode().
+ `leo_encode()`: Generate recovery data.


#### Decoder API:

For full documentation please read `leopard.h`.

+ `leo_init()` : Initialize library.
+ `leo_decode_work_count()` : Calculate the number of work_data buffers to provide to leo_decode().
+ `leo_decode()` : Recover original data.


#### Benchmarks:

On the the MacBook Pro 15-Inch (Mid-2015) featuring a 22 nm "Haswell/Crystalwell" 2.8 GHz Intel Core i7-4980HQ processor, compiled with Visual Studio 2017 RC:

~~~
Leopard Encoder(8.192 MB in 128 pieces, 128 losses): Input=2102.67 MB/s, Output=2102.67 MB/s
Leopard Decoder(8.192 MB in 128 pieces, 128 losses): Input=686.212 MB/s, Output=686.212 MB/s

Leopard Encoder(64 MB in 1000 pieces, 200 losses): Input=2194.94 MB/s, Output=438.988 MB/s
Leopard Decoder(64 MB in 1000 pieces, 200 losses): Input=455.633 MB/s, Output=91.1265 MB/s

Leopard Encoder(2097.15 MB in 32768 pieces, 32768 losses): Input=451.168 MB/s, Output=451.168 MB/s
Leopard Decoder(2097.15 MB in 32768 pieces, 32768 losses): Input=190.471 MB/s, Output=190.471 MB/s
~~~

2 GB of 64 KB pieces encoded in 4.6 seconds, and worst-case recovery in 11 seconds.

More benchmark results are available here:
[https://github.com/catid/leopard/blob/master/Benchmarks.md](https://github.com/catid/leopard/blob/master/Benchmarks.md)


#### Comparisons:

There is another library `FastECC` by Bulat-Ziganshin that should have similar performance:
[https://github.com/Bulat-Ziganshin/FastECC](https://github.com/Bulat-Ziganshin/FastECC).
Both libraries implement the same high-level algorithm in {3}, while Leopard implements the
newer polynomial basis GF(2^r) approach outlined in {1}, and FastECC uses complex finite fields
modulo special primes.  There are trade-offs that may make either approach preferable based
on the application:
+ Older processors do not support SSSE3 and FastECC supports these processors better.
+ FastECC supports data sets above 65,536 pieces as it uses 32-bit finite field math.  
+ Leopard does not require expanding the input or output data to make it fit in the field, so it can be more space efficient.


#### FFT Data Layout:

We pack the data into memory in this order:

~~~
[Recovery Data (Power of Two = M)] [Original Data] [Zero Padding out to 65536]
~~~

For encoding, the placement is implied instead of actual memory layout.
For decoding, the layout is explicitly used.


#### Encoder algorithm:

The encoder is described in {3}.  Operations are done O(K Log M),
where K is the original data size, and M is up to twice the
size of the recovery set.

Roughly in brief:

~~~
Recovery = FFT( IFFT(Data_0) xor IFFT(Data_1) xor ... )
~~~

It walks the original data M chunks at a time performing the IFFT.
Each IFFT intermediate result is XORed together into the first M chunks of
the data layout.  Finally the FFT is performed.


Encoder optimizations:
* The first IFFT can be performed directly in the first M chunks.
* The zero padding can be skipped while performing the final IFFT.
Unrolling is used in the code to accomplish both these optimizations.
* The final FFT can be truncated also if recovery set is not a power of 2.
It is easy to truncate the FFT by ending the inner loop early.
* The decimation-in-time (DIT) FFT is employed to calculate two layers at a time, rather than writing each layer out and reading it back in for the next layer of the FFT.


#### Decoder algorithm:

The decoder is described in {1}.  Operations are done O(N Log N), where N is up
to twice the size of the original data as described below.

Roughly in brief:

~~~
Original = -ErrLocator * FFT( Derivative( IFFT( ErrLocator * ReceivedData ) ) )
~~~


#### Precalculations:

At startup initialization, FFTInitialize() precalculates FWT(L) as
described by equation (92) in {1}, where L = Log[i] for i = 0..Order,
Order = 256 or 65536 for FF8/16.  This is stored in the LogWalsh vector.

It also precalculates the FFT skew factors (s_i) as described by
equation (28).  This is stored in the FFTSkew vector.

For memory workspace N data chunks are needed, where N is a power of two
at or above M + K.  K is the original data size and M is the next power
of two above the recovery data size.  For example for K = 200 pieces of
data and 10% redundancy, there are 20 redundant pieces, which rounds up
to 32 = M.  M + K = 232 pieces, so N rounds up to 256.


#### Online calculations:

At runtime, the error locator polynomial is evaluated using the
Fast Walsh-Hadamard transform as described in {1} equation (92).

At runtime the data is explicit laid out in workspace memory like this:
~~~
[Recovery Data (Power of Two = M)] [Original Data (K)] [Zero Padding out to N]
~~~

Data that was lost is replaced with zeroes.
Data that was received, including recovery data, is multiplied by the error
locator polynomial as it is copied into the workspace.

The IFFT is applied to the entire workspace of N chunks.
Since the IFFT starts with pairs of inputs and doubles in width at each
iteration, the IFFT is optimized by skipping zero padding at the end until
it starts mixing with non-zero data.

The formal derivative is applied to the entire workspace of N chunks.

The FFT is applied to the entire workspace of N chunks.
The FFT is optimized by only performing intermediate calculations required
to recover lost data.  Since it starts wide and ends up working on adjacent
pairs, at some point the intermediate results are not needed for data that
will not be read by the application.  This optimization is implemented by
the ErrorBitfield class.

Finally, only recovered data is multiplied by the negative of the
error locator polynomial as it is copied into the front of the
workspace for the application to retrieve.


#### Finite field arithmetic optimizations:

For faster finite field multiplication, large tables are precomputed and
applied during encoding/decoding on 64 bytes of data at a time using
SSSE3 or AVX2 vector instructions and the ALTMAP approach from Jerasure.

Addition in this finite field is XOR, and a vectorized memory XOR routine
is also used.


#### References:

This library implements an MDS erasure code introduced in this paper:

~~~
    {1} S.-J. Lin, T. Y. Al-Naffouri, Y. S. Han, and W.-H. Chung,
    "Novel Polynomial Basis with Fast Fourier Transform
	and Its Application to Reed-Solomon Erasure Codes"
    IEEE Trans. on Information Theory, pp. 6284-6299, November, 2016.
~~~

~~~
    {2} D. G. Cantor, "On arithmetical algorithms over finite fields",
    Journal of Combinatorial Theory, Series A, vol. 50, no. 2, pp. 285-300, 1989.
~~~

~~~
    {3} Sian-Jheng Lin, Wei-Ho Chung, "An Efficient (n, k) Information
    Dispersal Algorithm for High Code Rate System over Fermat Fields,"
    IEEE Commun. Lett., vol.16, no.12, pp. 2036-2039, Dec. 2012.
~~~

~~~
    {4} Plank, J. S., Greenan, K. M., Miller, E. L., "Screaming fast Galois Field
    arithmetic using Intel SIMD instructions."  In: FAST-2013: 11th Usenix
    Conference on File and Storage Technologies, San Jose, 2013
~~~
	
Some papers are mirrored in the /docs/ folder.


#### Credits

Inspired by discussion with:

+ Sian-Jhen Lin <sjhenglin@gmail.com> : Author of {1} {3}, basis for Leopard
+ Bulat Ziganshin <bulat.ziganshin@gmail.com> : Author of FastECC
+ Yutaka Sawada <tenfon@outlook.jp> : Author of MultiPar

Software by Christopher A. Taylor <mrcatid@gmail.com>

Please reach out if you need support or would like to collaborate on a project.

# Leopard-RS
## MDS Reed-Solomon Error Correction Codes for Large Data in C

#### Update: Data up to 256 pieces is working!  Next: Implementing 16-bit finite fields to enable data up to 65536 pieces.

Leopard-RS is a fast library for Forward Error Correction.
From a block of equally sized original data pieces, it generates recovery
symbols that can be used to recover lost original data.


#### Motivation:

It gets slower as O(N Log N) in the input data size, and its inner loops are
vectorized using the best approaches available on modern processors, using the
fastest finite fields (8-bit or 16-bit Galois fields with Cantor basis {2}).

It sets new speed records for MDS encoding and decoding of large data,
achieving over 1.2 GB/s to encode with the AVX2 instruction set on a single core.

There is another library `FastECC` by Bulat-Ziganshin that should have similar performance:
[https://github.com/Bulat-Ziganshin/FastECC](https://github.com/Bulat-Ziganshin/FastECC).
Both libraries implement the same high-level algorithm in {3}, while Leopard implements the
newer polynomial basis GF(2^r) approach outlined in {1}, and FastECC uses complex finite fields
modulo special primes.  There are trade-offs that may make either approach preferable based
on the application:
+ Older processors do not support SSSE3 and FastECC supports these processors better.
+ FastECC supports data sets above 64,000 pieces as it uses 32-bit finite field math.  
+ Leopard does not require expanding the input or output data to make it fit in the field, so it can be more space efficient.

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

On my laptop:

```
Leopard Encoder(8.192 MB in 128 pieces, 128 losses): Input=1266.13 MB/s, Output=1266.13 MB/s
Leopard Decoder(8.192 MB in 128 pieces, 128 losses): Input=482.243 MB/s, Output=482.243 MB/s
```


#### Comparisons:

Comparing performance from all my error correction code libraries, on my laptop:

To summarize, a set of 128 of 64 KB data files are supplemented by about 128 redundant code pieces (encoded) meaning a code rate of 1/2.  From those redundant code pieces the original set is recovered (decoded).

The results are all from libraries I've written over the past few years.  They all have the same vector-optimized inner loops, but the types of error correction codes are different.

```
For 64KB data chunks:

CM256 Encoder: 64000 bytes k = 128 m = 128 : 82194.7 usec, 99.6658 MBps
CM256 Decoder: 64000 bytes k = 128 m = 128 : 78279.5 usec, 104.651 MBps

Longhair Encoded k=128 data blocks with m=128 recovery blocks in 81641.2 usec : 100.342 MB/s
Longhair Decoded 128 erasures in 85000.7 usec : 96.3757 MB/s

WH256 wirehair_encode(N = 128) in 12381.3 usec, 661.644 MB/s after 127.385 avg losses
WH256 wirehair_decode(N = 128) average overhead = 0.025 blocks, average reconstruct time = 9868.65 usec, 830.103 MB/s

FEC-AL Encoder(8.192 MB in 128 pieces, 128 losses): Input=518.545 MB/s, Output=518.545 MB/s, (Encode create: 3762.73 MB/s)
FEC-AL Decoder(8.192 MB in 128 pieces, 128 losses): Input=121.093 MB/s, Output=121.093 MB/s, (Overhead = 0 pieces)
```

For 128 data pieces of input and 128 data pieces of redundancy:

+ Fastest to encode: Leopard (1.26 GB/s)
+ Distant second-place: WH256 (660 MB/s), FEC-AL (515 MB/s)
+ Slowest encoders: Longhair, CM256

+ Fastest to decode: WH256 (830 MB/s)
+ Distant second-place: Leopard (480 MB/s)
+ Slowest decoders: FEC-AL, CM256, Longhair

There are a lot of variables that affect when each of these libraries should be used.
Each one is ideal in a different situation, and no one library can be called the best overall.
The situation tested mainly helps explore the trade-offs of WH256, FEC-AL and Leopard for code rate 1/2.


##### CM256: Traditional O(N^2) Cauchy matrix MDS Reed-Solomon codec

Runs at about 100 MB/s encode and decode for this case.
This is an MDS code that uses a Cauchy matrix for structure.
Other examples of this type would be most MDS Reed-Solomon codecs online: Jerasure, Zfec, ISA-L, etc.
It requires SSSE3 or newer Intel instruction sets for this speed.  Otherwise it runs much slower.
This type of software gets slower as O(K*M) where K = input count and M = recovery count.
It is practical for either small data or small recovery set up to 255 pieces.

It is available for production use under BSD license here:
http://github.com/catid/cm256
(Note that the inner loops can be optimized more by applying the GF256 library.)


##### Longhair: Binary O(N^2) Cauchy matrix MDS Reed-Solomon codec

Runs at about 100 MB/s encode and decode for this case.
This is an MDS code that uses a Cauchy matrix for structure.
This one only requires XOR operations so it can run fast on low-end processors.
Requires data is a multiple of 8 bytes.
This type of software gets slower as O(K*M) where K = input count and M = recovery count.
It is practical for either small data or small recovery set up to 255 pieces.
There is no other optimized software available online for this type of error correction code.  There is a slow version available in the Jerasure software library.

It is available for production use under BSD license here:
http://github.com/catid/longhair
(Note that the inner loops can be optimized more by applying the GF256 library.)


##### Wirehair: O(N) Hybrid LDPC Erasure Code

Encodes at 660 MB/s, and decodes at 830 MB/s for ALL cases.
This is not an MDS code.  It has about a 3% chance of failing to recover and requiring one extra block of data.
It uses mostly XOR so it only gets a little slower on lower-end processors.
This type of software gets slower as O(K) where K = input count.
This library incorporates some novel ideas that are unpublished.  The new ideas are described in the source code.
It is practical for data up to 64,000 pieces and can be used as a "fountain" code.
There is no other optimized software available online for this type of error correction code.  I believe there are some public (slow) implementations of Raptor codes available online for study.

It is available for production use under BSD license here:
http://github.com/catid/wirehair

There's a pre-production version that needs more work here using GF256 for more speed,
which is what I used for the benchmark:
http://github.com/catid/wh256


##### FEC-AL *new*: O(N^2/8) XOR Structured Convolutional Matrix Code

Encodes at 510 MB/s.  Decodes at 121 MB/s.
This is not an MDS code.  It has about a 1% chance of failing to recover and requiring one extra block of data.
This library incorporates some novel ideas that are unpublished.  The new ideas are described in the README.
It uses mostly XOR operations so only gets about 2-4x slower on lower-end processors.
It gets slower as O(K*M/8) for larger data, bounded by the speed of XOR.
This new approach is ideal for streaming erasure codes; two implementations are offered one for files and another for real-time streaming reliable data.
It is practical for data up to about 4,000 pieces and can be used as a "fountain" code.
There is no other software available online for this type of error correction code.

It is available for production use under BSD license here:
http://github.com/catid/fecal

It can also be used as a convolutional streaming code here for e.g. rUDP:
http://github.com/catid/siamese


##### Leopard-RS *new*: O(K Log M) FFT MDS Reed-Solomon codec

Encodes at 1.2 GB/s, and decodes at 480 MB/s for this case.
12x faster than existing MDS approaches to encode, and almost 5x faster to decode.
This uses a recent result from 2014 introducing a novel polynomial basis permitting FFT over fast Galois fields.
This is an MDS Reed-Solomon similar to Jerasure, Zfec, ISA-L, etc, but much faster.
It requires SSSE3 or newer Intel instruction sets for this speed.  Otherwise it runs much slower.
Requires data is a multiple of 64 bytes.
This type of software gets slower as O(K Log M) where K = input count, M = recovery count.
It is practical for extremely large data.
There is no other software available online for this type of error correction code.


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


#### Future directions:

Note that a faster decoder is described in {3} that is O(K Log M) instead,
which should be 2x faster than the current one.  However I do not fully
understand how to implement it for this field and could use some help.


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
    "Novel Polynomial Basis with Fast Fourier Transform and Its Application to Reed-Solomon Erasure Codes"
    IEEE Trans. on Information Theory, pp. 6284-6299, November, 2016.
~~~

~~~
    {2} D. G. Cantor, "On arithmetical algorithms over finite fields",
    Journal of Combinatorial Theory, Series A, vol. 50, no. 2, pp. 285-300, 1989.
~~~

~~~
    {3} Sian-Jheng Lin, Wei-Ho Chung, “An Efficient (n, k) Information
    Dispersal Algorithm for High Code Rate System over Fermat Fields,”
    IEEE Commun. Lett., vol.16, no.12, pp. 2036-2039, Dec. 2012.
~~~
	
Some papers are mirrored in the /docs/ folder.


#### Credits

Inspired by discussion with:

+ Sian-Jhen Lin <sjhenglin@gmail.com> : Author of {1} {3}, basis for Leopard
+ Bulat Ziganshin <bulat.ziganshin@gmail.com> : Author of FastECC
+ Yutaka Sawada <tenfon@outlook.jp> : Author of MultiPar

This software was written entirely by myself ( Christopher A. Taylor mrcatid@gmail.com ). If you find it useful and would like to buy me a coffee, consider tipping.

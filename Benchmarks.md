# Benchmarks:

On my (few year old) laptop using AVX2 instruction set.  I'm not being very rigorous here.  The point is that it's really fast.

Some example performance measurements:

```
Leopard Encoder(0.256 MB in 100 pieces, 10 losses): Input=5333.33 MB/s, Output=533.333 MB/s
Leopard Decoder(0.256 MB in 100 pieces, 10 losses): Input=1695.36 MB/s, Output=169.536 MB/s

Leopard Encoder(0.256 MB in 100 pieces, 20 losses): Input=3878.79 MB/s, Output=775.758 MB/s
Leopard Decoder(0.256 MB in 100 pieces, 20 losses): Input=833.876 MB/s, Output=166.775 MB/s

Leopard Encoder(8.192 MB in 128 pieces, 128 losses): Input=1964.98 MB/s, Output=1964.98 MB/s
Leopard Decoder(8.192 MB in 128 pieces, 128 losses): Input=600.542 MB/s, Output=600.542 MB/s

Leopard Encoder(2.56 MB in 1000 pieces, 200 losses): Input=1942.34 MB/s, Output=388.467 MB/s
Leopard Decoder(2.56 MB in 1000 pieces, 200 losses): Input=367.109 MB/s, Output=73.4219 MB/s

Leopard Encoder(2.56 MB in 1000 pieces, 1000 losses): Input=1038.54 MB/s, Output=1038.54 MB/s
Leopard Decoder(2.56 MB in 1000 pieces, 1000 losses): Input=365.876 MB/s, Output=365.876 MB/s

Leopard Encoder(83.8861 MB in 32768 pieces, 32768 losses): Input=471.209 MB/s, Output=471.209 MB/s
Leopard Decoder(83.8861 MB in 32768 pieces, 32768 losses): Input=164.957 MB/s, Output=164.957 MB/s

Leopard Encoder(83.8861 MB in 32768 pieces, 2048 losses): Input=1359.71 MB/s, Output=84.982 MB/s
Leopard Decoder(83.8861 MB in 32768 pieces, 2048 losses): Input=169.359 MB/s, Output=10.585 MB/s
```

The benchmark project in the software can be configured to test any situation you are interested in.


#### Comparisons:

Comparing performance from all my error correction code libraries, on my laptop:

To summarize, a set of 128 of 64 KB data files are supplemented by about 128 redundant code pieces (encoded) meaning a code rate of 1/2.  From those redundant code pieces the original set is recovered (decoded).

The results are all from libraries I've written over the past few years.  They all have the same vector-optimized inner loops, but the types of error correction codes are different.

```
For 64KB data chunks:

Leopard Encoder(8.192 MB in 128 pieces, 128 losses): Input=1964.98 MB/s, Output=1964.98 MB/s
Leopard Decoder(8.192 MB in 128 pieces, 128 losses): Input=600.542 MB/s, Output=600.542 MB/s

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

+ Fastest to encode: Leopard (1.96 GB/s)
+ Distant second-place: WH256 (660 MB/s), FEC-AL (515 MB/s)
+ Slowest encoders: Longhair, CM256

+ Fastest to decode: WH256 (830 MB/s)
+ Distant second-place: Leopard (600 MB/s)
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

Encodes at 2 GB/s, and decodes at 600 MB/s for this case.
20x faster than existing MDS approaches to encode, and 6x faster to decode.
This uses a recent result from 2014 introducing a novel polynomial basis permitting FFT over fast Galois fields.
This is an MDS Reed-Solomon similar to Jerasure, Zfec, ISA-L, etc, but much faster.
It requires SSSE3 or newer Intel instruction sets for this speed.  Otherwise it runs much slower.
Requires data is a multiple of 64 bytes.
This type of software gets slower as O(K Log M) where K = input count, M = recovery count.
It is practical for extremely large data.
There is no other software available online for this type of error correction code.



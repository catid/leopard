# Leopard-RS
## Leopard Reed-Solomon Error Correction Codes in C

Leopard-RS is a portable, fast library for Forward Error Correction.
From a block of equally sized original data pieces, it generates recovery
symbols that can be used to recover lost original data.

* It requires that data pieces are all a fixed size, a multiple of 64 bytes.
* The original and recovery data must not exceed 65536 pieces.


#### Motivation:

It gets slower as O(N Log N) in the input data size, and its inner loops are
vectorized using the best approaches available on modern processors, using the
fastest finite fields (8-bit or 16-bit Galois fields) for bulk data.

It sets new speed records for MDS encoding and decoding of large data.
It is also the only open-source production ready software for this purpose
available today.

Example applications are data recovery software and data center replication.


#### Encoder API:

```
#include "leopard.h"
```

For full documentation please read `leopard.h`.

+ `leo_init()` : Initialize library.
+ `leo_encode_work_count()` : Calculate the number of work_data buffers to provide to leo_encode().
+ `leo_encode()`: Generate recovery data.


#### Decoder API:

```
#include "leopard.h"
```

For full documentation please read `leopard.h`.

+ `leo_init()` : Initialize library.
+ `leo_decode_work_count()` : Calculate the number of work_data buffers to provide to leo_decode().
+ `leo_decode()` : Generate recovery data.


#### Benchmarks:

```
TODO
```


#### Comparisons:

```
TODO
```


#### Background

This library implements an MDS erasure code introduced in this paper:

~~~
    S.-J. Lin,  T. Y. Al-Naffouri, Y. S. Han, and W.-H. Chung,
    "Novel Polynomial Basis with Fast Fourier Transform and Its Application to Reed-Solomon Erasure Codes"
    IEEE Trans. on Information Theory, pp. 6284-6299, November, 2016.
~~~

The paper is available here: [http://ct.ee.ntust.edu.tw/it2016-2.pdf](http://ct.ee.ntust.edu.tw/it2016-2.pdf)
And also mirrored in the /docs/ folder.

The high-level summary is that instead of using complicated fields,
an additive FFT was introduced that works with familiar Galois fields for the first time.
This is actually a huge new result that will change how Reed-Solomon codecs will be written.

My contribution is extending the ALTMAP approach from Jerasure
for 16-bit Galois fields out to 64 bytes to enable AVX2 speedups,
and marry it with the row parallelism introduced by ISA-L.


#### Credits

The idea is the brain-child of S.-J. Lin.  He is a super bright guy who should be recognized more widely!

This software was written entirely by myself ( Christopher A. Taylor mrcatid@gmail.com ). If you find it useful and would like to buy me a coffee, consider tipping.

# Lin-Han-Chung RS Codes
This is an attempt at implementing a fast version of the algorithm described here:

~~~
    S.-J. Lin,  T. Y. Al-Naffouri, Y. S. Han, and W.-H. Chung,
    "Novel Polynomial Basis with Fast Fourier Transform and Its Application to Reed-Solomon Erasure Codes"
    IEEE Trans. on Information Theory, pp. 6284-6299, November, 2016.
~~~
Available here: [http://ct.ee.ntust.edu.tw/it2016-2.pdf](http://ct.ee.ntust.edu.tw/it2016-2.pdf)

I wasn't able to figure out how to speed up the critical piece of code:

~~~
  // z = x + y (mod Q)
  static inline GFSymbol AddModQ(GFSymbol a, GFSymbol b)
  {
      const unsigned sum = (unsigned)a + b;

      // Partial reduction step, allowing for Q to be returned
      return static_cast<GFSymbol>(sum + (sum >> kGFBits));
  }

  // vx[] += vy[] * z
  static void muladd_mem(GFSymbol * vx, const GFSymbol * vy, GFSymbol z, unsigned symbolCount)
  {
      for (unsigned i = 0; i < symbolCount; ++i)
      {
          const GFSymbol a = vy[i];
          if (a == 0)
              continue;

          const GFSymbol sum = static_cast<GFSymbol>(AddModQ(GFLog[a], z));

          vx[i] ^= GFExp[sum];
      }
  }
~~~

This routine basically takes all the runtime.  Since it involves normal addition via `AddModQ()`, it is pretty slow and is hard to vectorize.

Math over GF(2^^16) is a lot faster because multiplies can be decomposed linearly into separate operations that can be done via vector instructions.

So unfortunately it looks like while this field does have an FFT, the field does not permit fast math.

Maybe someone else can figure out how to speed this up?

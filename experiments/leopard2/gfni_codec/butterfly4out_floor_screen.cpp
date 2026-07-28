// AVX2FF8Butterfly4Out isolation bench.
//
// That kernel owns ~40% of legacy-high GF8 encode and runs at a flat
// 116-188 GB/s of combined read+write from 1 KiB to 64 KiB per shard, which is
// well under its port bound and does not look bandwidth limited.  This harness
// reproduces its exact shape and compares variants, plus two floors:
//   memfloor  - four loads, four stores, no arithmetic at all
//   xorfloor  - four loads, four stores, one xor each (pure XOR cost)
// If the real kernel is close to memfloor there is nothing to win.
#include <immintrin.h>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <chrono>
#include <vector>
#include <algorithm>

namespace {
typedef std::chrono::steady_clock Clock;

uint64_t Rng(uint64_t& s){ s^=s<<13; s^=s>>7; s^=s<<17; return s; }

struct Map { uint8_t low[16]; uint8_t high[16]; uint8_t affine[32]; };

Map MakeMap(uint64_t& seed)
{
    Map m; uint8_t col[8];
    for (unsigned i=0;i<8;++i) col[i]=(uint8_t)Rng(seed);
    for (unsigned n=0;n<16;++n){ uint8_t lo=0,hi=0;
        for (unsigned b=0;b<4;++b) if(n&(1u<<b)){lo^=col[b];hi^=col[b+4];}
        m.low[n]=lo; m.high[n]=hi; }
    uint64_t mat=0;
    for (unsigned o=0;o<8;++o) for (unsigned i=0;i<8;++i)
        if (col[i]&(1u<<o)) mat |= 1ULL<<(8*(7-o)+i);
    for (unsigned r=0;r<4;++r) memcpy(m.affine+r*8,&mat,8);
    return m;
}

#define LOADS(P) \
    __m256i x0=_mm256_loadu_si256((const __m256i*)(i0+P)); \
    __m256i x1=_mm256_loadu_si256((const __m256i*)(i1+P)); \
    __m256i x2=_mm256_loadu_si256((const __m256i*)(i2+P)); \
    __m256i x3=_mm256_loadu_si256((const __m256i*)(i3+P));
#define STORES(P) \
    _mm256_storeu_si256((__m256i*)(o0+P),x0); \
    _mm256_storeu_si256((__m256i*)(o1+P),x1); \
    _mm256_storeu_si256((__m256i*)(o2+P),x2); \
    _mm256_storeu_si256((__m256i*)(o3+P),x3);
#define BODY(M01,M23,M02) \
    x1=_mm256_xor_si256(x1,x0); \
    x0=_mm256_xor_si256(x0,_mm256_gf2p8affine_epi64_epi8(x1,M01,0)); \
    x3=_mm256_xor_si256(x3,x2); \
    x2=_mm256_xor_si256(x2,_mm256_gf2p8affine_epi64_epi8(x3,M23,0)); \
    x2=_mm256_xor_si256(x2,x0); \
    x3=_mm256_xor_si256(x3,x1); \
    x0=_mm256_xor_si256(x0,_mm256_gf2p8affine_epi64_epi8(x2,M02,0)); \
    x1=_mm256_xor_si256(x1,_mm256_gf2p8affine_epi64_epi8(x3,M02,0));

// (0) pure load/store floor
__attribute__((target("avx2")))
void MemFloor(const uint8_t* i0,const uint8_t* i1,const uint8_t* i2,const uint8_t* i3,
              uint8_t* o0,uint8_t* o1,uint8_t* o2,uint8_t* o3,
              const Map&,const Map&,const Map&,size_t n)
{ for(size_t p=0;p<n;p+=32){ LOADS(p) STORES(p) } }

// (1) load/store plus the eight XORs the algorithm needs, no multiplies
__attribute__((target("avx2")))
void XorFloor(const uint8_t* i0,const uint8_t* i1,const uint8_t* i2,const uint8_t* i3,
              uint8_t* o0,uint8_t* o1,uint8_t* o2,uint8_t* o3,
              const Map&,const Map&,const Map&,size_t n)
{
    for(size_t p=0;p<n;p+=32){ LOADS(p)
        x1=_mm256_xor_si256(x1,x0); x3=_mm256_xor_si256(x3,x2);
        x2=_mm256_xor_si256(x2,x0); x3=_mm256_xor_si256(x3,x1);
        x0=_mm256_xor_si256(x0,x2); x1=_mm256_xor_si256(x1,x3);
        x2=_mm256_xor_si256(x2,x1); x3=_mm256_xor_si256(x3,x0);
        STORES(p) }
}

// (2) production shape: advancing pointers, 32 bytes per iteration
__attribute__((target("avx2,gfni")))
void Production(const uint8_t* i0,const uint8_t* i1,const uint8_t* i2,const uint8_t* i3,
                uint8_t* o0,uint8_t* o1,uint8_t* o2,uint8_t* o3,
                const Map& m01,const Map& m23,const Map& m02,size_t n)
{
    const __m256i M01=_mm256_loadu_si256((const __m256i*)m01.affine);
    const __m256i M23=_mm256_loadu_si256((const __m256i*)m23.affine);
    const __m256i M02=_mm256_loadu_si256((const __m256i*)m02.affine);
    while(n>=32){
        __m256i x0=_mm256_loadu_si256((const __m256i*)i0);
        __m256i x1=_mm256_loadu_si256((const __m256i*)i1);
        __m256i x2=_mm256_loadu_si256((const __m256i*)i2);
        __m256i x3=_mm256_loadu_si256((const __m256i*)i3);
        BODY(M01,M23,M02)
        _mm256_storeu_si256((__m256i*)o0,x0);
        _mm256_storeu_si256((__m256i*)o1,x1);
        _mm256_storeu_si256((__m256i*)o2,x2);
        _mm256_storeu_si256((__m256i*)o3,x3);
        i0+=32;i1+=32;i2+=32;i3+=32;o0+=32;o1+=32;o2+=32;o3+=32;n-=32;
    }
}

// (3) single offset, one tile
__attribute__((target("avx2,gfni")))
void Offset1(const uint8_t* i0,const uint8_t* i1,const uint8_t* i2,const uint8_t* i3,
             uint8_t* o0,uint8_t* o1,uint8_t* o2,uint8_t* o3,
             const Map& m01,const Map& m23,const Map& m02,size_t n)
{
    const __m256i M01=_mm256_loadu_si256((const __m256i*)m01.affine);
    const __m256i M23=_mm256_loadu_si256((const __m256i*)m23.affine);
    const __m256i M02=_mm256_loadu_si256((const __m256i*)m02.affine);
    for(size_t p=0;p<n;p+=32){ LOADS(p) BODY(M01,M23,M02) STORES(p) }
}

// (4) single offset, two tiles per iteration
__attribute__((target("avx2,gfni")))
void Offset2(const uint8_t* i0,const uint8_t* i1,const uint8_t* i2,const uint8_t* i3,
             uint8_t* o0,uint8_t* o1,uint8_t* o2,uint8_t* o3,
             const Map& m01,const Map& m23,const Map& m02,size_t n)
{
    const __m256i M01=_mm256_loadu_si256((const __m256i*)m01.affine);
    const __m256i M23=_mm256_loadu_si256((const __m256i*)m23.affine);
    const __m256i M02=_mm256_loadu_si256((const __m256i*)m02.affine);
    size_t p=0;
    for(;p+64<=n;p+=64){
        { LOADS(p) BODY(M01,M23,M02) STORES(p) }
        { LOADS(p+32) BODY(M01,M23,M02) STORES(p+32) }
    }
    for(;p<n;p+=32){ LOADS(p) BODY(M01,M23,M02) STORES(p) }
}

// (5) non-temporal stores (outputs are not re-read by this pass)
__attribute__((target("avx2,gfni")))
void NonTemporal(const uint8_t* i0,const uint8_t* i1,const uint8_t* i2,const uint8_t* i3,
                 uint8_t* o0,uint8_t* o1,uint8_t* o2,uint8_t* o3,
                 const Map& m01,const Map& m23,const Map& m02,size_t n)
{
    const __m256i M01=_mm256_loadu_si256((const __m256i*)m01.affine);
    const __m256i M23=_mm256_loadu_si256((const __m256i*)m23.affine);
    const __m256i M02=_mm256_loadu_si256((const __m256i*)m02.affine);
    for(size_t p=0;p<n;p+=32){ LOADS(p) BODY(M01,M23,M02)
        _mm256_stream_si256((__m256i*)(o0+p),x0);
        _mm256_stream_si256((__m256i*)(o1+p),x1);
        _mm256_stream_si256((__m256i*)(o2+p),x2);
        _mm256_stream_si256((__m256i*)(o3+p),x3); }
    _mm_sfence();
}


// ---------------------------------------------------------------- radix eight
//
// Three transform layers per load/store round instead of two.  Inverse DIT over
// eight coordinates:
//   layer 1 (dist 1): (0,1) (2,3) (4,5) (6,7)
//   layer 2 (dist 2): (0,2) (1,3) (4,6) (5,7)
//   layer 3 (dist 4): (0,4) (1,5) (2,6) (3,7)
// with inverse butterfly(x,y,m): y ^= x; x ^= mul(y,m).
//
// Eight values are live across layer three, but out-of-place means the inputs
// are consumed into registers up front and the outputs written at the end, so
// only eight data vectors are live, not sixteen.  All seven multiplier
// matrices are memory operands: the kernel is memory bound, so the extra load
// uops are free and the register file stays clear.
#define IBFLY(X,Y,M) \
    Y = _mm256_xor_si256(Y, X); \
    X = _mm256_xor_si256(X, _mm256_gf2p8affine_epi64_epi8( \
            Y, _mm256_loadu_si256((const __m256i*)(M)), 0));

__attribute__((target("avx2,gfni")))
void Radix8Out(const uint8_t* const* in, uint8_t* const* out,
               const Map* const* m, size_t n)
{
    for (size_t p = 0; p < n; p += 32)
    {
        __m256i x0=_mm256_loadu_si256((const __m256i*)(in[0]+p));
        __m256i x1=_mm256_loadu_si256((const __m256i*)(in[1]+p));
        __m256i x2=_mm256_loadu_si256((const __m256i*)(in[2]+p));
        __m256i x3=_mm256_loadu_si256((const __m256i*)(in[3]+p));
        __m256i x4=_mm256_loadu_si256((const __m256i*)(in[4]+p));
        __m256i x5=_mm256_loadu_si256((const __m256i*)(in[5]+p));
        __m256i x6=_mm256_loadu_si256((const __m256i*)(in[6]+p));
        __m256i x7=_mm256_loadu_si256((const __m256i*)(in[7]+p));
        IBFLY(x0,x1,m[0]->affine) IBFLY(x2,x3,m[1]->affine)
        IBFLY(x4,x5,m[2]->affine) IBFLY(x6,x7,m[3]->affine)
        IBFLY(x0,x2,m[4]->affine) IBFLY(x1,x3,m[4]->affine)
        IBFLY(x4,x6,m[5]->affine) IBFLY(x5,x7,m[5]->affine)
        IBFLY(x0,x4,m[6]->affine) IBFLY(x1,x5,m[6]->affine)
        IBFLY(x2,x6,m[6]->affine) IBFLY(x3,x7,m[6]->affine)
        _mm256_storeu_si256((__m256i*)(out[0]+p),x0);
        _mm256_storeu_si256((__m256i*)(out[1]+p),x1);
        _mm256_storeu_si256((__m256i*)(out[2]+p),x2);
        _mm256_storeu_si256((__m256i*)(out[3]+p),x3);
        _mm256_storeu_si256((__m256i*)(out[4]+p),x4);
        _mm256_storeu_si256((__m256i*)(out[5]+p),x5);
        _mm256_storeu_si256((__m256i*)(out[6]+p),x6);
        _mm256_storeu_si256((__m256i*)(out[7]+p),x7);
    }
}

// In-place radix-two layer, the extra pass radix-four needs to reach three
// layers.  Same shape as AVX2FF8Butterfly2.
__attribute__((target("avx2,gfni")))
void Radix2InPlace(uint8_t* x, uint8_t* y, const Map& m, size_t n)
{
    const __m256i M=_mm256_loadu_si256((const __m256i*)m.affine);
    for (size_t p=0;p<n;p+=32){
        __m256i a=_mm256_loadu_si256((const __m256i*)(x+p));
        __m256i b=_mm256_loadu_si256((const __m256i*)(y+p));
        b=_mm256_xor_si256(b,a);
        a=_mm256_xor_si256(a,_mm256_gf2p8affine_epi64_epi8(b,M,0));
        _mm256_storeu_si256((__m256i*)(x+p),a);
        _mm256_storeu_si256((__m256i*)(y+p),b);
    }
}

typedef void (*Kern)(const uint8_t*,const uint8_t*,const uint8_t*,const uint8_t*,
                     uint8_t*,uint8_t*,uint8_t*,uint8_t*,
                     const Map&,const Map&,const Map&,size_t);

struct Buf { std::vector<uint8_t*> p; size_t bytes;
    ~Buf(){ for(size_t i=0;i<p.size();++i) free(p[i]); } };

void Alloc(Buf& b,size_t count,size_t bytes,uint64_t& seed)
{
    b.bytes=bytes; b.p.resize(count);
    for(size_t i=0;i<count;++i){ void* q=NULL;
        if(posix_memalign(&q,64,bytes)!=0) abort();
        uint8_t* d=(uint8_t*)q;
        for(size_t j=0;j<bytes;++j) d[j]=(uint8_t)Rng(seed);
        b.p[i]=d; }
}

// One legacy-high encode block: 8 radix-4 groups reading 32 source shards
// out-of-place into 32 work shards, repeated for 7 message blocks.
double Sweep(Kern k, Buf& src, Buf& work, const Map& a,const Map& b,const Map& c,
             unsigned reps)
{
    double best=1e30;
    for(unsigned r=0;r<reps;++r){
        const Clock::time_point t0=Clock::now();
        for(unsigned blk=0; blk<7; ++blk)
            for(unsigned g=0; g<8; ++g){
                const unsigned s=blk*32+g*4;
                k(src.p[s],src.p[s+1],src.p[s+2],src.p[s+3],
                  work.p[g*4],work.p[g*4+1],work.p[g*4+2],work.p[g*4+3],
                  a,b,c,src.bytes);
            }
        best=std::min(best,std::chrono::duration<double,std::micro>(
            Clock::now()-t0).count());
    }
    return best;
}


// In-place radix-eight: same eight-live-vector profile as the out-of-place
// form, because the transform happens entirely in registers between the eight
// loads and the eight stores.  Distance-strided, matching the range kernels.
__attribute__((target("avx2,gfni")))
void Radix8InPlaceRange(uint8_t* const* work, unsigned distance,
                        const Map* const* m, size_t n)
{
    for (unsigned i = 0; i < distance; ++i)
    {
        uint8_t* v0 = work[i];
        uint8_t* v1 = work[i + distance];
        uint8_t* v2 = work[i + distance * 2];
        uint8_t* v3 = work[i + distance * 3];
        uint8_t* v4 = work[i + distance * 4];
        uint8_t* v5 = work[i + distance * 5];
        uint8_t* v6 = work[i + distance * 6];
        uint8_t* v7 = work[i + distance * 7];
        for (size_t p = 0; p < n; p += 32)
        {
            __m256i x0=_mm256_loadu_si256((const __m256i*)(v0+p));
            __m256i x1=_mm256_loadu_si256((const __m256i*)(v1+p));
            __m256i x2=_mm256_loadu_si256((const __m256i*)(v2+p));
            __m256i x3=_mm256_loadu_si256((const __m256i*)(v3+p));
            __m256i x4=_mm256_loadu_si256((const __m256i*)(v4+p));
            __m256i x5=_mm256_loadu_si256((const __m256i*)(v5+p));
            __m256i x6=_mm256_loadu_si256((const __m256i*)(v6+p));
            __m256i x7=_mm256_loadu_si256((const __m256i*)(v7+p));
            IBFLY(x0,x1,m[0]->affine) IBFLY(x2,x3,m[1]->affine)
            IBFLY(x4,x5,m[2]->affine) IBFLY(x6,x7,m[3]->affine)
            IBFLY(x0,x2,m[4]->affine) IBFLY(x1,x3,m[4]->affine)
            IBFLY(x4,x6,m[5]->affine) IBFLY(x5,x7,m[5]->affine)
            IBFLY(x0,x4,m[6]->affine) IBFLY(x1,x5,m[6]->affine)
            IBFLY(x2,x6,m[6]->affine) IBFLY(x3,x7,m[6]->affine)
            _mm256_storeu_si256((__m256i*)(v0+p),x0);
            _mm256_storeu_si256((__m256i*)(v1+p),x1);
            _mm256_storeu_si256((__m256i*)(v2+p),x2);
            _mm256_storeu_si256((__m256i*)(v3+p),x3);
            _mm256_storeu_si256((__m256i*)(v4+p),x4);
            _mm256_storeu_si256((__m256i*)(v5+p),x5);
            _mm256_storeu_si256((__m256i*)(v6+p),x6);
            _mm256_storeu_si256((__m256i*)(v7+p),x7);
        }
    }
}

// In-place radix-four range, the shape production uses.
__attribute__((target("avx2,gfni")))
void Radix4InPlaceRange(uint8_t* const* work, unsigned distance,
                        const Map& m01, const Map& m23, const Map& m02, size_t n)
{
    const __m256i M01=_mm256_loadu_si256((const __m256i*)m01.affine);
    const __m256i M23=_mm256_loadu_si256((const __m256i*)m23.affine);
    const __m256i M02=_mm256_loadu_si256((const __m256i*)m02.affine);
    for (unsigned i = 0; i < distance; ++i)
    {
        uint8_t* i0=work[i]; uint8_t* i1=work[i+distance];
        uint8_t* i2=work[i+distance*2]; uint8_t* i3=work[i+distance*3];
        uint8_t* o0=i0; uint8_t* o1=i1; uint8_t* o2=i2; uint8_t* o3=i3;
        for (size_t p=0;p<n;p+=32){ LOADS(p) BODY(M01,M23,M02) STORES(p) }
    }
}

// Three in-place layers over 32 shards: radix-four then radix-two, versus one
// radix-eight round.
double SweepInPlace3(Buf& work, const Map& a,const Map& b,const Map& c,
                     const Map* mm[7], bool use8, unsigned reps)
{
    double best=1e30;
    for(unsigned r=0;r<reps;++r){
        const Clock::time_point t0=Clock::now();
        for(unsigned blk=0; blk<7; ++blk){
            if (use8) {
                Radix8InPlaceRange(work.p.data(), 4, mm, work.bytes);
            } else {
                Radix4InPlaceRange(work.p.data(), 8, a, b, c, work.bytes);
                for(unsigned i=0;i<16;++i)
                    Radix2InPlace(work.p[i], work.p[i+16], c, work.bytes);
            }
        }
        best=std::min(best,std::chrono::duration<double,std::micro>(
            Clock::now()-t0).count());
    }
    return best;
}

// Three layers over 32 shards, current shape: radix-four out-of-place then one
// in-place radix-two sweep.
double SweepThreeLayersRadix4(Buf& src, Buf& work, const Map& a,const Map& b,
                              const Map& c, unsigned reps)
{
    double best=1e30;
    for(unsigned r=0;r<reps;++r){
        const Clock::time_point t0=Clock::now();
        for(unsigned blk=0; blk<7; ++blk){
            for(unsigned g=0; g<8; ++g){
                const unsigned s=blk*32+g*4;
                Production(src.p[s],src.p[s+1],src.p[s+2],src.p[s+3],
                           work.p[g*4],work.p[g*4+1],work.p[g*4+2],work.p[g*4+3],
                           a,b,c,src.bytes);
            }
            for(unsigned i=0;i<16;++i)
                Radix2InPlace(work.p[i], work.p[i+16], c, work.bytes);
        }
        best=std::min(best,std::chrono::duration<double,std::micro>(
            Clock::now()-t0).count());
    }
    return best;
}

// Three layers over 32 shards in one out-of-place radix-eight round.
double SweepThreeLayersRadix8(Buf& src, Buf& work, const Map* mm[7], unsigned reps)
{
    double best=1e30;
    for(unsigned r=0;r<reps;++r){
        const Clock::time_point t0=Clock::now();
        for(unsigned blk=0; blk<7; ++blk)
            for(unsigned g=0; g<4; ++g){
                const uint8_t* in[8]; uint8_t* out[8];
                for(unsigned j=0;j<8;++j){
                    in[j]=src.p[blk*32+g*8+j];
                    out[j]=work.p[g*8+j];
                }
                Radix8Out(in,out,mm,src.bytes);
            }
        best=std::min(best,std::chrono::duration<double,std::micro>(
            Clock::now()-t0).count());
    }
    return best;
}

} // namespace

int main(int argc,char** argv)
{
    const unsigned reps=argc>1?(unsigned)atoi(argv[1]):30;
    uint64_t seed=0x9E3779B97F4A7C15ULL;
    const Map a=MakeMap(seed),b=MakeMap(seed),c=MakeMap(seed);
    const char* names[]={"memfloor","xorfloor","production","offset1","offset2","nontemporal"};
    Kern kerns[]={MemFloor,XorFloor,Production,Offset1,Offset2,NonTemporal};
    const size_t sizes[]={1024,4096,16384,65536};
    printf("K=224 R=32 shape: 7 blocks x 8 radix-4 groups, 224 shards in / 32 out\n");
    printf("%9s %11s %11s %11s %11s %11s %11s\n","bytes",
           names[0],names[1],names[2],names[3],names[4],names[5]);
    for(size_t si=0;si<sizeof(sizes)/sizeof(sizes[0]);++si){
        Buf src,work; Alloc(src,224,sizes[si],seed); Alloc(work,32,sizes[si],seed);
        double t[6];
        for(int ki=0;ki<6;++ki) t[ki]=Sweep(kerns[ki],src,work,a,b,c,reps);
        printf("%9zu",sizes[si]);
        for(int ki=0;ki<6;++ki) printf(" %11.2f",t[ki]);
        printf("   prod/mem=%.2fx  best=%s\n", t[2]/t[0],
               names[(int)(std::min_element(t+2,t+6)-t)]);
        const double gb=7.0*32.0*sizes[si]*2/1e9;
        printf("%9s"," GB/s:");
        for(int ki=0;ki<6;++ki) printf(" %11.1f",gb/(t[ki]*1e-6));
        printf("\n");
    }

    printf("\nThree transform layers over 32 shards per block, 7 blocks\n");
    printf("%9s %14s %14s %9s\n","bytes","radix4+radix2","radix8-out","gain");
    const Map* mm[7]={&a,&b,&c,&a,&b,&c,&a};
    for(size_t si=0;si<sizeof(sizes)/sizeof(sizes[0]);++si){
        Buf src,work; Alloc(src,224,sizes[si],seed); Alloc(work,32,sizes[si],seed);
        const double t4=SweepThreeLayersRadix4(src,work,a,b,c,reps);
        const double t8=SweepThreeLayersRadix8(src,work,mm,reps);
        printf("%9zu %14.2f %14.2f %8.3fx\n",sizes[si],t4,t8,t4/t8);
    }
    printf("\nThree IN-PLACE layers over 32 shards, 7 rounds\n");
    printf("%9s %16s %14s %9s\n","bytes","radix4+radix2","radix8","gain");
    for(size_t si=0;si<sizeof(sizes)/sizeof(sizes[0]);++si){
        Buf work; Alloc(work,32,sizes[si],seed);
        const double t4=SweepInPlace3(work,a,b,c,mm,false,reps);
        const double t8=SweepInPlace3(work,a,b,c,mm,true,reps);
        printf("%9zu %16.2f %14.2f %8.3fx\n",sizes[si],t4,t8,t4/t8);
    }
    return 0;
}

// Delegation/counter tests only: fake callbacks do not touch payload bytes.
#define LEO_GF16_CALLBACK_NO_MAIN 1
#include "gf16_callback_probe.cpp"

namespace {
unsigned char storage[8][80] = {};
void* p[] = {storage[0],storage[1],storage[2],storage[3],storage[4],storage[5],storage[6],storage[7]};
unsigned delegated = 0;
bool expected_hint = false;
void FixedFake(void* d,const void* s,uint16_t log,uint64_t bytes)
{ Require(d==p[0] && s==p[1] && log==5 && bytes==66,"fixed arguments"); ++delegated; }
void MemoryFake(void* d,const void* s,uint64_t bytes)
{ Require(d==p[0] && s==p[1] && bytes==66,"memory arguments"); ++delegated; }
void Memory2Fake(void* d,const void* s0,const void* s1,uint64_t bytes)
{ Require(d==p[0] && s0==p[1] && s1==p[2] && bytes==66,"xor2 arguments"); ++delegated; }
void Memory4Fake(void* d0,const void* s0,void* d1,const void* s1,
                 void* d2,const void* s2,void* d3,const void* s3,uint64_t bytes)
{ Require(d0==p[0] && s0==p[1] && d1==p[2] && s1==p[3] && d2==p[4] && s2==p[5] &&
    d3==p[6] && s3==p[7] && bytes==66,"xor4 arguments"); ++delegated; }
void PairFake(void* x,void* y,uint16_t log,uint64_t bytes)
{ FixedFake(x,y,log,bytes); }
void PairOutFake(const void* x,const void* y,void* u,void* v,uint16_t log,uint64_t bytes)
{ Require(x==p[0] && y==p[1] && u==p[2] && v==p[3] && log==65535 && bytes==66,
    "pair-out arguments"); ++delegated; }
void QuadFake(void* a,void* b,void* c,void* d,uint16_t l0,uint16_t l1,uint16_t l2,uint64_t bytes)
{ Require(a==p[0] && b==p[1] && c==p[2] && d==p[3] && l0==17 && l1==65535 && l2==19 &&
    bytes==66,"quad arguments"); ++delegated; }
void QuadOutFake(const void* a,const void* b,const void* c,const void* d,
    void* e,void* f,void* g,void* h,uint16_t l0,uint16_t l1,uint16_t l2,uint64_t bytes)
{ Require(a==p[0] && b==p[1] && c==p[2] && d==p[3] && e==p[4] && f==p[5] && g==p[6] &&
    h==p[7] && l0==7 && l1==8 && l2==9 && bytes==66,"quad-out arguments"); ++delegated; }
void RangeFake(void* const* work,unsigned distance,uint16_t l0,uint16_t l1,uint16_t l2,
    uint64_t bytes,bool hint)
{ Require(work==p && distance==4 && l0==65535 && l1==5 && l2==65535 && bytes==66 &&
    hint==expected_hint,"range arguments"); ++delegated; }
template<class F> void Reject(F action)
{
    bool rejected = false;
    try { action(); } catch (const std::runtime_error&) { rejected = true; }
    Require(rejected,"expected diagnostic rejection");
}
}

int main()
{
    using namespace callback_probe;
    omp_set_dynamic(0); omp_set_num_threads(1);
    Ops original = {};
    original.kind = LEO2_BACKEND_GFNI; original.name = "fake-gfni";
    original.ff16_multiply = original.ff16_multiply_add = FixedFake;
    original.xor_memory = original.copy_memory = MemoryFake;
    original.xor_memory_2to1 = Memory2Fake; original.xor_memory4 = Memory4Fake;
    original.ff16_ifft_butterfly2 = original.ff16_fft_butterfly2 = PairFake;
    original.ff16_fft_butterfly2_out = original.ff16_ifft_butterfly2_xor = PairOutFake;
    original.ff16_ifft_butterfly4 = original.ff16_fft_butterfly4 = QuadFake;
    original.ff16_ifft_butterfly4_out = original.ff16_fft_butterfly4_out = QuadOutFake;
    original.ff16_ifft_butterfly4_range = original.ff16_fft_butterfly4_range = RangeFake;
    // An unrelated live GF8 entry must stay byte-for-byte identical too.
    original.ff8_multiply = FixedFake;
    Ops view = View(original);
    CheckView(view, original);
    Require(view.kind == original.kind && view.name == original.name && view.ff8_multiply == FixedFake,
            "nonobserved identity");
    Reset();
    {
        Scope scope(original);
        view.ff16_multiply(p[0],p[1],5,66); view.ff16_multiply_add(p[0],p[1],5,66);
        view.xor_memory(p[0],p[1],66); view.xor_memory_2to1(p[0],p[1],p[2],66);
        view.xor_memory4(p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[7],66);
        view.copy_memory(p[0],p[1],66);
        view.ff16_ifft_butterfly2(p[0],p[1],5,66); view.ff16_fft_butterfly2(p[0],p[1],5,66);
        view.ff16_fft_butterfly2_out(p[0],p[1],p[2],p[3],65535,66);
        view.ff16_ifft_butterfly2_xor(p[0],p[1],p[2],p[3],65535,66);
        view.ff16_ifft_butterfly4(p[0],p[1],p[2],p[3],17,65535,19,66);
        view.ff16_fft_butterfly4(p[0],p[1],p[2],p[3],17,65535,19,66);
        view.ff16_ifft_butterfly4_out(p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[7],7,8,9,66);
        view.ff16_fft_butterfly4_out(p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[7],7,8,9,66);
        expected_hint = true; view.ff16_ifft_butterfly4_range(p,4,65535,5,65535,66,true);
        expected_hint = false; view.ff16_fft_butterfly4_range(p,4,65535,5,65535,66,false);
        Require(state.calls == 16 && state.bucket_count == 16 && delegated == 16,"all callbacks delegated");
        for (unsigned i = 0; i < Count; ++i)
            Require(state.buckets[i].operation == i && state.buckets[i].calls == 1 &&
                    state.buckets[i].bytes == 66,"bucket order/count");
        Require(state.buckets[Ifft4].zero_mask == 2 && state.buckets[Fft2Out].zero_mask == 1 &&
                state.buckets[Ifft4Range].zero_mask == 5 && state.buckets[Ifft4Range].distance == 4 &&
                state.buckets[Ifft4Range].hint && !state.buckets[Fft4Range].hint,"mask/distance/hint");
        view.xor_memory(p[0],p[1],66);
        Require(state.calls == 17 && state.bucket_count == 16 && state.buckets[Xor].calls == 2,"aggregation");
        Reject([&]() { Scope nested(original); });
        Reject([]() { Reset(); });
        Require(active == &original && delegated == 17,"reentry preservation");
    }
    Reject([]() { Record(Xor,64); });
    Ops absent = {}; absent.kind = LEO2_BACKEND_GFNI;
    Ops absent_view = View(absent); CheckView(absent_view, absent);
    Require(absent_view.copy_memory == NULL && absent_view.ff16_ifft_butterfly4_range == NULL,"null callbacks");
    Reset();
    {
        Scope scope(original);
        for (unsigned i = 0; i < 256; ++i) Record(Multiply,2U*(i+1));
        Reject([]() { Record(Multiply,514); });
        Require(state.calls == 256 && state.bucket_count == 256,"bucket overflow atomicity");
        state.calls = 1000000;
        Reject([]() { Record(Multiply,2); });
        Require(state.calls == 1000000 && state.buckets[0].calls == 1,"call overflow atomicity");
    }
    Reset();
    Require(state.calls == 0 && state.bucket_count == 0 && state.pass_count == 0,"reset");
    std::puts("callback observer: 16 exact delegations, metadata/null preservation, masks/hints, aggregation, reentry and both bounds passed");
}

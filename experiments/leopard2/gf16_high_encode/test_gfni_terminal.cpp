// Narrow hook and readonly-input accumulating kernel checks for .38.5.4.14.
#include "gfni_terminal.h"
#include "leopard.h"
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <memory>
#include <stdexcept>
#include <vector>
#include <omp.h>

namespace {
void Require(bool value, const char* message)
{ if (!value) throw std::runtime_error(message); }

void Predicate()
{
    using namespace gfni_terminal;
    const Call base = {LEO2_BACKEND_GFNI,1000,200,200,256,0,32768,65536,true};
    Require(Matches(base), "target predicate");
    unsigned negatives = 0;
    unsigned Call::* fields[] = {&Call::kind,&Call::k,&Call::r,&Call::requested,&Call::side,&Call::sparse_blocks};
    for (auto field : fields)
        for (int delta : {-1,1})
        {
            Call c = base; c.*field += delta;
            Require(!Matches(c), "predicate neighbor"); ++negatives;
        }
    uint64_t Call::* sizes[] = {&Call::bytes,&Call::source_policy};
    for (auto field : sizes)
        for (int delta : {-1,1})
        {
            Call c = base; c.*field += delta;
            Require(!Matches(c), "size predicate neighbor"); ++negatives;
        }
    Call absent = base; absent.sparse_present = false;
    Require(!Matches(absent), "null descriptor is not public target"); ++negatives;
    leopard::backend::Ops ops = {}; ops.kind = LEO2_BACKEND_GFNI;
    leopard2_internal::SparseForwardPlanBatchView sparse = {};
    for (bool enabled : {false,true})
    {
        Reset(enabled);
        for (unsigned i=0; i<16; ++i)
            Require(LeoGFNITerminalExperiment(ops,32768,65536,1000,200,200,256,&sparse) == enabled,
                "hook selection");
        Require(Get().calls == 16 && Get().matches == 16 && Get().changed == (enabled ? 16U : 0U),
            "hook counts");
        bool rejected = false;
        try { LeoGFNITerminalExperiment(ops,32768,65536,1000,200,200,256,&sparse); }
        catch (const std::runtime_error&) { rejected = true; }
        Require(rejected && Get().calls == 16 && Get().matches == 16 && Get().changed == (enabled ? 16U : 0U),
            "hook overflow atomicity");
    }
    Reset(true);
    Require(!LeoGFNITerminalExperiment(ops,32768,65536,1000,200,200,256,NULL), "null hook guard");
    sparse.block_count = 1;
    Require(!LeoGFNITerminalExperiment(ops,32768,65536,1000,200,200,256,&sparse), "sparse hook guard");
    Reset(false);
    Require(!Get().enabled && Get().calls == 0 && Get().matches == 0 && Get().changed == 0, "reset");
    std::printf("terminal hook: %u negative predicates, both modes, descriptor identity and atomic 16-pass bound passed\n", negatives);
}

struct Rows {
    size_t prefix;
    std::vector<std::vector<uint8_t> > storage;
    std::vector<void*> pointers;
    Rows(unsigned count,size_t bytes,size_t offset,uint32_t seed)
        : prefix(offset ? offset : bytes ? 0 : 17),storage(count),pointers(count)
    {
        for (unsigned row=0; row<count; ++row)
        {
            storage[row] = std::vector<uint8_t>(prefix+bytes,0xa5);
            Require(storage[row].capacity() == storage[row].size(), "exact allocation end");
            for (size_t i=prefix; i<storage[row].size(); ++i)
            { seed^=seed<<13; seed^=seed>>17; seed^=seed<<5; storage[row][i]=static_cast<uint8_t>(seed); }
            pointers[row]=storage[row].data()+prefix;
        }
    }
    void Prefix() const
    { for (const auto& row : storage) for (size_t i=0; i<prefix; ++i) Require(row[i]==0xa5,"prefix canary"); }
};

void Kernels()
{
    Require(leo_init()==Leopard_Success,"initialize");
    const auto* gfni=leopard::backend::GetQualifiedOps(LEO2_BACKEND_GFNI);
    const auto* scalar=leopard::backend::GetQualifiedOps(LEO2_BACKEND_SCALAR);
    Require(gfni && scalar,"GFNI/scalar required, no skip");
    LeoGFNIFinalAccumulate(NULL,NULL,64,0,0,0,0);
    LeoGFNIFinalAccumulate(NULL,NULL,0,0,0,0,65536);
    const size_t lengths[]={0,2,62,64,66,128,130,8190,8192,16384,16386,32768,32770,65536,65598};
    unsigned cases=0;
    const auto check=[&](unsigned distance,size_t bytes,size_t offset,unsigned mask) {
        const unsigned rows=4*distance;
        Rows input(rows,bytes,offset,29+mask),work(rows,bytes,offset,29+mask);
        Rows sums(rows,bytes,offset,91+mask),expected(rows,bytes,offset,91+mask);
        const auto source_before=input.storage;
        const auto sums_before=sums.storage;
        std::vector<const void*> inputs(rows);
        for (unsigned i=0; i<rows; ++i) inputs[i]=input.pointers[i];
        const uint16_t a=(mask&1) ? 65535 : 0;
        const uint16_t b=(mask&2) ? 65535 : 32768;
        const uint16_t c=(mask&4) ? 65535 : 65534;
        const auto pair=[&](unsigned x,unsigned y,uint16_t log) {
            if (log==65535) scalar->xor_memory(work.pointers[y],work.pointers[x],bytes);
            else scalar->ff16_ifft_butterfly2(work.pointers[x],work.pointers[y],log,bytes);
        };
        for (unsigned i=0; i<distance; ++i)
        {
            pair(i,i+distance,a); pair(i+2*distance,i+3*distance,b);
            pair(i,i+2*distance,c); pair(i+distance,i+3*distance,c);
        }
        for (unsigned i=0; i<rows; ++i) scalar->xor_memory(expected.pointers[i],work.pointers[i],bytes);
        LeoGFNIFinalAccumulate(inputs.data(),sums.pointers.data(),distance,a,b,c,bytes);
        Require(sums.storage==expected.storage && input.storage==source_before,"terminal scalar layers plus XOR");
        // A second identical accumulation cancels in GF(2); this also detects
        // accidentally overwriting the accumulator or reading mutated inputs.
        LeoGFNIFinalAccumulate(inputs.data(),sums.pointers.data(),distance,a,b,c,bytes);
        Require(sums.storage==sums_before && input.storage==source_before,"repeat accumulation cancellation");
        input.Prefix(); work.Prefix(); sums.Prefix(); expected.Prefix(); ++cases;
    };
    for (unsigned distance : {1U,4U})
    for (size_t bytes : lengths)
    for (size_t offset : {size_t(0),size_t(17)})
    for (unsigned mask=0; mask<8; ++mask) check(distance,bytes,offset,mask);
    for (size_t bytes : {size_t(32768),size_t(32770)})
    for (size_t offset : {size_t(0),size_t(17)})
    for (unsigned mask=0; mask<8; ++mask) check(64,bytes,offset,mask);
    std::printf("terminal GFNI kernel: %u scalar-layer/XOR, zero-skew, exact-end, unaligned, readonly-input and cancellation cases passed\n",cases);
}

void Public(unsigned index)
{
    using namespace gfni_terminal;
    struct Shape { unsigned k,r; size_t bytes; leo2_backend backend; bool partial; unsigned changes; };
    const Shape shapes[] = {
        {1000,200,65536,LEO2_BACKEND_AUTO,false,2},
        {1000,200,65535,LEO2_BACKEND_AUTO,false,0},
        {1000,200,65537,LEO2_BACKEND_AUTO,false,0},
        {1000,200,65536,LEO2_BACKEND_GFNI,false,2},
        {1000,200,65538,LEO2_BACKEND_GFNI,false,2},
        {1000,200,65536,LEO2_BACKEND_AUTO,true,0},
        {1000,200,65536,LEO2_BACKEND_GFNI,true,0},
        {1000,199,65536,LEO2_BACKEND_AUTO,false,0},
        {999,200,65536,LEO2_BACKEND_GFNI,false,0},
        {1000,201,65536,LEO2_BACKEND_GFNI,false,0},
        {1000,200,32768,LEO2_BACKEND_GFNI,false,0},
        {1000,200,65536,LEO2_BACKEND_AVX512,false,0},
        {1000,200,65534,LEO2_BACKEND_AUTO,false,0},
        {1000,200,65538,LEO2_BACKEND_AUTO,false,0}
    };
    Require(index<14,"public index");
    const Shape& shape=shapes[index];
    leo2_context_options options={}; options.struct_size=sizeof(options);
    options.backend=shape.backend; options.thread_count=1;
    leo2_context* raw_context=NULL;
    Require(leo2_context_create(&options,&raw_context)==LEO2_SUCCESS,"context");
    std::unique_ptr<leo2_context,decltype(&leo2_context_destroy)> context(raw_context,leo2_context_destroy);
    leo2_codec* raw_codec=NULL;
    Require(leo2_codec_create(context.get(),shape.k,shape.r,LEO2_PROFILE_LEGACY_HIGH_V1,
        LEO2_FIELD_GF16,NULL,&raw_codec)==LEO2_SUCCESS,"codec");
    std::unique_ptr<leo2_codec,decltype(&leo2_codec_destroy)> codec(raw_codec,leo2_codec_destroy);
    if (shape.bytes&1)
    {
        for (bool enabled : {false,true})
        {
            Reset(enabled); size_t size=77;
            Require(leo2_encode_scratch_size(codec.get(),shape.bytes,&size)==LEO2_UNSUPPORTED && size==0,
                "native odd query");
            Require(leo2_encode(codec.get(),shape.bytes,NULL,NULL,NULL,0)==LEO2_UNSUPPORTED && Get().calls==0,
                "native odd rejected before transform");
        }
        std::printf("terminal public %u: native odd rejected in both modes\n",index); return;
    }
    Rows source(17,shape.bytes,17,20260906),control(shape.r,shape.bytes,17,29),candidate(shape.r,shape.bytes,17,29);
    const auto source_before=source.storage;
    const auto omitted=candidate.storage[17];
    std::vector<const void*> inputs(shape.k);
    for (unsigned i=0; i<shape.k; ++i) inputs[i]=source.pointers[i%17];
    if (shape.partial) control.pointers[17]=candidate.pointers[17]=NULL;
    size_t size=0;
    Require(leo2_encode_scratch_size(codec.get(),shape.bytes,&size)==LEO2_SUCCESS,"scratch size");
    void* raw_scratch=NULL;
    Require(posix_memalign(&raw_scratch,64,size)==0,"scratch allocation");
    std::unique_ptr<void,decltype(&std::free)> scratch(raw_scratch,std::free);
    if (index==0)
    {
        Reset(true);
        Require(leo2_encode(codec.get(),shape.bytes,inputs.data(),candidate.pointers.data(),scratch.get(),size-1)
            ==LEO2_SCRATCH_TOO_SMALL && Get().calls==0,"short scratch guard");
        void* saved=candidate.pointers[0]; candidate.pointers[0]=source.pointers[0];
        Require(leo2_encode(codec.get(),shape.bytes,inputs.data(),candidate.pointers.data(),scratch.get(),size)
            ==LEO2_OVERLAP && Get().calls==0,"overlap guard");
        candidate.pointers[0]=saved;
        Require(candidate.storage==control.storage && source.storage==source_before,"rejected output unchanged");
    }
    unsigned passes=0;
    for (bool enabled : {false,true})
    {
        Reset(enabled); std::memset(scratch.get(),enabled ? 0xa7 : 0x5a,size);
        Rows& output=enabled ? candidate : control;
        Require(leo2_encode(codec.get(),shape.bytes,inputs.data(),output.pointers.data(),scratch.get(),size)
            ==LEO2_SUCCESS,"public encode");
        Require(Get().calls>0 && Get().matches==shape.changes && Get().changed==(enabled ? shape.changes : 0U),
            "directed predicate counts");
        if (!enabled) passes=Get().calls;
        Require(Get().calls==passes,"pass count preserved");
    }
    Require(control.storage==candidate.storage && source.storage==source_before,"parity/input identity");
    if (shape.partial) Require(candidate.storage[17]==omitted,"omitted output unchanged");
    source.Prefix(); control.Prefix(); candidate.Prefix();
    std::printf("terminal public %u: K%u/R%u/B%zu backend%d partial%d calls%u changed%u parity exact\n",
        index,shape.k,shape.r,shape.bytes,shape.backend,shape.partial,passes,Get().changed);
}
}

int main(int argc,char** argv)
{
    try {
        omp_set_dynamic(0); omp_set_num_threads(1);
        if (argc==2 && std::strcmp(argv[1],"--kernel")==0) { Predicate(); Kernels(); }
        else {
            Require(argc==3 && std::strcmp(argv[1],"--public")==0,"--kernel or --public [0..13]");
            char* end=NULL; const unsigned long index=std::strtoul(argv[2],&end,10);
            Require(argv[2][0] && *end==0 && index<14,"public index"); Public(static_cast<unsigned>(index));
        }
        return 0;
    }
    catch (const std::exception& error) { std::fprintf(stderr,"terminal unit: %s\n",error.what()); return 1; }
}

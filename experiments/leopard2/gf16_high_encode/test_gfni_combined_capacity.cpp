// No transform or kernel execution; exact four-mode state for both capacities.
#include "gfni_combined.h"
#include <cstdio>
#include <stdexcept>

namespace {
void Require(bool value,const char* message)
{ if (!value) throw std::runtime_error(message); }
}

int main()
{
    try
    {
        using namespace gfni_combined;
        leopard::backend::Ops ops={}; ops.kind=LEO2_BACKEND_GFNI;
        leopard2_internal::SparseForwardPlanBatchView sparse={};
        for (unsigned mode=0; mode<4; ++mode)
        {
            Reset(mode);
            for (unsigned i=0; i<kCapacity; ++i)
            {
                Require(LeoGFNICombinedExperiment(ops,32768,65536,1000,200,200,256,&sparse)==mode,
                    "mode selection");
                const Call& c=Get().records[i];
                Require(c.kind==LEO2_BACKEND_GFNI && c.k==1000 && c.r==200 && c.requested==200 &&
                    c.side==256 && c.sparse_blocks==0 && c.sparse_present && c.bytes==32768 &&
                    c.source_policy==65536,"stored record");
            }
            bool refused=false;
            try { LeoGFNICombinedExperiment(ops,32768,65536,1000,200,200,256,&sparse); }
            catch (const std::runtime_error&) { refused=true; }
            Require(refused && Get().mode==mode && Get().calls==kCapacity && Get().matches==kCapacity &&
                Get().first==((mode&1U) ? kCapacity : 0U) &&
                Get().terminal==((mode&2U) ? kCapacity : 0U),"atomic capacity refusal");
            Reset(mode);
            Require(Get().mode==mode && Get().calls==0 && Get().matches==0 &&
                Get().first==0 && Get().terminal==0,"reset");
            Require(!LeoGFNICombinedExperiment(ops,32768,65536,1000,199,199,256,&sparse),"neighbor");
            Require(Get().calls==1 && Get().matches==0 && Get().first==0 && Get().terminal==0,"neighbor record");
        }
        std::printf("combined capacity %u: four modes, overflow, exact records, reset, neighbor passed\n",kCapacity);
        return 0;
    }
    catch (const std::exception& error) { std::fprintf(stderr,"combined capacity: %s\n",error.what()); return 1; }
}

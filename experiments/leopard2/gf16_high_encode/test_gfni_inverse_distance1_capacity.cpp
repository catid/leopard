// No transform/kernel execution; bound and exact state checks for both capacities.
#include "gfni_inverse_distance1.h"
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
        using namespace gfni_inverse_distance1;
        leopard::backend::Ops ops={}; ops.kind=LEO2_BACKEND_GFNI;
        leopard2_internal::SparseForwardPlanBatchView sparse={};
        for (bool enabled : {false,true})
        {
            Reset(enabled);
            for (unsigned i=0; i<kCapacity; ++i)
            {
                Require(LeoGFNIInverseDistance1Experiment(ops,32768,65536,1000,200,200,256,&sparse)==enabled,
                    "mode selection");
                const Call& c=Get().records[i];
                Require(c.kind==LEO2_BACKEND_GFNI && c.k==1000 && c.r==200 && c.requested==200 &&
                    c.side==256 && c.sparse_blocks==0 && c.sparse_present && c.bytes==32768 &&
                    c.source_policy==65536,"stored record");
            }
            bool refused=false;
            try { LeoGFNIInverseDistance1Experiment(ops,32768,65536,1000,200,200,256,&sparse); }
            catch (const std::runtime_error&) { refused=true; }
            Require(refused && Get().calls==kCapacity && Get().matches==kCapacity &&
                Get().changed==(enabled ? kCapacity : 0U),"atomic capacity refusal");
            Reset(enabled);
            Require(Get().calls==0 && Get().matches==0 && Get().changed==0,"reset");
            Require(!LeoGFNIInverseDistance1Experiment(ops,32768,65536,1000,199,199,256,&sparse),"neighbor");
            Require(Get().calls==1 && Get().matches==0 && Get().changed==0,"neighbor record");
        }
        std::printf("distance-one capacity %u: both modes, overflow, exact records, reset, neighbor passed\n",kCapacity);
        return 0;
    }
    catch (const std::exception& error) { std::fprintf(stderr,"distance-one capacity: %s\n",error.what()); return 1; }
}

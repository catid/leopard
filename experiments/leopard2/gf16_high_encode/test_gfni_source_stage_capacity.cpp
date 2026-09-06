// No codec linked or executed: test fail-before-delegation at both capacities.
#include "gfni_source_stage_wrap.cpp"

namespace {
unsigned delegated = 0;
uint64_t observed_policy = 0;
void Check(bool value)
{
    if (!value) throw std::runtime_error("capacity check failed");
}
}

extern "C" void RealSourceStage(const leopard::backend::Ops&, uint64_t,
    uint64_t policy, unsigned, unsigned, unsigned, unsigned,
    const void* const*, void**, const leopard2_internal::SparseForwardPlanBatchView*)
{
    ++delegated;
    observed_policy = policy;
}

int main()
{
    using namespace gfni_source_stage_probe;
    leopard::backend::Ops ops = {};
    ops.kind = LEO2_BACKEND_GFNI;
    for (bool enabled : {false, true})
    {
        Reset(enabled);
        delegated = 0;
        for (unsigned i = 0; i < kCapacity; ++i)
        {
            WrappedSourceStage(ops, 32768, 65536, 1000, 200, 200, 256,
                NULL, NULL, NULL);
            Check(Get().calls == i + 1 && delegated == i + 1 &&
                  observed_policy == (enabled ? 16384U : 65536U));
        }
        bool rejected = false;
        try
        {
            WrappedSourceStage(ops, 32768, 65536, 1000, 200, 200, 256,
                NULL, NULL, NULL);
        }
        catch (const std::runtime_error&) { rejected = true; }
        Check(rejected && delegated == kCapacity && Get().calls == kCapacity &&
              Get().matches == kCapacity && Get().changed == (enabled ? kCapacity : 0));
        Reset(enabled);
        Check(Get().calls == 0 && Get().matches == 0 && Get().changed == 0);
        for (unsigned i = 0; i < kCapacity; ++i)
            Check(Get().records[i].source_policy == 0);
        WrappedSourceStage(ops, 32768, 65536, 1000, 199, 199, 256,
            NULL, NULL, NULL);
        Check(Get().calls == 1 && Get().matches == 0 && Get().changed == 0 &&
              observed_policy == 65536);
    }
    std::printf("source-stage capacity %u: both modes, overflow delegation, reset, neighbor passed\n", kCapacity);
}

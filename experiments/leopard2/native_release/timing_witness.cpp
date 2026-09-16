// Only synthetic-clock links interpose the public call. Production links do not.
#include <cstdint>
#ifdef LEO_NATIVE_RELEASE_BASELINE
#include "leopard.h"
#else
#include "leopard2.h"
#endif
namespace { uint64_t count = 0; }
extern "C" uint64_t NativeWitnessCalls() { return count; }
#ifdef LEO_NATIVE_RELEASE_BASELINE
extern "C" LeopardResult __real_leo_encode(uint64_t, unsigned, unsigned, unsigned,
    const void* const*, void**);
extern "C" LeopardResult __wrap_leo_encode(uint64_t bytes, unsigned k, unsigned r,
    unsigned n, const void* const* originals, void** work)
{
    ++count;
    return __real_leo_encode(bytes, k, r, n, originals, work);
}
#else
extern "C" leo2_result __real_leo2_encode(const leo2_codec*, uint64_t,
    const void* const*, void* const*, void*, size_t);
extern "C" leo2_result __wrap_leo2_encode(const leo2_codec* codec, uint64_t bytes,
    const void* const* originals, void* const* recovery, void* scratch, size_t scratch_bytes)
{
    ++count;
    return __real_leo2_encode(codec, bytes, originals, recovery, scratch, scratch_bytes);
}
#endif

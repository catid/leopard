/*
    Experimental page-relative cache coloring for FF8 XOR benchmark buffers.

    This is deliberately an allocation-side helper: It never rewrites caller
    pointers and does not change the native eight-plane shard format.  Each
    exposed shard is still one contiguous buffer of exactly BufferBytes usable
    bytes.  The extra allocation space only permits selecting its address bits
    6..11 (the 64 conventional 64-byte L1D set colors within a 4 KiB period).
*/

#ifndef LEOPARD_TEST_FF8XOR_CACHE_COLORED_BUFFERS_H
#define LEOPARD_TEST_FF8XOR_CACHE_COLORED_BUFFERS_H

#include "../LeopardCommon.h"

#include <limits>
#include <new>
#include <stdint.h>
#include <vector>

namespace leopard_ff8xor_test {

static const unsigned kCacheLineBytes = 64;
static const unsigned kCacheColorCount = 64;
static const unsigned kCacheColorPeriod =
    kCacheLineBytes * kCacheColorCount;

// Decode first copies original shard i into transform work position m + i,
// where m is a power of two and m + original_count <= 256.  This salt was
// selected by exhaustive search over that complete domain so source and work
// colors differ for every copy.  Keep the exhaustive invariant check in the
// benchmark in sync with this value.
static const unsigned kDecodeWorkCacheColorSalt = 5;

// A rank-six linear map from the eight transform-index bits to the six L1D
// color bits.  Each column is nonzero, so flipping any one index bit (exactly
// what a radix-2 FFT partner does) must change the color.  The first six bits
// cover every color; the final two masks keep the property through n=256.
static inline unsigned TransformCacheColor(unsigned transform_index)
{
    static const unsigned kIndexBitMasks[8] = {
        1, 2, 4, 8, 16, 32, 21, 42
    };
    unsigned color = 0;
    for (unsigned bit = 0; bit < 8; ++bit)
    {
        if (transform_index & (1U << bit))
            color ^= kIndexBitMasks[bit];
    }
    return color;
}

static inline unsigned SaltedTransformCacheColor(
    unsigned transform_index,
    unsigned color_salt)
{
    return TransformCacheColor(transform_index) ^
        (color_salt & (kCacheColorCount - 1));
}

static inline size_t TransformCacheColorOffset(
    unsigned transform_index,
    unsigned color_salt)
{
    return static_cast<size_t>(SaltedTransformCacheColor(
        transform_index, color_salt)) * kCacheLineBytes;
}

// All eight plane starts select the same set exactly when the plane stride is
// a multiple of the 4 KiB color period.  Smaller/non-aliasing layouts benefit
// from compact allocator placement and should not pay the extra-page/TLB cost.
static inline bool PlaneStartsFullyCacheAlias(size_t buffer_bytes)
{
    return buffer_bytes >= 8 &&
        (buffer_bytes / 8) % kCacheColorPeriod == 0;
}

class BufferSet
{
public:
    BufferSet() {}
    ~BufferSet() { Clear(); }

    BufferSet(const BufferSet&) = delete;
    BufferSet& operator=(const BufferSet&) = delete;

    bool Allocate(
        unsigned count,
        size_t bytes,
        bool cache_colored = false,
        unsigned transform_index_base = 0,
        unsigned color_salt = 0)
    {
        Clear();
        try
        {
            Entries.reserve(count);
            Pointers.reserve(count);
            for (unsigned index = 0; index < count; ++index)
            {
                Entry entry;
                if (!AllocateOne(
                        bytes,
                        cache_colored,
                        transform_index_base + index,
                        color_salt,
                        entry))
                {
                    Clear();
                    return false;
                }
                Entries.push_back(entry);
                Pointers.push_back(entry.Exposed);
            }
        }
        catch (const std::bad_alloc&)
        {
            Clear();
            return false;
        }
        return true;
    }

    void Clear()
    {
        for (size_t index = 0; index < Entries.size(); ++index)
            leopard::SIMDSafeFree(Entries[index].Owner);
        Entries.clear();
        Pointers.clear();
    }

    bool TransferFirst(unsigned count, BufferSet& destination)
    {
        if (count > Entries.size())
            return false;
        destination.Clear();
        try
        {
            destination.Entries.reserve(count);
            destination.Pointers.reserve(count);
        }
        catch (const std::bad_alloc&)
        {
            return false;
        }
        for (unsigned index = 0; index < count; ++index)
        {
            destination.Entries.push_back(Entries[index]);
            destination.Pointers.push_back(Entries[index].Exposed);
            Entries[index].Owner = NULL;
            Entries[index].Exposed = NULL;
        }
        Clear();
        return true;
    }

    uint8_t* operator[](size_t index)
    {
        return Entries[index].Exposed;
    }

    const uint8_t* operator[](size_t index) const
    {
        return Entries[index].Exposed;
    }

    void** Data() { return Pointers.empty() ? NULL : &Pointers[0]; }
    unsigned Count() const { return static_cast<unsigned>(Entries.size()); }

    bool ValidateColoredAllocation(
        unsigned index,
        size_t expected_bytes,
        unsigned expected_transform_index,
        unsigned expected_color_salt) const
    {
        if (index >= Entries.size())
            return false;
        const Entry& entry = Entries[index];
        if (!entry.Colored || entry.Bytes != expected_bytes ||
            entry.Owner == NULL || entry.Exposed == NULL)
            return false;

        const uintptr_t owner = reinterpret_cast<uintptr_t>(entry.Owner);
        const uintptr_t exposed = reinterpret_cast<uintptr_t>(entry.Exposed);
        const uintptr_t end = owner + entry.OwnerBytes;
        if (end < owner || exposed < owner || exposed > end ||
            expected_bytes > end - exposed)
            return false;

        const size_t expected_offset = TransformCacheColorOffset(
            expected_transform_index, expected_color_salt);
        return exposed % kCacheColorPeriod == expected_offset &&
            exposed % kCacheLineBytes == 0;
    }

private:
    struct Entry
    {
        uint8_t* Owner;
        uint8_t* Exposed;
        size_t OwnerBytes;
        size_t Bytes;
        bool Colored;

        Entry()
            : Owner(NULL)
            , Exposed(NULL)
            , OwnerBytes(0)
            , Bytes(0)
            , Colored(false)
        {
        }
    };

    static bool AllocateOne(
        size_t bytes,
        bool cache_colored,
        unsigned transform_index,
        unsigned color_salt,
        Entry& entry)
    {
        if (!cache_colored)
        {
            entry.Owner = leopard::SIMDSafeAllocate(bytes);
            entry.Exposed = entry.Owner;
            entry.OwnerBytes = bytes;
            entry.Bytes = bytes;
            entry.Colored = false;
            return entry.Owner != NULL;
        }

        const size_t maximum_color_offset =
            static_cast<size_t>(kCacheColorCount - 1) * kCacheLineBytes;
        const size_t overhead =
            static_cast<size_t>(kCacheColorPeriod - 1) + maximum_color_offset;
        if (bytes > (std::numeric_limits<size_t>::max)() - overhead)
            return false;

        entry.OwnerBytes = bytes + overhead;
        entry.Owner = leopard::SIMDSafeAllocate(entry.OwnerBytes);
        if (!entry.Owner)
            return false;

        const uintptr_t owner = reinterpret_cast<uintptr_t>(entry.Owner);
        const uintptr_t period_aligned =
            (owner + kCacheColorPeriod - 1) &
            ~static_cast<uintptr_t>(kCacheColorPeriod - 1);
        entry.Exposed = reinterpret_cast<uint8_t*>(period_aligned +
            TransformCacheColorOffset(transform_index, color_salt));
        entry.Bytes = bytes;
        entry.Colored = true;
        return true;
    }

    std::vector<Entry> Entries;
    std::vector<void*> Pointers;
};

} // namespace leopard_ff8xor_test

#endif // LEOPARD_TEST_FF8XOR_CACHE_COLORED_BUFFERS_H

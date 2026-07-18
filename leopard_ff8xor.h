/*
    Experimental native plane-sliced FF8 XOR-circuit API.

    This API is intentionally separate from leopard.h.  It uses Leopard's
    existing FF8 algebra and transform schedule, but each B-byte shard is laid
    out as eight contiguous B/8-byte coordinate-bit planes.  It does not
    transpose packed Leopard shards at the API boundary.
*/

#ifndef CAT_LEOPARD_FF8_XOR_H
#define CAT_LEOPARD_FF8_XOR_H

#include "leopard.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
    Return the number of B-byte work buffers required by leo_ff8xor_encode().
    The count rules are identical to leo_encode_work_count().
*/
LEO_EXPORT unsigned leo_ff8xor_encode_work_count(
    unsigned original_count,
    unsigned recovery_count);

/*
    Encode native plane-sliced shards.

    Each shard is one contiguous B-byte allocation split into eight B/8-byte
    planes.  At group byte g, bit lane l of plane j is coordinate bit j of
    field element 8*g+l.  The first recovery_count work buffers receive parity.

    leo_init() must succeed first.  buffer_bytes must be a positive multiple
    of 64, and the padded FF8 transform size must not exceed 256.
*/
LEO_EXPORT LeopardResult leo_ff8xor_encode(
    uint64_t buffer_bytes,
    unsigned original_count,
    unsigned recovery_count,
    unsigned work_count,
    const void* const* original_data,
    void** work_data);

/*
    Return the number of B-byte work buffers required by leo_ff8xor_decode().
    The count rules are identical to leo_decode_work_count().
*/
LEO_EXPORT unsigned leo_ff8xor_decode_work_count(
    unsigned original_count,
    unsigned recovery_count);

/*
    Decode native plane-sliced shards.  Missing shard pointers are NULL.
    Recovered original i is returned in work_data[i], matching leo_decode().
*/
LEO_EXPORT LeopardResult leo_ff8xor_decode(
    uint64_t buffer_bytes,
    unsigned original_count,
    unsigned recovery_count,
    unsigned work_count,
    const void* const* original_data,
    const void* const* recovery_data,
    void** work_data);

#ifdef __cplusplus
}
#endif

#endif // CAT_LEOPARD_FF8_XOR_H

// Production GFNI backend member.  Re-includes the AVX2 implementation with
// the GFNI variant kernels enabled and its own backend identity, mirroring
// the Leopard2BackendAVX512.cpp pattern.  This translation unit is compiled
// with -mavx2 -mgfni and owns its affine multiplication tables outright; it
// never shares the AVX2 member's nibble tables.
#define LEO2_GFNI_MEMBER 1
#define LEO2_GFNI_VARIANT 1
#include "Leopard2BackendAVX2.cpp"

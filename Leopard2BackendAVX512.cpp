/*
    Copyright (c) 2017 Christopher A. Taylor.  All rights reserved.
    See Leopard2Backend.h for the full BSD license notice.
*/

// Recompile the reviewed AVX2 kernels with AVX-512F/BW/VL enabled.  Their
// external data width remains 256 bits; the wider ISA is used for the extra
// architectural vector registers, which prevent GF16 table spills.  Keeping
// one implementation source also keeps scalar tails and fixed-multiplier
// semantics identical between the two runtime-qualified backends.
#define LEO2_AVX512_VARIANT 1
#include "Leopard2BackendAVX2.cpp"

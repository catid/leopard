/*
    Copyright (c) 2017 Christopher A. Taylor.  All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the conditions in the Leopard-RS
    BSD license are met.  See the repository LICENSE file.
*/

#ifndef LEOPARD2_LEGACY_GOLDEN_VECTORS_H
#define LEOPARD2_LEGACY_GOLDEN_VECTORS_H

#include <stddef.h>
#include <stdint.h>

namespace leopard2_legacy_golden {

/*
    Frozen old-Leopard wire vectors generated with public leo_encode() at the
    repository baseline 3492b9740b55d85b05df97164d4a1f14c6f9e9ae.  The
    legacy encoder source's last commit there is
    36427dd25bf665f410525a03a1f0a0ea9275150b.  The compatibility test does not
    include leopard.h or invoke that encoder.

    Original data generator (unsigned 64-bit arithmetic wraps modulo 2^64):

        state = seed
        for shard = 0..K-1, then byte_offset = 0..shard_bytes-1:
            state += 0x9e3779b97f4a7c15
            z = state
            z = (z xor (z >> 30)) * 0xbf58476d1ce4e5b9
            z = (z xor (z >> 27)) * 0x94d049bb133111eb
            z = z xor (z >> 31)
            original[shard][byte_offset] = z >> 56

    parity_hex is the raw transmitted byte layout: recovery shard 0 bytes
    0..63, followed by recovery shard 1, and so on.  Hex pairs are in
    increasing memory-address order.  They are never parsed as host integers.
    In particular, GF16 vectors preserve legacy Leopard's 64-byte ALTMAP tile
    and require no host-endian conversion.

    SHA-256 over the 960 decoded parity bytes, concatenated in kGoldenVectors
    order, is:

        6973b1ec2b02ab7d171b8c7ed3a8281ea744896ca66805c2544aef4884bebafe

    Per-vector decoded parity SHA-256 values are:

        gf8-k1-r1:     c54604beb0af7771249f48a82031730bf14832bc608c6c3ac16b924d3082568c
        gf8-k3-r2:     866e2470b0c1a55ef28c2111d9a792f4d12f48ff04805b67f2d665e91db94666
        gf8-k5-r3:     46a6ad12bc4879129b29ca178c2f0da05ec1c562adc86e17b2449ebffc3ceb0b
        gf8-k9-r7:     0bed3f6b54759fa91d012efc09f9404b08452a8278b37045aa4769d7a4ef1024
        gf16-k255-r2:  11f506bfec4e243dd21eb076c6b5af98b749d872486732f5552846a97fed9788
*/

struct GoldenVector
{
    const char* name;
    uint32_t original_count;
    uint32_t recovery_count;
    uint32_t field_bits;
    uint32_t shard_bytes;
    uint64_t seed;
    const uint32_t* missing_originals;
    size_t missing_count;
    const char* parity_hex;
};

static const char kGf8K1R1Parity[] =
    "5f038a8dbcdf7838f9cf521bf291cf789fce5814f841afb72ada41b10699096741824682f77a826d7def1ac7252e4776"
    "4f2fe110f915da0441fd47d2bb09584a";

static const char kGf8K3R2Parity[] =
    "bff0a17a71e04eec2c8a017e4f3b18a4b4649621e8add0b5bee2649a110fed61271d8d19eab0132e85a02c8f97022a27"
    "605dcf69af6c9dd3bcd1fd0f36f900f2f12c5d13ad35fb3b3a08ba38d585b14487aaefb047efbb67ae09fba310c0e0a0"
    "fa853a9da933dbad2597c5aface090b444e78e03cf3082b9b250d0ca51748010";

static const char kGf8K5R3Parity[] =
    "3dce972f3f492cc5e481370f9099fb08f7671d4661c6f5ad098181ee68fb085f5977d5743b9a8a0df4d973994235eb0f"
    "c77f4b4c60984e5e0c8f87d107e970774da1e8013f609adac295245d478a263a15f3face56148c42d9765af7f952bdff"
    "2fcd806e41daeef0fda06612d84a29054ce4c2f38d1bfefb7a6ff8aaa3746c1fbbd22cf613de0e1d7f7e1f05e207d5b5"
    "e5a97f9747ed485f70a16034da07d015c613253b1eb54e88912ae9ae3b2e2025ef5b55c4d550526d1d076149ed7cd29e";

static const char kGf8K9R7Parity[] =
    "15090627474ef8d9569913d1b0c2de4f37a5b6614b44a309650289d87ab168aea844bab91aa58529651e5aec18cfb41a"
    "5607f9a253c5b6dd3a46622290d2a9d72976eea9ef11585be4514bc878d8b8cc16493837329fc94f4b3836beff01b221"
    "30d64b80beb3e18b7e6858ae82d1af8a12bd1aed0bea773afa7551f68817f372609995ca59a97cd852bebb6c085c748d"
    "a7c497ce86532b09a09f23bf02800c07d9b4e7499f39c87a289d0a33063bab77006073be793bbc630072388439f6a812"
    "43f62e1d1c4b182dc71b390a931a876fe908bbd125f1d5f02c1676e22d256aafeda65ddbc3ab39f85a200645630adba8"
    "ffcead6b4589bd280807e4ca63e2781eaf8e63e5cc35900dea7056ecb3ce4246ed8762ed94de96ee700c2a92ffc7f45e"
    "dc30da9f288b433e078c2595087506141e62d0f07f63e35fc8850ab4e089af518039d639521460cf3433184fe0ec5877"
    "775d7e364b9bca16c1c886153057da36cd29fb34744a5fad303de398cf4d2897bde687f43cb91d564053c74c35546371"
    "19bddfe947e753a980455976d3783f1a65aaf7db965587a097568703ba381de68bfa40999b67e9183b608229a63d11a2"
    "aaa7ba791f81f00ced9622330be9a903";

static const char kGf16K255R2Parity[] =
    "c7b7e4235a4e68d093464c14b982ae9b68e480f31bfa6c70c106ca5e516004bbc258189aaae91b2a8b0945fe3ed7f94b"
    "0ed5050fbfd5fd33899d226c7689172154eeb72d8c8feae481ae63004d5fe1332ae85e729dc99e9c8df5797fc53f334a"
    "a62e1ac08110ec0fa3ed8d52351d16132359e429ca1ecc7cc24b080898270c67";

static const uint32_t kLossK1R1[] = { 0 };
static const uint32_t kLossK3R2[] = { 0, 2 };
static const uint32_t kLossK5R3[] = { 0, 2, 4 };
static const uint32_t kLossK9R7[] = { 0, 1, 2, 3, 5, 6, 8 };
static const uint32_t kLossK255R2[] = { 0, 254 };

static const GoldenVector kGoldenVectors[] = {
    { "gf8-k1-r1", 1, 1, 8, 64, UINT64_C(0x4c32474638010101),
      kLossK1R1, sizeof(kLossK1R1) / sizeof(kLossK1R1[0]), kGf8K1R1Parity },
    { "gf8-k3-r2", 3, 2, 8, 64, UINT64_C(0x4c32474638030202),
      kLossK3R2, sizeof(kLossK3R2) / sizeof(kLossK3R2[0]), kGf8K3R2Parity },
    { "gf8-k5-r3", 5, 3, 8, 64, UINT64_C(0x4c32474638050303),
      kLossK5R3, sizeof(kLossK5R3) / sizeof(kLossK5R3[0]), kGf8K5R3Parity },
    { "gf8-k9-r7", 9, 7, 8, 64, UINT64_C(0x4c32474638090707),
      kLossK9R7, sizeof(kLossK9R7) / sizeof(kLossK9R7[0]), kGf8K9R7Parity },
    { "gf16-k255-r2", 255, 2, 16, 64, UINT64_C(0x4c324746ffff0202),
      kLossK255R2, sizeof(kLossK255R2) / sizeof(kLossK255R2[0]), kGf16K255R2Parity }
};

static const size_t kGoldenVectorCount =
    sizeof(kGoldenVectors) / sizeof(kGoldenVectors[0]);

} // namespace leopard2_legacy_golden

#endif // LEOPARD2_LEGACY_GOLDEN_VECTORS_H

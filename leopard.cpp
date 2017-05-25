/*
    Copyright (c) 2017 Christopher A. Taylor.  All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.
    * Neither the name of Leopard-RS nor the names of its contributors may be
      used to endorse or promote products derived from this software without
      specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
    ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
    LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
    POSSIBILITY OF SUCH DAMAGE.
*/

#include "leopard.h"
#include "FecalEncoder.h"
#include "FecalDecoder.h"

extern "C" {


//------------------------------------------------------------------------------
// Initialization API

static bool m_Initialized = false;

FECAL_EXPORT int fecal_init_(int version)
{
    if (version != FECAL_VERSION)
        return Fecal_InvalidInput;

    if (0 != gf256_init())
        return Fecal_Platform;

    m_Initialized = true;
    return Fecal_Success;
}


//------------------------------------------------------------------------------
// Encoder API

FECAL_EXPORT FecalEncoder fecal_encoder_create(unsigned input_count, void* const * const input_data, uint64_t total_bytes)
{
    if (input_count <= 0 || !input_data || total_bytes < input_count)
    {
        FECAL_DEBUG_BREAK; // Invalid input
        return nullptr;
    }

    FECAL_DEBUG_ASSERT(m_Initialized); // Must call fecal_init() first
    if (!m_Initialized)
        return nullptr;

    fecal::Encoder* encoder = new(std::nothrow) fecal::Encoder;
    if (!encoder)
    {
        FECAL_DEBUG_BREAK; // Out of memory
        return nullptr;
    }

    if (Fecal_Success != encoder->Initialize(input_count, input_data, total_bytes))
    {
        delete encoder;
        return nullptr;
    }

    return reinterpret_cast<FecalEncoder>( encoder );
}

FECAL_EXPORT int fecal_encode(FecalEncoder encoder_v, FecalSymbol* symbol)
{
    fecal::Encoder* encoder = reinterpret_cast<fecal::Encoder*>( encoder_v );
    if (!encoder || !symbol)
        return Fecal_InvalidInput;

    return encoder->Encode(*symbol);
}

FECAL_EXPORT void fecal_free(void* codec_v)
{
    if (codec_v)
    {
        fecal::ICodec* icodec = reinterpret_cast<fecal::ICodec*>( codec_v );
        delete icodec;
    }
}


//------------------------------------------------------------------------------
// Decoder API

FECAL_EXPORT FecalDecoder fecal_decoder_create(unsigned input_count, uint64_t total_bytes)
{
    if (input_count <= 0 || total_bytes < input_count)
    {
        FECAL_DEBUG_BREAK; // Invalid input
        return nullptr;
    }

    FECAL_DEBUG_ASSERT(m_Initialized); // Must call fecal_init() first
    if (!m_Initialized)
        return nullptr;

    fecal::Decoder* decoder = new(std::nothrow) fecal::Decoder;
    if (!decoder)
    {
        FECAL_DEBUG_BREAK; // Out of memory
        return nullptr;
    }

    if (Fecal_Success != decoder->Initialize(input_count, total_bytes))
    {
        delete decoder;
        return nullptr;
    }

    return reinterpret_cast<FecalDecoder>( decoder );
}

FECAL_EXPORT int fecal_decoder_add_original(FecalDecoder decoder_v, const FecalSymbol* symbol)
{
    fecal::Decoder* decoder = reinterpret_cast<fecal::Decoder*>( decoder_v );
    if (!decoder || !symbol)
        return Fecal_InvalidInput;

    return decoder->AddOriginal(*symbol);
}

FECAL_EXPORT int fecal_decoder_add_recovery(FecalDecoder decoder_v, const FecalSymbol* symbol)
{
    fecal::Decoder* decoder = reinterpret_cast<fecal::Decoder*>( decoder_v );
    if (!decoder || !symbol)
        return Fecal_InvalidInput;

    return decoder->AddRecovery(*symbol);
}

FECAL_EXPORT int fecal_decode(FecalDecoder decoder_v, RecoveredSymbols* symbols)
{
    fecal::Decoder* decoder = reinterpret_cast<fecal::Decoder*>( decoder_v );
    if (!decoder || !symbols)
        return Fecal_InvalidInput;

    return decoder->Decode(*symbols);
}

FECAL_EXPORT int fecal_decoder_get(FecalDecoder decoder_v, unsigned input_index, FecalSymbol* symbol)
{
    fecal::Decoder* decoder = reinterpret_cast<fecal::Decoder*>( decoder_v );
    if (!decoder || !symbol)
        return Fecal_InvalidInput;

    return decoder->GetOriginal(input_index, *symbol);
}


} // extern "C"

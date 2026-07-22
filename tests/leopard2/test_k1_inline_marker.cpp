#include "Leopard2Direct.h"

#include <iostream>

int main()
{
#if defined(LEO2_EXPECT_K1_DECODE_VALIDATOR_INLINE)
    const bool expected = true;
#else
    const bool expected = false;
#endif
    const bool enabled = leopard2_internal::
        K1DecodeValidatorInlineExperimentEnabled();
    if (enabled != expected)
    {
        std::cerr << "K1 decode-validator inline experiment marker mismatch\n";
        return 1;
    }
    return 0;
}

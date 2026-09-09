"""Exact, experiment-only source overlay; never modifies production inputs.

Consumes the already qualified static overlay, not a guessed source revision.
Output keeps its two arithmetic schedules and scalar tails; wrappers choose
template bodies once per pair outside the vector loop. No codec execution.
"""
import hashlib

STATIC_SHA = 'be5a07759cd0968299ebf243be77b299afda78065de704dcb4a319ca551ba205'
ORDINARY = '!defined(LEO2_GFNI_VARIANT) && !defined(LEO2_AVX512_VARIANT)'


def once(text, old, new):
    if text.count(old)!=1: raise ValueError('ambiguous or missing overlay anchor')
    return text.replace(old,new,1)


def overlay(text):
    if hashlib.sha256(text.encode()).hexdigest()!=STATIC_SHA:
        raise ValueError('unqualified static source')
    text = once(text,'#include <new>\n',
        '#include <new>\n#if '+ORDINARY+'\n#include "avx2_adjacent_control.h"\n#endif\n')
    text = once(text,'template<bool Inverse>\nstatic LEO_FORCE_INLINE void AVX2FF16Butterfly2Prepared(',
        'template<bool Inverse, bool Scheduled>\nstatic LEO_FORCE_INLINE void AVX2FF16Butterfly2PreparedImpl(')
    text = once(text,'        if (!Inverse)\n            AVX2FF16AdjacentProductAdd',
        '        if (!Inverse && Scheduled)\n            AVX2FF16AdjacentProductAdd')
    anchor = 'template<bool Inverse>\nstatic void AVX2FF16Butterfly2('
    wrapper = '''template<bool Inverse>
static LEO_FORCE_INLINE void AVX2FF16Butterfly2Prepared(
    void* x, void* y, uint16_t log, const __m256i low[4],
    const __m256i high[4], uint64_t bytes)
{
#if !defined(LEO2_GFNI_VARIANT) && !defined(LEO2_AVX512_VARIANT)
    const bool enabled = !Inverse && leo_adjacent_schedule_enabled;
#if defined(LEO_ADJACENT_TRACE) && LEO_ADJACENT_TRACE
    if (!Inverse) LeoAdjacentRecord(LeoAdjacentForward, enabled, bytes);
#endif
    if (enabled)
        AVX2FF16Butterfly2PreparedImpl<Inverse, true>(x, y, log, low, high, bytes);
    else
#endif
        AVX2FF16Butterfly2PreparedImpl<Inverse, false>(x, y, log, low, high, bytes);
}

'''
    text = once(text,anchor,wrapper+anchor)
    start = text.index('static void AVX2FF16IFFTButterfly2Xor(')
    end = text.index('\n#endif // LEO_HAS_FF16',start)
    body = text[start:end]
    body = once(body,'static void AVX2FF16IFFTButterfly2Xor(',
        'template<bool Scheduled>\nstatic LEO_FORCE_INLINE void AVX2FF16IFFTButterfly2XorImpl(')
    opening = '#if '+ORDINARY+' && (LEO_AVX2_ADJACENT_SCHEDULE & 2)\n'
    a = body.index(opening)
    b = body.index('#else\n',a)
    c = body.index('#endif\n',b)
    scheduled, original = body[a+len(opening):b],body[b+len('#else\n'):c]
    body = body[:a]+opening+'        if (Scheduled)\n        {\n'+scheduled+\
        '        }\n        else\n#endif\n        {\n'+original+'        }\n'+body[c+len('#endif\n'):]
    opening = '#if defined(LEO2_GFNI_VARIANT) || defined(LEO2_AVX512_VARIANT) || !(LEO_AVX2_ADJACENT_SCHEDULE & 2)\n'
    a = body.index(opening); b = body.index('#endif\n',a)
    stores = body[a+len(opening):b]
    body = body[:a]+'#if '+ORDINARY+' && (LEO_AVX2_ADJACENT_SCHEDULE & 2)\n'+\
        '        if (!Scheduled)\n#endif\n        {\n'+stores+'        }\n'+body[b+len('#endif\n'):]
    wrapper = '''
static void AVX2FF16IFFTButterfly2Xor(
    const void* x, const void* y, void* u, void* v, uint16_t log, uint64_t bytes)
{
#if !defined(LEO2_GFNI_VARIANT) && !defined(LEO2_AVX512_VARIANT)
    const bool enabled = leo_adjacent_schedule_enabled;
#if defined(LEO_ADJACENT_TRACE) && LEO_ADJACENT_TRACE
    LeoAdjacentRecord(LeoAdjacentAccumulating, enabled, bytes);
#endif
    if (enabled)
        AVX2FF16IFFTButterfly2XorImpl<true>(x, y, u, v, log, bytes);
    else
#endif
        AVX2FF16IFFTButterfly2XorImpl<false>(x, y, u, v, log, bytes);
}
'''
    return text[:start]+body+wrapper+text[end:]

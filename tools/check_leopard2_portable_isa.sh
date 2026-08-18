#!/bin/sh

set -eu

usage()
{
    echo "usage: $0 OBJDUMP STATIC_ARCHIVE [BUILD_DIRECTORY [EXPECTED_CLASSES] | [BUILD_DIRECTORY CC AR [EXPECTED_CLASSES [optional|required]]]]" >&2
    echo "       $0 --self-test-only OBJDUMP CC AR" >&2
    exit 2
}

self_test_only=0
if [ "${1:-}" = --self-test-only ]; then
    [ "$#" -eq 4 ] || usage
    self_test_only=1
    objdump_bin=$2
    cc_bin=$3
    ar_bin=$4
    archive=
    build_dir=
    expected_classes=
    metadata_mode=optional
else
    case "$#" in
        2|3|4|5|6|7) ;;
        *) usage ;;
    esac

    objdump_bin=$1
    archive=$2
    build_dir=${3:-}
    cc_bin=
    ar_bin=
    expected_classes=
    metadata_mode=optional
fi
if [ "$self_test_only" -eq 0 ]; then
    case "$#" in
        4)
            expected_classes=$4
            ;;
        5)
            cc_bin=$4
            ar_bin=$5
            ;;
        6)
            cc_bin=$4
            ar_bin=$5
            expected_classes=$6
            ;;
        7)
            cc_bin=$4
            ar_bin=$5
            expected_classes=$6
            metadata_mode=$7
            ;;
    esac
fi
case "$metadata_mode" in
    optional|required) ;;
    *) usage ;;
esac

if [ "$self_test_only" -eq 0 ] && [ ! -f "$archive" ]; then
    echo "portable ISA check: archive not found: $archive" >&2
    exit 2
fi

# This is intentionally an x86-64 post-baseline denylist boundary: SSE2
# instructions are allowed and known compiler-generated extension families are
# rejected. Matching is performed only on objdump's mnemonic field, so symbol
# names and operands cannot trigger a false positive.
forbidden_mnemonics='^(addsubp[ds]|haddp[ds]|hsubp[ds]|lddqu|movddup|movshdup|movsldup|fisttp[[:alnum:]]*|monitor|mwait|monitorx|mwaitx|prefetchw|prefetchwt1|lahf|sahf|pabs[bdw]|palignr|phadd(d|sw|w)|phsub(d|sw|w)|pmaddubsw|pmulhrsw|pshufb|psign[bdw]|blendp[ds]|blendvp[ds]|dpp[ds]|extractps|insertps|movntdqa|mpsadbw|packusdw|pblendvb|pblendw|pcmpeqq|pextr[bdq]|phminposuw|pinsr[bdq]|pmax(sb|sd|uw|ud)|pmin(sb|sd|uw|ud)|pmov(sx|zx)(bd|bq|bw|dq|dw|wd|wq)|pmuldq|pmulld|ptest|round(p[ds]|s[ds])|crc32[[:alnum:]]*|pcmp(e|i)str[im]|pcmpgtq|popcnt|extrq|insertq|movntsd|movntss|v[[:alnum:]_.]*|cmpxchg16b|lzcnt|tzcnt|andn|bextr|blcfill|blci|blcic|blcmsk|blcs|blsfill|blsic|t1mskc|tzmsk|blsi|blsmsk|blsr|bzhi|pdep|pext|mulx|rorx|sarx|shlx|shrx|adcx|adox|movbe|rdrand|rdseed|rdpid|rdtscp|clflushopt|clwb|clzero|cldemote|wbnoinvd|serialize|umonitor|umwait|tpause|movdiri|movdir64b|enqcmd|enqcmds|xgetbv|xsetbv|xsave[[:alnum:]]*|xrstor[[:alnum:]]*|aes(enc|enclast|dec|declast|imc|keygenassist)|pclmul[[:alnum:]]*|gf2p8[[:alnum:]]*|sha1(msg1|msg2|nexte|rnds4)|sha256(msg1|msg2|rnds2)|k(add|and|andn|mov|not|or|ortest|shiftl|shiftr|test|unpck|xnor|xor)[bdqw]*)$'

# The AVX2 runtime probe establishes AVX, OS-managed XMM/YMM state, and AVX2
# only.  It does not establish FMA, F16C, XOP, AES, GFNI, VNNI, or any AVX-512
# subset.  Keep this list intentionally exact and fail closed: when an AVX2
# kernel starts using another AVX/AVX2 mnemonic, reviewers must add that
# mnemonic after checking its architectural feature contract.  A broad `v*'
# exemption would silently admit instructions whose CPUID bits are not probed.
allowed_avx2_vex_mnemonics='^(vbroadcastf128|vbroadcasti128|vbroadcastss|vextracti128|vinserti128|vmovaps|vmovd|vmovdqa|vmovdqu|vmovq|vmovups|vpackuswb|vpaddq|vpaddw|vpand|vpandn|vpblendd|vpblendw|vpbroadcastb|vpbroadcastd|vpbroadcastq|vpbroadcastw|vpcmpeqb|vpcmpeqd|vpcmpgtw|vperm2i128|vpermq|vpextrq|vpextrw|vpinsrb|vpminub|vpmovsxbw|vpmovzxbw|vpmuludq|vpmullw|vpshufb|vpshufd|vpshufhw|vpshuflw|vpsrldq|vpsrlq|vpsrlw|vpsubw|vpunpckldq|vpunpcklqdq|vpunpcklwd|vpxor|vxorps|vzeroupper)$'

# The GFNI candidate is compiled with `-mavx2 -mgfni -mno-avx512f' and is gated
# on a runtime probe that establishes AVX2 *and* the separately enumerated GFNI
# bit.  It is therefore exactly the reviewed AVX2 VEX set plus the VEX-encoded
# affine transform; derive it from that set so the two lists cannot drift.
# VGF2P8MULB and VGF2P8AFFINEINVQB are deliberately absent: the backend sources
# emit no such form, and admitting an unused mnemonic would widen the reviewed
# contract without evidence.  The class stays VEX-only, so EVEX-encoded GFNI is
# rejected below even though its mnemonic is spelled identically.
# Beyond the reviewed AVX2 set, the GFNI translation unit alone compiles the
# affine-multiply helpers and fused radix-eight kernels, which additionally
# emit these plain AVX2 integer mnemonics.  All are baseline AVX2: no new ISA
# floor beyond the member's AVX2+GFNI runtime probe.
allowed_gfni_extra_mnemonics='^(vpackusdw|vpaddb|vpextrb|vpinsrd|vpinsrw|vpmovsxdq|vpor|vpslld|vpsllq|vpsllw|vpsrad|vpsrld|vpsrlvq|vpsrlw|vpunpcklbw)$'
allowed_gfni_vex_mnemonics="$allowed_avx2_vex_mnemonics|$allowed_gfni_extra_mnemonics|^vgf2p8affineqb\$"

# The AVX-512 candidate is separately gated on F/BW/VL and complete OS ZMM
# state.  It deliberately retains 256-bit data width and uses EVEX only to
# access the expanded vector register file.  Keep the allowlist exact.
allowed_avx512vl_mnemonics='^(valignq|vbroadcasti32x4|vextracti64x2|vmovdqa32|vmovdqa64|vmovdqu8|vmovdqu64|vpandq|vpinsrq|vpshufb|vpsrlq|vpternlogd|vpternlogq|vpxord|vpxorq)$'

# The exact T=128 member is separately gated on AVX2, AVX-512F/BW/VL,
# complete OS ZMM state, and GFNI.  Unlike the general AVX-512 member it is
# intentionally ZMM-wide, and unlike the GFNI backend it is intentionally
# EVEX-encoded.  Keep its compiler-emitted initialization forms explicit;
# every listed instruction is covered by F/BW/VL plus baseline AVX2.
allowed_avx512_gfni_t128_extra_mnemonics='^(vgf2p8affineqb|vinserti64x4|vpaddb|vpbroadcastb|vpbroadcastq|vpermi2d|vpermi2w|vpermt2d|vpermt2w|vpmovwb|vpsllq|vpsllw|vpsrlw|vpsubb)$'

# Reject target-raising options in Make, Ninja, or compilation-database
# metadata.  -mno-* and the x86-64 SSE2 baseline remain allowed.  All -march
# and -mcpu values are conservatively rejected because proving that an
# arbitrary spelling retains the baseline is toolchain-specific.
forbidden_flags='(^|[[:space:]"=,:])(-march=[^[:space:]",]+|-mcpu=[^[:space:]",]+|-m(sse3|ssse3|sse4|sse4a|sse4[.][12]|avx[0-9.]*|avx512[^[:space:]",]*|fma|fma4|f16c|bmi|bmi2|lzcnt|popcnt|cx16|lahf-sahf|prfchw|mwaitx|clzero|wbnoinvd|aes|pclmul|vpclmulqdq|gfni|vaes|sha|adx|movbe|xsave[^[:space:]",]*|rdrnd|rdseed|xop|tbm|3dnow)([[:space:]",]|$)|/arch:(avx|avx2|avx512[^[:space:]",]*|sse3|ssse3|sse4([.][12])?)([[:space:]",]|$)|[+](sse3|ssse3|sse4a?|sse4[.]?[12]|avx[0-9.]*|avx512[^[:space:]",]*|fma|f16c|bmi2?|lzcnt|popcnt|cx16|lahf-sahf|prfchw|mwaitx|clzero|wbnoinvd|aes|pclmul|gfni|vaes|sha|adx|movbe)([[:space:]",]|$))'
forbidden_metadata="$forbidden_flags|(^|[[:space:]\"=,:])-flto([^[:space:]\",]*)([[:space:]\",]|$)"

scratch_root=$(mktemp -d "${TMPDIR:-/tmp}/leopard2-portable-isa.XXXXXX")
trap 'rm -rf "$scratch_root"' EXIT HUP INT TERM
scan_index=0

scan_object()
{
    object_file=$1
    object_class=$2
    object_name=$3
    scan_index=$((scan_index + 1))
    disassembly_file="$scratch_root/disassembly.$scan_index"
    raw_disassembly_file="$scratch_root/disassembly-raw.$scan_index"
    mnemonics_file="$scratch_root/mnemonics.$scan_index"
    forbidden_file="$scratch_root/forbidden.$scan_index"
    violations_file="$scratch_root/violations.$scan_index"

    LC_ALL=C "$objdump_bin" -d --no-show-raw-insn "$object_file" > "$disassembly_file"
    LC_ALL=C "$objdump_bin" -d "$object_file" > "$raw_disassembly_file"
    sed -nE 's/^[[:space:]]*[0-9a-f]+:[[:space:]]+([[:alnum:]_.]+).*/\1/p' \
        "$disassembly_file" | sort -u > "$mnemonics_file"

    if [ ! -s "$mnemonics_file" ]; then
        echo "portable ISA check: objdump produced no instruction mnemonics" >&2
        return 1
    fi

    grep -E "$forbidden_mnemonics" "$mnemonics_file" > "$forbidden_file" || true
    case "$object_class" in
        baseline)
            cp "$forbidden_file" "$violations_file"
            ;;
        cpu_features)
            grep -Ev '^xgetbv$' "$forbidden_file" > "$violations_file" || true
            ;;
        ssse3)
            # The named SSSE3 member may contain SSE3/SSSE3 instructions, but
            # no AVX, SSE4, feature-probe, or unrelated later extension.
            grep -Ev '^(addsubp[ds]|haddp[ds]|hsubp[ds]|lddqu|movddup|movshdup|movsldup|fisttp[[:alnum:]]*|pabs[bdw]|palignr|phadd(d|sw|w)|phsub(d|sw|w)|pmaddubsw|pmulhrsw|pshufb|psign[bdw])$' \
                "$forbidden_file" > "$violations_file" || true
            if ! grep -Eq '^pshufb$' "$mnemonics_file"; then
                echo "portable ISA check: SSSE3 member has no pshufb: $object_name" >&2
                return 1
            fi
            ;;
        avx2)
            # Admit only the reviewed AVX/AVX2 VEX mnemonics above plus the
            # SSSE3 family.  In particular, AVX2 does not imply FMA/F16C/XOP or
            # other separately probed extensions even though their mnemonics
            # also begin with `v'.
            grep -Ev "$allowed_avx2_vex_mnemonics|^(addsubp[ds]|haddp[ds]|hsubp[ds]|lddqu|movddup|movshdup|movsldup|fisttp[[:alnum:]]*|pabs[bdw]|palignr|phadd(d|sw|w)|phsub(d|sw|w)|pmaddubsw|pmulhrsw|pshufb|psign[bdw])$" \
                "$forbidden_file" > "$violations_file" || true
            if ! grep -Eq '^v[[:alnum:]_.]*$' "$mnemonics_file"; then
                echo "portable ISA check: AVX2 member has no VEX instruction: $object_name" >&2
                return 1
            fi
            # EVEX can encode AVX-512-family instructions using only XMM/YMM
            # operands and no mask, so operand spelling alone is insufficient.
            # Prefix byte 0x62 is EVEX in the x86-64 objects checked here.
            if grep -Eq '^[[:space:]]*[0-9a-f]+:[[:space:]]+62([[:space:]]|$)' \
                "$raw_disassembly_file"
            then
                echo "portable ISA check: AVX2 member contains an EVEX instruction: $object_name" >&2
                return 1
            fi
            if grep -Eq '(%zmm[0-9]+|%k[0-7]|\{1to[0-9]+\}|\{z\})' "$disassembly_file"; then
                echo "portable ISA check: AVX2 member contains AVX-512 operands: $object_name" >&2
                return 1
            fi
            ;;
        gfni)
            # The named GFNI member is a strict VEX-only superset of the AVX2
            # class: the reviewed AVX/AVX2 VEX mnemonics plus VGF2P8AFFINEQB,
            # plus the SSSE3 family.  Every other separately probed extension
            # still fails closed, including the remaining GFNI mnemonics.
            grep -Ev "$allowed_gfni_vex_mnemonics|^(addsubp[ds]|haddp[ds]|hsubp[ds]|lddqu|movddup|movshdup|movsldup|fisttp[[:alnum:]]*|pabs[bdw]|palignr|phadd(d|sw|w)|phsub(d|sw|w)|pmaddubsw|pmulhrsw|pshufb|psign[bdw])$" \
                "$forbidden_file" > "$violations_file" || true
            if ! grep -Eq '^v[[:alnum:]_.]*$' "$mnemonics_file"; then
                echo "portable ISA check: GFNI member has no VEX instruction: $object_name" >&2
                return 1
            fi
            # GFNI is also enumerated as an EVEX-encoded AVX-512 form whose
            # mnemonic is spelled identically, and EVEX can address only
            # XMM/YMM with no mask.  Reject on the 0x62 prefix byte so the
            # allowlist above cannot be used to smuggle in an AVX-512 encoding.
            if grep -Eq '^[[:space:]]*[0-9a-f]+:[[:space:]]+62([[:space:]]|$)' \
                "$raw_disassembly_file"
            then
                echo "portable ISA check: GFNI member contains an EVEX instruction: $object_name" >&2
                return 1
            fi
            if grep -Eq '(%zmm[0-9]+|%k[0-7]|\{1to[0-9]+\}|\{z\})' "$disassembly_file"; then
                echo "portable ISA check: GFNI member contains AVX-512 operands: $object_name" >&2
                return 1
            fi
            ;;
        avx512)
            grep -Ev "$allowed_avx2_vex_mnemonics|$allowed_avx512vl_mnemonics|^(addsubp[ds]|haddp[ds]|hsubp[ds]|lddqu|movddup|movshdup|movsldup|fisttp[[:alnum:]]*|pabs[bdw]|palignr|phadd(d|sw|w)|phsub(d|sw|w)|pmaddubsw|pmulhrsw|pshufb|psign[bdw])$" \
                "$forbidden_file" > "$violations_file" || true
            if ! grep -Eq '^[[:space:]]*[0-9a-f]+:[[:space:]]+62([[:space:]]|$)' \
                "$raw_disassembly_file"
            then
                echo "portable ISA check: AVX-512 member has no EVEX instruction: $object_name" >&2
                return 1
            fi
            if grep -Eq '(%zmm[0-9]+|\{1to[0-9]+\}|\{z\})' "$disassembly_file"; then
                echo "portable ISA check: AVX-512VL member widened beyond YMM: $object_name" >&2
                return 1
            fi
            ;;
        avx512_gfni_t128)
            grep -Ev "$allowed_avx2_vex_mnemonics|$allowed_avx512vl_mnemonics|$allowed_avx512_gfni_t128_extra_mnemonics|^(addsubp[ds]|haddp[ds]|hsubp[ds]|lddqu|movddup|movshdup|movsldup|fisttp[[:alnum:]]*|pabs[bdw]|palignr|phadd(d|sw|w)|phsub(d|sw|w)|pmaddubsw|pmulhrsw|pshufb|psign[bdw])$" \
                "$forbidden_file" > "$violations_file" || true
            if ! grep -Eq '^vgf2p8affineqb$' "$mnemonics_file"; then
                echo "portable ISA check: AVX-512/GFNI T128 member has no affine instruction: $object_name" >&2
                return 1
            fi
            if ! grep -Eq '^[[:space:]]*[0-9a-f]+:[[:space:]]+62([[:space:]]|$)' \
                "$raw_disassembly_file" ||
               ! grep -Eq '%zmm[0-9]+' "$disassembly_file"
            then
                echo "portable ISA check: AVX-512/GFNI T128 member is not ZMM/EVEX: $object_name" >&2
                return 1
            fi
            ;;
        *)
            echo "portable ISA check: unknown object class: $object_class" >&2
            return 1
            ;;
    esac

    if [ -s "$violations_file" ]; then
        cat "$violations_file" >&2
        echo "portable ISA check: forbidden instruction in $object_class member $object_name" >&2
        return 1
    fi
    return 0
}

require_expected_members()
{
    members_file=$1
    expected=$2
    old_ifs=$IFS
    IFS=,
    set -- $expected
    IFS=$old_ifs
    for object_class in "$@"
    do
        case "$object_class" in
            cpu_features) expected_member=Leopard2CpuFeatures.cpp.o ;;
            ssse3) expected_member=Leopard2BackendSSSE3.cpp.o ;;
            avx2) expected_member=Leopard2BackendAVX2.cpp.o ;;
            avx2_t16_q2) expected_member=Leopard2BackendAVX2T16Q2.cpp.o ;;
            avx2_t16_k66) expected_member=Leopard2BackendAVX2T16K66.cpp.o ;;
            avx2_t8_k62) expected_member=Leopard2BackendAVX2T8K62.cpp.o ;;
            avx2_xor) expected_member=Leopard2BackendAVX2Xor.cpp.o ;;
            avx2_t2_k4) expected_member=Leopard2BackendAVX2T2K4.cpp.o ;;
            avx2_t8_k8_b1024) expected_member=Leopard2BackendAVX2T8K8B1024.cpp.o ;;
            avx2_t16_b64) expected_member=Leopard2BackendAVX2T16B64.cpp.o ;;
            avx2_t32_b256) expected_member=Leopard2BackendAVX2T32B256.cpp.o ;;
            avx2_p32) expected_member=Leopard2LowP32B64AVX2.cpp.o ;;
            gfni) expected_member=Leopard2BackendGFNI.cpp.o ;;
            avx512) expected_member=Leopard2BackendAVX512.cpp.o ;;
            avx512_gfni_t128)
                expected_member=Leopard2BackendAVX512GFNIT128.cpp.o ;;
            '') continue ;;
            *)
                echo "portable ISA check: unknown expected class: $object_class" >&2
                return 1
                ;;
        esac
        member_count=$(awk -v expected="$expected_member" '
            BEGIN {
                expected_obj = expected
                sub(/[.]o$/, ".obj", expected_obj)
            }
            $0 == expected || $0 == expected_obj {
                ++count
            }
            END { print count + 0 }
        ' "$members_file")
        if [ "$member_count" -ne 1 ]; then
            echo "portable ISA check: expected exactly one $object_class member, found $member_count" >&2
            return 1
        fi
    done
    return 0
}

scan_archive()
{
    archive_file=$1
    archive_ar=${2:-${ar_bin:-}}
    expected=${3:-}
    if [ -z "$archive_ar" ]; then
        archive_ar=$(command -v ar || true)
    fi
    if [ -z "$archive_ar" ]; then
        echo "portable ISA check: ar is required for member classification" >&2
        return 1
    fi

    members_file="$scratch_root/members.$scan_index"
    scan_index=$((scan_index + 1))
    if ! "$archive_ar" t "$archive_file" > "$members_file"; then
        echo "portable ISA check: failed to list archive members: $archive_file" >&2
        return 1
    fi
    if [ ! -s "$members_file" ]; then
        echo "portable ISA check: archive contains no members: $archive_file" >&2
        return 1
    fi
    duplicate_members="$scratch_root/duplicate-members.$scan_index"
    LC_ALL=C sort "$members_file" | uniq -d > "$duplicate_members"
    if [ -s "$duplicate_members" ]; then
        cat "$duplicate_members" >&2
        echo "portable ISA check: duplicate archive member name" >&2
        return 1
    fi
    if [ -n "$expected" ] && ! require_expected_members "$members_file" "$expected"; then
        return 1
    fi

    member_index=0
    while IFS= read -r member
    do
        member_index=$((member_index + 1))
        extracted="$scratch_root/member.$member_index.o"
        if ! "$archive_ar" p "$archive_file" "$member" > "$extracted"; then
            echo "portable ISA check: failed to extract archive member: $member" >&2
            return 1
        fi
        case "$member" in
            Leopard2BackendSSSE3.cpp.o|Leopard2BackendSSSE3.cpp.obj)
                object_class=ssse3 ;;
            Leopard2BackendAVX2.cpp.o|Leopard2BackendAVX2.cpp.obj|\
            Leopard2BackendAVX2T16Q2.cpp.o|Leopard2BackendAVX2T16Q2.cpp.obj|\
            Leopard2BackendAVX2T16K66.cpp.o|Leopard2BackendAVX2T16K66.cpp.obj|\
            Leopard2BackendAVX2T8K62.cpp.o|Leopard2BackendAVX2T8K62.cpp.obj|\
            Leopard2BackendAVX2Xor.cpp.o|Leopard2BackendAVX2Xor.cpp.obj|\
            Leopard2BackendAVX2T2K4.cpp.o|Leopard2BackendAVX2T2K4.cpp.obj|\
            Leopard2BackendAVX2T8K8B1024.cpp.o|Leopard2BackendAVX2T8K8B1024.cpp.obj|\
            Leopard2BackendAVX2T16B64.cpp.o|Leopard2BackendAVX2T16B64.cpp.obj|\
            Leopard2BackendAVX2T32B256.cpp.o|Leopard2BackendAVX2T32B256.cpp.obj|\
            Leopard2LowP32B64AVX2.cpp.o|Leopard2LowP32B64AVX2.cpp.obj)
                object_class=avx2 ;;
            Leopard2BackendGFNI.cpp.o|Leopard2BackendGFNI.cpp.obj)
                object_class=gfni ;;
            Leopard2BackendAVX512.cpp.o|Leopard2BackendAVX512.cpp.obj)
                object_class=avx512 ;;
            Leopard2BackendAVX512GFNIT128.cpp.o|Leopard2BackendAVX512GFNIT128.cpp.obj)
                object_class=avx512_gfni_t128 ;;
            Leopard2CpuFeatures.cpp.o|Leopard2CpuFeatures.cpp.obj)
                object_class=cpu_features ;;
            *) object_class=baseline ;;
        esac
        scan_object "$extracted" "$object_class" "$member" || return 1
    done < "$members_file"
}

scan_build_metadata()
{
    metadata_root=$1
    required=${2:-optional}
    baseline_flags="$metadata_root/CMakeFiles/leopard.dir/flags.make"
    if [ -f "$baseline_flags" ] &&
       LC_ALL=C grep -Ein -- "$forbidden_metadata" "$baseline_flags"
    then
        echo "portable ISA check: baseline compiler metadata raises the ISA floor: $baseline_flags" >&2
        return 1
    fi

    compile_commands="$metadata_root/compile_commands.json"
    if [ "$required" = required ]; then
        if [ ! -s "$compile_commands" ]; then
            echo "portable ISA check: required compile metadata is missing or empty: $compile_commands" >&2
            return 1
        fi
        if ! grep -Eq '(^|[/\\])leopard[.]cpp' "$compile_commands"; then
            echo "portable ISA check: required compile metadata has no baseline leopard.cpp entry" >&2
            return 1
        fi
    fi
    if [ -f "$compile_commands" ]; then
        # A diagnostic-only target compiles the historical whole-translation-
        # unit SIMD path with SSSE3 so allocation-failure handling can be
        # exercised.  It is not linked into the installed archive.  Admit that
        # one fixture only when a unique CMake definition, the exact target,
        # the reviewed source/object mapping, and the exact SSSE3-only flag
        # contract all agree.  This is intentionally narrower than filtering
        # by a target-name substring: a lookalike target or a newly added
        # source must fail closed until it is reviewed here.
        is_legacy_simd_fixture_command()
        {
            command_line=$1
            case "$command_line" in
                *-DLEO2_PORTABLE_ISA_PRIVILEGED_FIXTURE=1*) ;;
                *) return 1 ;;
            esac

            mssse3_count=$(printf '%s\n' "$command_line" |
                grep -o -- '-mssse3' | wc -l | tr -d '[:space:]')
            mno_avx_count=$(printf '%s\n' "$command_line" |
                grep -o -- '-mno-avx' | wc -l | tr -d '[:space:]')
            if [ "$mssse3_count" -ne 1 ] || [ "$mno_avx_count" -ne 1 ]; then
                return 1
            fi

            fixture_without_ssse3=$(printf '%s\n' "$command_line" |
                sed 's/[[:space:]]-mssse3[[:space:]]/ /')
            if printf '%s\n' "$fixture_without_ssse3" |
               LC_ALL=C grep -Eq -- "$forbidden_metadata"
            then
                return 1
            fi

            fixture_target=leopard_legacy_simd_init_test_lib
            for fixture_source in \
                leopard.cpp \
                leopard2.cpp \
                Leopard2Backend.cpp \
                Leopard2BackendScalar.cpp \
                Leopard2CpuFeatures.cpp \
                Leopard2Plan.cpp \
                LeopardCommon.cpp \
                LeopardFF8.cpp \
                LeopardFF16.cpp
            do
                fixture_object="CMakeFiles/$fixture_target.dir/$fixture_source.o"
                case "$command_line" in
                    *" -o $fixture_object -c "*"/$fixture_source\""*)
                        return 0 ;;
                esac
            done

            fixture_source=tests/leopard2/test_legacy_simd_init_failure.cpp
            fixture_object="CMakeFiles/leopard2_legacy_simd_init_failure_test.dir/$fixture_source.o"
            case "$command_line" in
                *" -o $fixture_object -c "*"/$fixture_source\""*)
                    return 0 ;;
            esac
            return 1
        }

        # Every use of the privileged marker must itself satisfy the exact
        # fixture contract, even if a future edit accidentally drops the
        # otherwise-forbidden -mssse3 flag from that command.
        fixture_lines="$scratch_root/fixture-lines"
        LC_ALL=C grep -F -- \
            '-DLEO2_PORTABLE_ISA_PRIVILEGED_FIXTURE=1' \
            "$compile_commands" > "$fixture_lines" || true
        while IFS= read -r fixture_line
        do
            [ -n "$fixture_line" ] || continue
            if ! is_legacy_simd_fixture_command "$fixture_line"; then
                printf '%s\n' "$fixture_line" >&2
                echo "portable ISA check: invalid privileged diagnostic compile command" >&2
                return 1
            fi
        done < "$fixture_lines"

        violating_lines="$scratch_root/metadata-violations"
        candidate_lines="$scratch_root/metadata-candidates"
        LC_ALL=C grep -Ein -- "$forbidden_metadata" "$compile_commands" |
            grep -Ev '(^|[/[:space:]"])(Leopard2Backend(SSSE3|AVX2|AVX2T16Q2|AVX2T16K66|AVX2T8K62|AVX2Xor|AVX2T2K4|AVX2T8K8B1024|AVX2T16B64|AVX2T32B256|AVX512|AVX512GFNIT128|GFNI)|Leopard2LowP32B64AVX2)[.]cpp([[:space:]",]|$)' > \
                "$candidate_lines" || true
        : > "$violating_lines"
        while IFS= read -r candidate_line
        do
            [ -n "$candidate_line" ] || continue
            if ! is_legacy_simd_fixture_command "$candidate_line"; then
                printf '%s\n' "$candidate_line" >> "$violating_lines"
            fi
        done < "$candidate_lines"
        if [ -s "$violating_lines" ]; then
            cat "$violating_lines" >&2
            echo "portable ISA check: baseline compile command raises the ISA floor" >&2
            return 1
        fi
        if grep -q 'Leopard2BackendSSSE3[.]cpp' "$compile_commands" &&
           ! grep 'Leopard2BackendSSSE3[.]cpp' "$compile_commands" |
               grep -q -- '-mssse3'
        then
            echo "portable ISA check: SSSE3 object lacks its ISA flag" >&2
            return 1
        fi
        ssse3_command="$scratch_root/ssse3-command"
        grep 'Leopard2BackendSSSE3[.]cpp' "$compile_commands" |
            sed 's/-mssse3//g' > "$ssse3_command" || true
        if [ -s "$ssse3_command" ] &&
           grep -Ein -- "$forbidden_metadata" "$ssse3_command"
        then
            echo "portable ISA check: SSSE3 object has an unrelated ISA/LTO flag" >&2
            return 1
        fi
        for avx2_source in \
            Leopard2BackendAVX2.cpp \
            Leopard2BackendAVX2T16Q2.cpp \
            Leopard2BackendAVX2T16K66.cpp \
            Leopard2BackendAVX2T8K62.cpp \
            Leopard2BackendAVX2Xor.cpp \
            Leopard2LowP32B64AVX2.cpp
        do
            avx2_stem=${avx2_source%.cpp}
            avx2_pattern="(^|[/[:space:]\"])$avx2_stem[.]cpp([[:space:]\",]|$)"
            if grep -Eq "$avx2_pattern" "$compile_commands" &&
               ! grep -E "$avx2_pattern" "$compile_commands" |
                   grep -q -- '-mavx2'
            then
                echo "portable ISA check: $avx2_source lacks its ISA flag" >&2
                return 1
            fi
            avx2_command="$scratch_root/$avx2_source-command"
            grep -E "$avx2_pattern" "$compile_commands" |
                sed 's/-mavx2//g' > "$avx2_command" || true
            if [ -s "$avx2_command" ] &&
               grep -Ein -- "$forbidden_metadata" "$avx2_command"
            then
                echo "portable ISA check: $avx2_source has an unrelated ISA/LTO flag" >&2
                return 1
            fi
        done
        t2_k4_source=Leopard2BackendAVX2T2K4.cpp
        t2_k4_stem=${t2_k4_source%.cpp}
        t2_k4_source_pattern="(^|[/[:space:]\"])$t2_k4_stem[.]cpp([[:space:]\",]|$)"
        if grep -Eq "$t2_k4_source_pattern" "$compile_commands"; then
            t2_k4_lines="$scratch_root/t2-k4-lines"
            grep -E "$t2_k4_source_pattern" "$compile_commands" |
                grep -F ' -c ' > "$t2_k4_lines" || true
            if [ ! -s "$t2_k4_lines" ]; then
                echo "portable ISA check: $t2_k4_source has no compile command" >&2
                return 1
            fi
            t2_k4_object_pattern='(^|[/[:space:]"])CMakeFiles/leopard2_backend_avx2_t2_k4[.]dir/Leopard2BackendAVX2T2K4[.]cpp[.](o|obj)([[:space:]",]|$)'
            while IFS= read -r t2_k4_line
            do
                if ! printf '%s\n' "$t2_k4_line" |
                    grep -Eq "$t2_k4_object_pattern"
                then
                    echo "portable ISA check: invalid $t2_k4_source object mapping" >&2
                    return 1
                fi
                for required_flag in -mavx2 -mno-avx512f
                do
                    required_flag_pattern="(^|[[:space:]\"=,:])$required_flag([[:space:]\",]|$)"
                    if ! printf '%s\n' "$t2_k4_line" |
                        grep -Eq -- "$required_flag_pattern"
                    then
                        echo "portable ISA check: $t2_k4_source lacks $required_flag" >&2
                        return 1
                    fi
                done
                t2_k4_command="$scratch_root/t2-k4-command"
                printf '%s\n' "$t2_k4_line" |
                    sed 's/-mavx2//g' > "$t2_k4_command"
                if grep -Ein -- "$forbidden_metadata" "$t2_k4_command"; then
                    echo "portable ISA check: $t2_k4_source has an unrelated ISA/LTO flag" >&2
                    return 1
                fi
            done < "$t2_k4_lines"
        fi
        t8_k8_b1024_source=Leopard2BackendAVX2T8K8B1024.cpp
        t8_k8_b1024_stem=${t8_k8_b1024_source%.cpp}
        t8_k8_b1024_source_pattern="(^|[/[:space:]\"])$t8_k8_b1024_stem[.]cpp([[:space:]\",]|$)"
        if grep -Eq "$t8_k8_b1024_source_pattern" "$compile_commands"; then
            t8_k8_b1024_lines="$scratch_root/t8-k8-b1024-lines"
            grep -E "$t8_k8_b1024_source_pattern" "$compile_commands" |
                grep -F ' -c ' > "$t8_k8_b1024_lines" || true
            if [ ! -s "$t8_k8_b1024_lines" ]; then
                echo "portable ISA check: $t8_k8_b1024_source has no compile command" >&2
                return 1
            fi
            t8_k8_b1024_object_pattern='(^|[/[:space:]"])CMakeFiles/leopard2_backend_avx2_t8_k8_b1024[.]dir/Leopard2BackendAVX2T8K8B1024[.]cpp[.](o|obj)([[:space:]",]|$)'
            while IFS= read -r t8_k8_b1024_line
            do
                if ! printf '%s\n' "$t8_k8_b1024_line" |
                    grep -Eq "$t8_k8_b1024_object_pattern"
                then
                    echo "portable ISA check: invalid $t8_k8_b1024_source object mapping" >&2
                    return 1
                fi
                for required_flag in -mavx2 -mno-avx512f
                do
                    required_flag_pattern="(^|[[:space:]\"=,:])$required_flag([[:space:]\",]|$)"
                    if ! printf '%s\n' "$t8_k8_b1024_line" |
                        grep -Eq -- "$required_flag_pattern"
                    then
                        echo "portable ISA check: $t8_k8_b1024_source lacks $required_flag" >&2
                        return 1
                    fi
                done
                t8_k8_b1024_command="$scratch_root/t8-k8-b1024-command"
                printf '%s\n' "$t8_k8_b1024_line" |
                    sed 's/-mavx2//g' > "$t8_k8_b1024_command"
                if grep -Ein -- "$forbidden_metadata" "$t8_k8_b1024_command"; then
                    echo "portable ISA check: $t8_k8_b1024_source has an unrelated ISA/LTO flag" >&2
                    return 1
                fi
            done < "$t8_k8_b1024_lines"
        fi
        t16_source=Leopard2BackendAVX2T16B64.cpp
        t16_stem=${t16_source%.cpp}
        t16_source_pattern="(^|[/[:space:]\"])$t16_stem[.]cpp([[:space:]\",]|$)"
        if grep -Eq "$t16_source_pattern" "$compile_commands"; then
            t16_lines="$scratch_root/t16-b64-lines"
            grep -E "$t16_source_pattern" "$compile_commands" |
                grep -F ' -c ' > "$t16_lines" || true
            if [ ! -s "$t16_lines" ]; then
                echo "portable ISA check: $t16_source has no compile command" >&2
                return 1
            fi
            t16_object_pattern='(^|[/[:space:]"])CMakeFiles/leopard2_backend_avx2_t16_b64[.]dir/Leopard2BackendAVX2T16B64[.]cpp[.](o|obj)([[:space:]",]|$)'
            while IFS= read -r t16_line
            do
                if ! printf '%s\n' "$t16_line" |
                    grep -Eq "$t16_object_pattern"
                then
                    echo "portable ISA check: invalid $t16_source object mapping" >&2
                    return 1
                fi
                for required_flag in -mavx2 -mno-avx512f
                do
                    required_flag_pattern="(^|[[:space:]\"=,:])$required_flag([[:space:]\",]|$)"
                    if ! printf '%s\n' "$t16_line" |
                        grep -Eq -- "$required_flag_pattern"
                    then
                        echo "portable ISA check: $t16_source lacks $required_flag" >&2
                        return 1
                    fi
                done
                t16_command="$scratch_root/t16-b64-command"
                printf '%s\n' "$t16_line" |
                    sed 's/-mavx2//g' > "$t16_command"
                if grep -Ein -- "$forbidden_metadata" "$t16_command"; then
                    echo "portable ISA check: $t16_source has an unrelated ISA/LTO flag" >&2
                    return 1
                fi
            done < "$t16_lines"
        fi
        t16_k66_source=Leopard2BackendAVX2T16K66.cpp
        t16_k66_stem=${t16_k66_source%.cpp}
        t16_k66_source_pattern="(^|[/[:space:]\"])$t16_k66_stem[.]cpp([[:space:]\",]|$)"
        if grep -Eq "$t16_k66_source_pattern" "$compile_commands"; then
            t16_k66_lines="$scratch_root/t16-k66-lines"
            grep -E "$t16_k66_source_pattern" "$compile_commands" |
                grep -F ' -c ' > "$t16_k66_lines" || true
            if [ ! -s "$t16_k66_lines" ]; then
                echo "portable ISA check: $t16_k66_source has no compile command" >&2
                return 1
            fi
            t16_k66_object_pattern='(^|[/[:space:]"])CMakeFiles/leopard2_backend_avx2_t16_k66[.]dir/Leopard2BackendAVX2T16K66[.]cpp[.](o|obj)([[:space:]",]|$)'
            while IFS= read -r t16_k66_line
            do
                if ! printf '%s\n' "$t16_k66_line" |
                    grep -Eq "$t16_k66_object_pattern"
                then
                    echo "portable ISA check: invalid $t16_k66_source object mapping" >&2
                    return 1
                fi
                for required_flag in -mavx2 -mno-avx512f
                do
                    required_flag_pattern="(^|[[:space:]\"=,:])$required_flag([[:space:]\",]|$)"
                    if ! printf '%s\n' "$t16_k66_line" |
                        grep -Eq -- "$required_flag_pattern"
                    then
                        echo "portable ISA check: $t16_k66_source lacks $required_flag" >&2
                        return 1
                    fi
                done
                t16_k66_command="$scratch_root/t16-k66-command"
                printf '%s\n' "$t16_k66_line" |
                    sed 's/-mavx2//g' > "$t16_k66_command"
                if grep -Ein -- "$forbidden_metadata" "$t16_k66_command"; then
                    echo "portable ISA check: $t16_k66_source has an unrelated ISA/LTO flag" >&2
                    return 1
                fi
            done < "$t16_k66_lines"
        fi
        t32_source=Leopard2BackendAVX2T32B256.cpp
        t32_stem=${t32_source%.cpp}
        t32_source_pattern="(^|[/[:space:]\"])$t32_stem[.]cpp([[:space:]\",]|$)"
        if grep -Eq "$t32_source_pattern" "$compile_commands"; then
            t32_lines="$scratch_root/t32-b256-lines"
            grep -E "$t32_source_pattern" "$compile_commands" |
                grep -F ' -c ' > "$t32_lines" || true
            if [ ! -s "$t32_lines" ]; then
                echo "portable ISA check: $t32_source has no compile command" >&2
                return 1
            fi
            t32_object_pattern='(^|[/[:space:]"])CMakeFiles/leopard2_backend_avx2_t32_b256[.]dir/Leopard2BackendAVX2T32B256[.]cpp[.](o|obj)([[:space:]",]|$)'
            while IFS= read -r t32_line
            do
                if ! printf '%s\n' "$t32_line" |
                    grep -Eq "$t32_object_pattern"
                then
                    echo "portable ISA check: invalid $t32_source object mapping" >&2
                    return 1
                fi
                for required_flag in -mavx2 -mno-avx512f
                do
                    required_flag_pattern="(^|[[:space:]\"=,:])$required_flag([[:space:]\",]|$)"
                    if ! printf '%s\n' "$t32_line" |
                        grep -Eq -- "$required_flag_pattern"
                    then
                        echo "portable ISA check: $t32_source lacks $required_flag" >&2
                        return 1
                    fi
                done
                t32_command="$scratch_root/t32-b256-command"
                printf '%s\n' "$t32_line" | sed 's/-mavx2//g' > "$t32_command"
                if grep -Ein -- "$forbidden_metadata" "$t32_command"; then
                    echo "portable ISA check: $t32_source has an unrelated ISA/LTO flag" >&2
                    return 1
                fi
            done < "$t32_lines"
        fi
        if grep -q 'Leopard2BackendGFNI[.]cpp' "$compile_commands"; then
            gfni_lines="$scratch_root/gfni-lines"
            grep 'Leopard2BackendGFNI[.]cpp' "$compile_commands" > "$gfni_lines"
            # The GFNI translation unit is admitted above by name, so its exact
            # flag contract has to be re-proved here.  -mno-avx512f keeps the
            # unit VEX-only even if a future toolchain default enables EVEX.
            for required_flag in -mavx2 -mgfni -mno-avx512f
            do
                if ! grep -q -- "$required_flag" "$gfni_lines"; then
                    echo "portable ISA check: GFNI object lacks $required_flag" >&2
                    return 1
                fi
            done
            gfni_command="$scratch_root/gfni-command"
            # -mgfni is itself in the forbidden set, so strip exactly the
            # reviewed flags before re-checking the remainder.  -mno-avx512f is
            # a -mno-* spelling and is never forbidden, so it is left in place.
            sed 's/-mavx2//g; s/-mgfni//g' "$gfni_lines" > "$gfni_command"
            if grep -Ein -- "$forbidden_metadata" "$gfni_command"; then
                echo "portable ISA check: GFNI object has an unrelated ISA/LTO flag" >&2
                return 1
            fi
        fi
        if grep -q 'Leopard2BackendAVX512[.]cpp' "$compile_commands"; then
            avx512_lines="$scratch_root/avx512-lines"
            grep 'Leopard2BackendAVX512[.]cpp' "$compile_commands" > \
                "$avx512_lines"
            for required_flag in -mavx2 -mavx512f -mavx512bw -mavx512vl
            do
                if ! grep -q -- "$required_flag" "$avx512_lines"; then
                    echo "portable ISA check: AVX-512 object lacks $required_flag" >&2
                    return 1
                fi
            done
            avx512_command="$scratch_root/avx512-command"
            sed 's/-mavx2//g; s/-mavx512f//g; s/-mavx512bw//g; s/-mavx512vl//g' \
                "$avx512_lines" > "$avx512_command"
            if grep -Ein -- "$forbidden_metadata" "$avx512_command"; then
                echo "portable ISA check: AVX-512 object has an unrelated ISA/LTO flag" >&2
                return 1
            fi
        fi
        if grep -q 'Leopard2BackendAVX512GFNIT128[.]cpp' \
            "$compile_commands"
        then
            t128_lines="$scratch_root/avx512-gfni-t128-lines"
            grep 'Leopard2BackendAVX512GFNIT128[.]cpp' \
                "$compile_commands" > "$t128_lines"
            for required_flag in \
                -mavx2 -mavx512f -mavx512bw -mavx512vl -mgfni
            do
                required_flag_pattern="(^|[[:space:]\"=,:])$required_flag([[:space:]\",]|$)"
                if ! grep -Eq -- "$required_flag_pattern" "$t128_lines"; then
                    echo "portable ISA check: AVX-512/GFNI T128 object lacks $required_flag" >&2
                    return 1
                fi
            done
            t128_command="$scratch_root/avx512-gfni-t128-command"
            sed -E \
                's/(^|[[:space:]"=,:])-mavx2([[:space:]",]|$)/\1\2/g; s/(^|[[:space:]"=,:])-mavx512f([[:space:]",]|$)/\1\2/g; s/(^|[[:space:]"=,:])-mavx512bw([[:space:]",]|$)/\1\2/g; s/(^|[[:space:]"=,:])-mavx512vl([[:space:]",]|$)/\1\2/g; s/(^|[[:space:]"=,:])-mgfni([[:space:]",]|$)/\1\2/g' \
                "$t128_lines" > "$t128_command"
            if grep -Ein -- "$forbidden_metadata" "$t128_command"; then
                echo "portable ISA check: AVX-512/GFNI T128 object has an unrelated ISA/LTO flag" >&2
                return 1
            fi
        fi
    fi
    return 0
}

write_assembly_archive()
{
    fixture_name=$1
    fixture_instruction=$2
    fixture_source="$scratch_root/$fixture_name.s"
    fixture_object="$scratch_root/$fixture_name.o"
    fixture_archive="$scratch_root/lib$fixture_name.a"

    printf '%s\n' \
        '.text' \
        ".globl $fixture_name" \
        "$fixture_name:" \
        "    $fixture_instruction" \
        '    ret' > "$fixture_source"
    "$cc_bin" -c "$fixture_source" -o "$fixture_object"
    "$ar_bin" rcs "$fixture_archive" "$fixture_object"
    echo "$fixture_archive"
}

write_classified_archive()
{
    fixture_name=$1
    fixture_member=$2
    fixture_instruction=$3
    fixture_source="$scratch_root/$fixture_name.s"
    fixture_object="$scratch_root/$fixture_member"
    fixture_archive="$scratch_root/lib$fixture_name.a"

    printf '%s\n' \
        '.text' \
        ".globl $fixture_name" \
        "$fixture_name:" \
        "    $fixture_instruction" \
        '    ret' > "$fixture_source"
    "$cc_bin" -c "$fixture_source" -o "$fixture_object"
    "$ar_bin" rcs "$fixture_archive" "$fixture_object"
    echo "$fixture_archive"
}

write_duplicate_classified_archive()
{
    fixture_name=$1
    fixture_member=$2
    first_instruction=$3
    second_instruction=$4
    first_dir="$scratch_root/$fixture_name-first"
    second_dir="$scratch_root/$fixture_name-second"
    fixture_archive="$scratch_root/lib$fixture_name.a"
    mkdir -p "$first_dir" "$second_dir"

    printf '%s\n' \
        '.text' \
        ".globl ${fixture_name}_first" \
        "${fixture_name}_first:" \
        "    $first_instruction" \
        '    ret' > "$first_dir/$fixture_name.s"
    printf '%s\n' \
        '.text' \
        ".globl ${fixture_name}_second" \
        "${fixture_name}_second:" \
        "    $second_instruction" \
        '    ret' > "$second_dir/$fixture_name.s"
    "$cc_bin" -c "$first_dir/$fixture_name.s" -o "$first_dir/$fixture_member"
    "$cc_bin" -c "$second_dir/$fixture_name.s" -o "$second_dir/$fixture_member"
    (cd "$first_dir" && "$ar_bin" qc "$fixture_archive" "$fixture_member")
    (cd "$second_dir" && "$ar_bin" q "$fixture_archive" "$fixture_member")
    echo "$fixture_archive"
}

expect_classified_archive_accepted()
{
    fixture_name=$1
    fixture_member=$2
    fixture_instruction=$3
    fixture_archive=$(write_classified_archive \
        "$fixture_name" "$fixture_member" "$fixture_instruction")
    if ! scan_archive "$fixture_archive" > "$scratch_root/$fixture_name.log" 2>&1; then
        cat "$scratch_root/$fixture_name.log" >&2
        echo "portable ISA checker self-test: classified $fixture_name was rejected" >&2
        return 1
    fi
}

expect_classified_archive_rejected()
{
    fixture_name=$1
    fixture_member=$2
    fixture_instruction=$3
    fixture_archive=$(write_classified_archive \
        "$fixture_name" "$fixture_member" "$fixture_instruction")
    if scan_archive "$fixture_archive" > "$scratch_root/$fixture_name.log" 2>&1; then
        echo "portable ISA checker self-test: classified $fixture_name was accepted" >&2
        return 1
    fi
}

expect_archive_rejected()
{
    fixture_name=$1
    fixture_instruction=$2
    fixture_archive=$(write_assembly_archive "$fixture_name" "$fixture_instruction")
    fixture_log="$scratch_root/$fixture_name.log"
    if scan_archive "$fixture_archive" > "$fixture_log" 2>&1; then
        echo "portable ISA checker self-test: $fixture_name was not rejected" >&2
        return 1
    fi
    return 0
}

expect_duplicate_archive_rejected()
{
    fixture_name=$1
    fixture_member=$2
    first_instruction=$3
    second_instruction=$4
    fixture_archive=$(write_duplicate_classified_archive \
        "$fixture_name" "$fixture_member" "$first_instruction" "$second_instruction")
    fixture_log="$scratch_root/$fixture_name.log"
    if scan_archive "$fixture_archive" "$ar_bin" > "$fixture_log" 2>&1; then
        echo "portable ISA checker self-test: duplicate-member archive was accepted" >&2
        return 1
    fi
    if ! grep -q 'duplicate archive member name' "$fixture_log"; then
        cat "$fixture_log" >&2
        echo "portable ISA checker self-test: duplicate-member rejection reason missing" >&2
        return 1
    fi
}

expect_missing_expected_member_rejected()
{
    fixture_archive=$(write_classified_archive missing_expected_avx2 \
        'Leopard2BackendAVX2.cpp.o' 'vpxor %ymm0, %ymm0, %ymm0')
    fixture_log="$scratch_root/missing-expected-member.log"
    if scan_archive "$fixture_archive" "$ar_bin" \
        'cpu_features,avx2' > "$fixture_log" 2>&1
    then
        echo "portable ISA checker self-test: missing expected member was accepted" >&2
        return 1
    fi
    if ! grep -q 'expected exactly one cpu_features member, found 0' "$fixture_log"; then
        cat "$fixture_log" >&2
        echo "portable ISA checker self-test: missing-member rejection reason missing" >&2
        return 1
    fi
}

expect_missing_avx2_xor_member_rejected()
{
    fixture_archive=$(write_classified_archive missing_expected_avx2_xor \
        'Leopard2BackendAVX2.cpp.o' 'vpxor %ymm0, %ymm0, %ymm0')
    fixture_log="$scratch_root/missing-avx2-xor-member.log"
    if scan_archive "$fixture_archive" "$ar_bin" \
        'avx2,avx2_xor' > "$fixture_log" 2>&1
    then
        echo "portable ISA checker self-test: missing AVX2 XOR member was accepted" >&2
        return 1
    fi
    if ! grep -q 'expected exactly one avx2_xor member, found 0' \
        "$fixture_log"
    then
        cat "$fixture_log" >&2
        echo "portable ISA checker self-test: missing AVX2 XOR rejection reason missing" >&2
        return 1
    fi
}

expect_missing_avx2_t16_b64_member_rejected()
{
    fixture_archive=$(write_classified_archive missing_expected_avx2_t16_b64 \
        'Leopard2BackendAVX2.cpp.o' 'vpxor %ymm0, %ymm0, %ymm0')
    fixture_log="$scratch_root/missing-avx2-t16-b64-member.log"
    if scan_archive "$fixture_archive" "$ar_bin" \
        'avx2,avx2_t16_b64' > "$fixture_log" 2>&1
    then
        echo "portable ISA checker self-test: missing AVX2 T16/B64 member was accepted" >&2
        return 1
    fi
    if ! grep -q 'expected exactly one avx2_t16_b64 member, found 0' \
        "$fixture_log"
    then
        cat "$fixture_log" >&2
        echo "portable ISA checker self-test: missing AVX2 T16/B64 rejection reason missing" >&2
        return 1
    fi
}

expect_missing_avx2_t16_k66_member_rejected()
{
    fixture_archive=$(write_classified_archive missing_expected_avx2_t16_k66 \
        'Leopard2BackendAVX2.cpp.o' 'vpxor %ymm0, %ymm0, %ymm0')
    fixture_log="$scratch_root/missing-avx2-t16-k66-member.log"
    if scan_archive "$fixture_archive" "$ar_bin" \
        'avx2,avx2_t16_k66' > "$fixture_log" 2>&1
    then
        echo "portable ISA checker self-test: missing AVX2 T16/K66 member was accepted" >&2
        return 1
    fi
    if ! grep -q 'expected exactly one avx2_t16_k66 member, found 0' \
        "$fixture_log"
    then
        cat "$fixture_log" >&2
        echo "portable ISA checker self-test: missing AVX2 T16/K66 rejection reason missing" >&2
        return 1
    fi
}

expect_missing_avx2_t8_k8_b1024_member_rejected()
{
    fixture_archive=$(write_classified_archive missing_expected_avx2_t8_k8_b1024 \
        'Leopard2BackendAVX2.cpp.o' 'vpxor %ymm0, %ymm0, %ymm0')
    fixture_log="$scratch_root/missing-avx2-t8-k8-b1024-member.log"
    if scan_archive "$fixture_archive" "$ar_bin" \
        'avx2,avx2_t8_k8_b1024' > "$fixture_log" 2>&1
    then
        echo "portable ISA checker self-test: missing AVX2 T8/K8/B1024 member was accepted" >&2
        return 1
    fi
    if ! grep -q 'expected exactly one avx2_t8_k8_b1024 member, found 0' \
        "$fixture_log"
    then
        cat "$fixture_log" >&2
        echo "portable ISA checker self-test: missing AVX2 T8/K8/B1024 rejection reason missing" >&2
        return 1
    fi
}

expect_missing_avx2_t8_k62_member_rejected()
{
    fixture_archive=$(write_classified_archive missing_expected_avx2_t8_k62 \
        'Leopard2BackendAVX2.cpp.o' 'vpxor %ymm0, %ymm0, %ymm0')
    fixture_log="$scratch_root/missing-avx2-t8-k62-member.log"
    if scan_archive "$fixture_archive" "$ar_bin" \
        'avx2,avx2_t8_k62' > "$fixture_log" 2>&1
    then
        echo "portable ISA checker self-test: missing AVX2 T8/K62 member was accepted" >&2
        return 1
    fi
    if ! grep -q 'expected exactly one avx2_t8_k62 member, found 0' \
        "$fixture_log"
    then
        cat "$fixture_log" >&2
        echo "portable ISA checker self-test: missing AVX2 T8/K62 rejection reason missing" >&2
        return 1
    fi
}

expect_missing_avx2_t32_b256_member_rejected()
{
    fixture_archive=$(write_classified_archive missing_expected_avx2_t32_b256 \
        'Leopard2BackendAVX2.cpp.o' 'vpxor %ymm0, %ymm0, %ymm0')
    fixture_log="$scratch_root/missing-avx2-t32-b256-member.log"
    if scan_archive "$fixture_archive" "$ar_bin" \
        'avx2,avx2_t32_b256' > "$fixture_log" 2>&1
    then
        echo "portable ISA checker self-test: missing AVX2 T32/B256 member was accepted" >&2
        return 1
    fi
    if ! grep -q 'expected exactly one avx2_t32_b256 member, found 0' \
        "$fixture_log"
    then
        cat "$fixture_log" >&2
        echo "portable ISA checker self-test: missing AVX2 T32/B256 rejection reason missing" >&2
        return 1
    fi
}

expect_missing_avx2_t2_k4_member_rejected()
{
    fixture_archive=$(write_classified_archive missing_expected_avx2_t2_k4 \
        'Leopard2BackendAVX2.cpp.o' 'vpxor %ymm0, %ymm0, %ymm0')
    fixture_log="$scratch_root/missing-avx2-t2-k4-member.log"
    if scan_archive "$fixture_archive" "$ar_bin" \
        'avx2,avx2_t2_k4' > "$fixture_log" 2>&1
    then
        echo "portable ISA checker self-test: missing AVX2 T2/K4 member was accepted" >&2
        return 1
    fi
    if ! grep -q 'expected exactly one avx2_t2_k4 member, found 0' \
        "$fixture_log"
    then
        cat "$fixture_log" >&2
        echo "portable ISA checker self-test: missing AVX2 T2/K4 rejection reason missing" >&2
        return 1
    fi
}

expect_missing_gfni_member_rejected()
{
    fixture_archive=$(write_classified_archive missing_expected_gfni \
        'Leopard2BackendAVX2.cpp.o' 'vpxor %ymm0, %ymm0, %ymm0')
    fixture_log="$scratch_root/missing-gfni-member.log"
    if scan_archive "$fixture_archive" "$ar_bin" \
        'avx2,gfni' > "$fixture_log" 2>&1
    then
        echo "portable ISA checker self-test: missing GFNI member was accepted" >&2
        return 1
    fi
    if ! grep -q 'expected exactly one gfni member, found 0' \
        "$fixture_log"
    then
        cat "$fixture_log" >&2
        echo "portable ISA checker self-test: missing GFNI rejection reason missing" >&2
        return 1
    fi
}

expect_missing_avx512_gfni_t128_member_rejected()
{
    fixture_archive=$(write_classified_archive missing_expected_t128 \
        'Leopard2BackendAVX2.cpp.o' 'vpxor %ymm0, %ymm0, %ymm0')
    fixture_log="$scratch_root/missing-avx512-gfni-t128-member.log"
    if scan_archive "$fixture_archive" "$ar_bin" \
        'avx2,avx512_gfni_t128' > "$fixture_log" 2>&1
    then
        echo "portable ISA checker self-test: missing AVX-512/GFNI T128 member was accepted" >&2
        return 1
    fi
    if ! grep -q \
        'expected exactly one avx512_gfni_t128 member, found 0' \
        "$fixture_log"
    then
        cat "$fixture_log" >&2
        echo "portable ISA checker self-test: missing AVX-512/GFNI T128 rejection reason missing" >&2
        return 1
    fi
}

expect_listing_failure_rejected()
{
    fixture_archive=$(write_assembly_archive listing_failure_archive \
        'pxor %xmm0, %xmm0')
    failing_ar="$scratch_root/failing-ar.sh"
    printf '%s\n' \
        '#!/bin/sh' \
        'if [ "$1" = t ]; then exit 17; fi' \
        "exec \"$ar_bin\" \"\$@\"" > "$failing_ar"
    chmod +x "$failing_ar"
    fixture_log="$scratch_root/listing-failure.log"
    if scan_archive "$fixture_archive" "$failing_ar" > "$fixture_log" 2>&1; then
        echo "portable ISA checker self-test: archive listing failure was accepted" >&2
        return 1
    fi
    if ! grep -q 'failed to list archive members' "$fixture_log"; then
        cat "$fixture_log" >&2
        echo "portable ISA checker self-test: listing-failure rejection reason missing" >&2
        return 1
    fi
}

expect_metadata_rejected()
{
    fixture_name=$1
    metadata_kind=$2
    fixture_flag=$3
    fixture_dir="$scratch_root/metadata-$fixture_name"
    mkdir -p "$fixture_dir/CMakeFiles/leopard.dir"

    case "$metadata_kind" in
        make)
            printf 'CXX_FLAGS = -O2 %s\n' "$fixture_flag" > \
                "$fixture_dir/CMakeFiles/leopard.dir/flags.make"
            ;;
        compile_commands)
            printf '[{"command":"c++ -O2 %s -c x.cpp","file":"x.cpp"}]\n' \
                "$fixture_flag" > "$fixture_dir/compile_commands.json"
            ;;
        *)
            echo "portable ISA checker self-test: unknown metadata fixture" >&2
            return 1
            ;;
    esac

    fixture_log="$scratch_root/metadata-$fixture_name.log"
    if scan_build_metadata "$fixture_dir" > "$fixture_log" 2>&1; then
        echo "portable ISA checker self-test: $fixture_name flag was not rejected" >&2
        return 1
    fi
    return 0
}

expect_required_metadata_rejected()
{
    fixture_name=$1
    fixture_kind=$2
    fixture_dir="$scratch_root/required-metadata-$fixture_name"
    mkdir -p "$fixture_dir"
    case "$fixture_kind" in
        missing) ;;
        empty) : > "$fixture_dir/compile_commands.json" ;;
        unrelated)
            printf '[{"command":"c++ -O2 -c unrelated.cpp","file":"unrelated.cpp"}]\n' > \
                "$fixture_dir/compile_commands.json"
            ;;
        *)
            echo "portable ISA checker self-test: unknown required metadata fixture" >&2
            return 1
            ;;
    esac
    fixture_log="$scratch_root/required-metadata-$fixture_name.log"
    if scan_build_metadata "$fixture_dir" required > "$fixture_log" 2>&1; then
        echo "portable ISA checker self-test: $fixture_name metadata was accepted" >&2
        return 1
    fi
    if ! grep -q 'required compile metadata' "$fixture_log"; then
        cat "$fixture_log" >&2
        echo "portable ISA checker self-test: required-metadata rejection reason missing" >&2
        return 1
    fi
}

write_compile_command_fixture()
{
    fixture_dir=$1
    fixture_command=$2
    fixture_file=$3
    fixture_output=$4
    mkdir -p "$fixture_dir"
    printf '%s\n' \
        '[' \
        '  {' \
        '    "directory": "/tmp/leopard-build",' \
        "    \"command\": \"$fixture_command\"," \
        "    \"file\": \"$fixture_file\"," \
        "    \"output\": \"$fixture_output\"" \
        '  }' \
        ']' > "$fixture_dir/compile_commands.json"
}

expect_fixture_metadata_accepted()
{
    fixture_name=$1
    fixture_command=$2
    fixture_file=$3
    fixture_output=$4
    fixture_dir="$scratch_root/fixture-metadata-$fixture_name"
    write_compile_command_fixture "$fixture_dir" "$fixture_command" \
        "$fixture_file" "$fixture_output"
    fixture_log="$scratch_root/fixture-metadata-$fixture_name.log"
    if ! scan_build_metadata "$fixture_dir" > "$fixture_log" 2>&1; then
        cat "$fixture_log" >&2
        echo "portable ISA checker self-test: valid fixture $fixture_name was rejected" >&2
        return 1
    fi
}

expect_fixture_metadata_rejected()
{
    fixture_name=$1
    fixture_command=$2
    fixture_file=$3
    fixture_output=$4
    fixture_dir="$scratch_root/fixture-metadata-$fixture_name"
    write_compile_command_fixture "$fixture_dir" "$fixture_command" \
        "$fixture_file" "$fixture_output"
    fixture_log="$scratch_root/fixture-metadata-$fixture_name.log"
    if scan_build_metadata "$fixture_dir" > "$fixture_log" 2>&1; then
        echo "portable ISA checker self-test: invalid fixture $fixture_name was accepted" >&2
        return 1
    fi
}

run_negative_controls()
{
    baseline_archive=$(write_assembly_archive baseline_sse2 'pxor %xmm0, %xmm0')
    if ! scan_archive "$baseline_archive" > "$scratch_root/baseline.log" 2>&1; then
        echo "portable ISA checker self-test: SSE2 baseline was rejected" >&2
        return 1
    fi

    expect_archive_rejected bad_ssse3 'pshufb %xmm0, %xmm0'
    expect_archive_rejected bad_sse41 'pminud %xmm0, %xmm0'
    expect_archive_rejected bad_avx2 'vpaddd %ymm0, %ymm0, %ymm0'
    expect_archive_rejected bad_vextracti128 \
        'vextracti128 $1, %ymm0, %xmm0'
    expect_archive_rejected bad_vpextrq 'vpextrq $1, %xmm0, %rax'

    expect_classified_archive_accepted good_probe \
        'Leopard2CpuFeatures.cpp.o' 'xgetbv'
    expect_classified_archive_accepted good_ssse3 \
        'Leopard2BackendSSSE3.cpp.o' 'pshufb %xmm0, %xmm0'
    expect_classified_archive_accepted good_avx2 \
        'Leopard2BackendAVX2.cpp.o' 'vpxor %ymm0, %ymm0, %ymm0'
    expect_classified_archive_accepted good_avx2_xor \
        'Leopard2BackendAVX2Xor.cpp.o' 'vpxor %ymm0, %ymm0, %ymm0'
    expect_classified_archive_accepted good_avx2_t2_k4 \
        'Leopard2BackendAVX2T2K4.cpp.o' \
        'vpxor %ymm0, %ymm0, %ymm0'
    expect_classified_archive_accepted good_avx2_t8_k8_b1024 \
        'Leopard2BackendAVX2T8K8B1024.cpp.o' \
        'vpxor %ymm0, %ymm0, %ymm0'
    expect_classified_archive_accepted good_avx2_t8_k62 \
        'Leopard2BackendAVX2T8K62.cpp.o' \
        'vpxor %ymm0, %ymm0, %ymm0'
    expect_classified_archive_accepted good_avx2_t16_b64 \
        'Leopard2BackendAVX2T16B64.cpp.o' \
        'vpxor %ymm0, %ymm0, %ymm0'
    expect_classified_archive_accepted good_avx2_t16_k66 \
        'Leopard2BackendAVX2T16K66.cpp.o' \
        'vpxor %ymm0, %ymm0, %ymm0'
    expect_classified_archive_accepted good_avx2_t32_b256 \
        'Leopard2BackendAVX2T32B256.cpp.o' \
        'vpxor %ymm0, %ymm0, %ymm0'
    # GCC can vectorize the in-range tail pointer-map construction as packed
    # 64-bit additions.  VPADDQ is baseline AVX2, and the classified-object
    # scan below still rejects EVEX encodings and AVX-512 operands.
    expect_classified_archive_accepted good_avx2_pointer_add \
        'Leopard2BackendAVX2.cpp.o' 'vpaddq %ymm0, %ymm0, %ymm0'
    # GCC 13 uses VPEXTRW for an exact two-byte tiny tail and vectorizes the
    # bounded source-row index arithmetic with VPMULUDQ.  Both VEX forms are
    # covered by the AVX2 runtime gate; the generic baseline fixture above
    # still rejects them and the classified scan still rejects EVEX.
    expect_classified_archive_accepted good_avx2_gcc13_tails \
        'Leopard2BackendAVX2.cpp.o' \
        'vpextrw $0, %xmm0, (%rax); vpmuludq %ymm0, %ymm0, %ymm0'
    # The active GF8 Walsh locator performs packed 16-bit modulo-255
    # butterflies and byte/word expansion entirely within baseline AVX2.
    # Keep every newly admitted mnemonic in one classified positive control;
    # the generic baseline negative controls still reject all VEX code.
    expect_classified_archive_accepted good_avx2_walsh_locator \
        'Leopard2BackendAVX2.cpp.o' \
        'vpackuswb %ymm0, %ymm0, %ymm0; vpaddw %ymm0, %ymm0, %ymm0; vpandn %ymm0, %ymm0, %ymm0; vpblendd $0, %ymm0, %ymm0, %ymm0; vpblendw $0, %ymm0, %ymm0, %ymm0; vpbroadcastw %xmm0, %ymm0; vpcmpeqb %ymm0, %ymm0, %ymm0; vpcmpgtw %ymm0, %ymm0, %ymm0; vperm2i128 $0, %ymm0, %ymm0, %ymm0; vpermq $0, %ymm0, %ymm0; vpmovzxbw %xmm0, %ymm0; vpmullw %ymm0, %ymm0, %ymm0; vpshufd $0, %ymm0, %ymm0; vpshufhw $0, %ymm0, %ymm0; vpshuflw $0, %ymm0, %ymm0; vpsrlw $1, %ymm0, %ymm0; vpsubw %ymm0, %ymm0, %ymm0'
    # The GFNI member carries the affine transform alongside ordinary AVX2 VEX
    # code, so the fixture must prove both are admitted in the same object.
    expect_classified_archive_accepted good_gfni \
        'Leopard2BackendGFNI.cpp.o' \
        'vgf2p8affineqb $0, %ymm0, %ymm0, %ymm0; vpxor %ymm0, %ymm0, %ymm0'
    # Clang 18's GFNI translation unit uses these additional ordinary
    # VEX-encoded AVX/AVX2 operations around the affine primitive.  They do
    # not widen the runtime contract beyond the already-required AVX2 bit.
    expect_classified_archive_accepted good_clang18_gfni_vectorization \
        'Leopard2BackendGFNI.cpp.o' \
        'vpackusdw %ymm0, %ymm0, %ymm0; vpextrb $0, %xmm0, %eax; vpinsrd $0, %eax, %xmm0, %xmm0; vpinsrw $0, %eax, %xmm0, %xmm0; vpmovsxdq %xmm0, %ymm0; vpslld $1, %ymm0, %ymm0; vpsllq $1, %ymm0, %ymm0; vpsrad $1, %ymm0, %ymm0; vpsrld $1, %ymm0, %ymm0; vpsrlvq %ymm0, %ymm0, %ymm0; vpunpcklbw %ymm0, %ymm0, %ymm0; vgf2p8affineqb $0, %ymm0, %ymm0, %ymm0'
    expect_classified_archive_accepted good_avx512vl \
        'Leopard2BackendAVX512.cpp.o' \
        'vpternlogq $0, %ymm0, %ymm0, %ymm0'
    expect_classified_archive_accepted good_avx512vl_align \
        'Leopard2BackendAVX512.cpp.o' \
        'valignq $1, %ymm0, %ymm0, %ymm0'
    expect_classified_archive_accepted good_avx512vl_extract \
        'Leopard2BackendAVX512.cpp.o' \
        'vextracti64x2 $1, %ymm0, %xmm0'
    expect_classified_archive_accepted good_avx512_gfni_t128 \
        'Leopard2BackendAVX512GFNIT128.cpp.o' \
        'vgf2p8affineqb $0, %zmm0, %zmm0, %zmm0; vpternlogd $0, %zmm0, %zmm0, %zmm0'
    expect_classified_archive_rejected avx512_gfni_t128_leaks_fma \
        'Leopard2BackendAVX512GFNIT128.cpp.o' \
        'vfmadd132ps %zmm0, %zmm0, %zmm0'
    expect_classified_archive_rejected avx512_gfni_t128_leaks_mulb \
        'Leopard2BackendAVX512GFNIT128.cpp.o' \
        'vgf2p8mulb %zmm0, %zmm0, %zmm0'
    expect_classified_archive_rejected avx512_gfni_t128_stays_ymm \
        'Leopard2BackendAVX512GFNIT128.cpp.o' \
        'vgf2p8affineqb $0, %ymm0, %ymm0, %ymm0'
    expect_classified_archive_rejected avx512_gfni_t128_leaks_probe \
        'Leopard2BackendAVX512GFNIT128.cpp.o' 'xgetbv'
    # Sanitizer instrumentation can materialize an integer argument with the
    # VEX-encoded AVX vmovd form.  AVX is already required by the runtime
    # probe, so admit this exact move while retaining the fail-closed v* list.
    expect_classified_archive_accepted good_avx_move \
        'Leopard2BackendAVX2.cpp.o' 'vmovd %eax, %xmm0'
    # Clang spells a 128-bit integer-table broadcast as the AVX
    # VBROADCASTF128 encoding; the operation is bitwise and AVX is already
    # established together with OS-managed YMM state by the AVX2 probe.
    expect_classified_archive_accepted good_avx_broadcast \
        'Leopard2BackendAVX2.cpp.o' 'vbroadcastf128 (%rax), %ymm0'
    # Clang 18 selects these VEX-encoded AVX/AVX2 forms while vectorizing the
    # same fixed-table and tail-map code for which GCC uses the reviewed forms
    # above.  Each instruction remains inside the AVX2 probe contract.  The
    # raw-prefix and operand checks still reject EVEX, ZMM, and opmask forms.
    expect_classified_archive_accepted good_clang18_avx2_vectorization \
        'Leopard2BackendAVX2.cpp.o' \
        'vbroadcastss (%rax), %ymm0; vmovaps %ymm0, %ymm1; vpbroadcastd %xmm0, %ymm0; vpcmpeqd %ymm0, %ymm0, %ymm0; vpminub %ymm0, %ymm0, %ymm0; vpmovsxbw %xmm0, %ymm0'
    # GCC can use these AVX/AVX2 extraction forms while scalarizing a local
    # pointer vector at -O2.  The AVX2 runtime probe already establishes their
    # architectural and OS-state requirements.  Keep them class-scoped so a
    # baseline or SSSE3 object still fails closed.
    expect_classified_archive_accepted good_avx2_extract_lane \
        'Leopard2BackendAVX2.cpp.o' \
        'vextracti128 $1, %ymm0, %xmm0'
    expect_classified_archive_accepted good_avx_extract_qword \
        'Leopard2BackendAVX2.cpp.o' 'vpextrq $1, %xmm0, %rax'
    expect_classified_archive_rejected ssse3_leaks_vextracti128 \
        'Leopard2BackendSSSE3.cpp.o' \
        'vextracti128 $1, %ymm0, %xmm0'
    expect_classified_archive_rejected ssse3_leaks_vpextrq \
        'Leopard2BackendSSSE3.cpp.o' 'vpextrq $1, %xmm0, %rax'
    expect_classified_archive_rejected avx2_leaks_fma \
        'Leopard2BackendAVX2.cpp.o' \
        'vfmadd132ps %ymm0, %ymm0, %ymm0'
    expect_classified_archive_rejected avx2_xor_leaks_fma \
        'Leopard2BackendAVX2Xor.cpp.o' \
        'vfmadd132ps %ymm0, %ymm0, %ymm0'
    expect_classified_archive_rejected avx2_t2_k4_leaks_fma \
        'Leopard2BackendAVX2T2K4.cpp.o' \
        'vfmadd132ps %ymm0, %ymm0, %ymm0'
    expect_classified_archive_rejected avx2_t2_k4_leaks_evex \
        'Leopard2BackendAVX2T2K4.cpp.o' \
        'vpternlogd $0, %ymm0, %ymm0, %ymm0'
    expect_classified_archive_rejected avx2_t2_k4_leaks_zmm \
        'Leopard2BackendAVX2T2K4.cpp.o' \
        'vpxord %zmm0, %zmm0, %zmm0'
    expect_classified_archive_rejected avx2_t2_k4_leaks_probe \
        'Leopard2BackendAVX2T2K4.cpp.o' 'xgetbv'
    expect_classified_archive_rejected avx2_t8_k8_b1024_leaks_fma \
        'Leopard2BackendAVX2T8K8B1024.cpp.o' \
        'vfmadd132ps %ymm0, %ymm0, %ymm0'
    expect_classified_archive_rejected avx2_t8_k8_b1024_leaks_evex \
        'Leopard2BackendAVX2T8K8B1024.cpp.o' \
        'vpternlogd $0, %ymm0, %ymm0, %ymm0'
    expect_classified_archive_rejected avx2_t8_k8_b1024_leaks_zmm \
        'Leopard2BackendAVX2T8K8B1024.cpp.o' \
        'vpxord %zmm0, %zmm0, %zmm0'
    expect_classified_archive_rejected avx2_t8_k8_b1024_leaks_probe \
        'Leopard2BackendAVX2T8K8B1024.cpp.o' 'xgetbv'
    expect_classified_archive_rejected avx2_t8_k62_leaks_fma \
        'Leopard2BackendAVX2T8K62.cpp.o' \
        'vfmadd132ps %ymm0, %ymm0, %ymm0'
    expect_classified_archive_rejected avx2_t8_k62_leaks_evex \
        'Leopard2BackendAVX2T8K62.cpp.o' \
        'vpternlogd $0, %ymm0, %ymm0, %ymm0'
    expect_classified_archive_rejected avx2_t8_k62_leaks_zmm \
        'Leopard2BackendAVX2T8K62.cpp.o' \
        'vpxord %zmm0, %zmm0, %zmm0'
    expect_classified_archive_rejected avx2_t8_k62_leaks_probe \
        'Leopard2BackendAVX2T8K62.cpp.o' 'xgetbv'
    expect_classified_archive_rejected avx2_t16_b64_leaks_fma \
        'Leopard2BackendAVX2T16B64.cpp.o' \
        'vfmadd132ps %ymm0, %ymm0, %ymm0'
    expect_classified_archive_rejected avx2_t16_b64_leaks_evex \
        'Leopard2BackendAVX2T16B64.cpp.o' \
        'vpternlogd $0, %ymm0, %ymm0, %ymm0'
    expect_classified_archive_rejected avx2_t16_b64_leaks_zmm \
        'Leopard2BackendAVX2T16B64.cpp.o' \
        'vpxord %zmm0, %zmm0, %zmm0'
    expect_classified_archive_rejected avx2_t16_b64_leaks_probe \
        'Leopard2BackendAVX2T16B64.cpp.o' 'xgetbv'
    expect_classified_archive_rejected avx2_t16_k66_leaks_fma \
        'Leopard2BackendAVX2T16K66.cpp.o' \
        'vfmadd132ps %ymm0, %ymm0, %ymm0'
    expect_classified_archive_rejected avx2_t16_k66_leaks_evex \
        'Leopard2BackendAVX2T16K66.cpp.o' \
        'vpternlogd $0, %ymm0, %ymm0, %ymm0'
    expect_classified_archive_rejected avx2_t16_k66_leaks_zmm \
        'Leopard2BackendAVX2T16K66.cpp.o' \
        'vpxord %zmm0, %zmm0, %zmm0'
    expect_classified_archive_rejected avx2_t16_k66_leaks_probe \
        'Leopard2BackendAVX2T16K66.cpp.o' 'xgetbv'
    expect_classified_archive_rejected avx2_t32_b256_leaks_fma \
        'Leopard2BackendAVX2T32B256.cpp.o' \
        'vfmadd132ps %ymm0, %ymm0, %ymm0'
    expect_classified_archive_rejected avx2_t32_b256_leaks_evex \
        'Leopard2BackendAVX2T32B256.cpp.o' \
        'vpternlogd $0, %ymm0, %ymm0, %ymm0'
    expect_classified_archive_rejected avx2_t32_b256_leaks_zmm \
        'Leopard2BackendAVX2T32B256.cpp.o' \
        'vpxord %zmm0, %zmm0, %zmm0'
    expect_classified_archive_rejected avx2_t32_b256_leaks_probe \
        'Leopard2BackendAVX2T32B256.cpp.o' 'xgetbv'
    expect_classified_archive_rejected avx2_leaks_f16c \
        'Leopard2BackendAVX2.cpp.o' 'vcvtph2ps %xmm0, %ymm0'
    expect_classified_archive_rejected avx2_leaks_evex \
        'Leopard2BackendAVX2.cpp.o' \
        'vpternlogd $0, %ymm0, %ymm0, %ymm0'
    expect_classified_archive_rejected avx2_leaks_avx512_align \
        'Leopard2BackendAVX2.cpp.o' \
        'valignq $1, %ymm0, %ymm0, %ymm0'
    expect_classified_archive_rejected avx2_leaks_avx512_extract \
        'Leopard2BackendAVX2.cpp.o' \
        'vextracti64x2 $1, %ymm0, %xmm0'
    expect_classified_archive_rejected probe_leaks_avx2 \
        'Leopard2CpuFeatures.cpp.o' 'vpaddd %ymm0, %ymm0, %ymm0'
    expect_classified_archive_rejected ssse3_leaks_avx2 \
        'Leopard2BackendSSSE3.cpp.o' 'vpaddd %ymm0, %ymm0, %ymm0'
    expect_classified_archive_rejected avx2_leaks_probe \
        'Leopard2BackendAVX2.cpp.o' 'xgetbv'
    expect_classified_archive_rejected avx2_leaks_avx512 \
        'Leopard2BackendAVX2.cpp.o' 'vpxord %zmm0, %zmm0, %zmm0'
    expect_classified_archive_rejected avx512vl_leaks_zmm \
        'Leopard2BackendAVX512.cpp.o' \
        'vpternlogq $0, %zmm0, %zmm0, %zmm0'
    expect_classified_archive_rejected avx512vl_leaks_gfni \
        'Leopard2BackendAVX512.cpp.o' \
        'vgf2p8affineqb $0, %ymm0, %ymm0, %ymm0'
    # Admitting VGF2P8AFFINEQB in the GFNI class must not relax any other
    # class, and must not relax the GFNI class beyond VEX or beyond the one
    # reviewed GFNI mnemonic.
    expect_classified_archive_rejected avx2_still_rejects_gfni \
        'Leopard2BackendAVX2.cpp.o' \
        'vgf2p8affineqb $0, %ymm0, %ymm0, %ymm0'
    expect_classified_archive_rejected gfni_leaks_fma \
        'Leopard2BackendGFNI.cpp.o' \
        'vfmadd132ps %ymm0, %ymm0, %ymm0'
    expect_classified_archive_rejected gfni_leaks_gf2p8mulb \
        'Leopard2BackendGFNI.cpp.o' 'vgf2p8mulb %ymm0, %ymm0, %ymm0'
    expect_classified_archive_rejected gfni_leaks_evex \
        'Leopard2BackendGFNI.cpp.o' \
        'vpternlogd $0, %ymm0, %ymm0, %ymm0'
    # The EVEX encoding of the allowed mnemonic: same spelling, YMM operands,
    # but an AVX-512 encoding and an opmask the AVX2+GFNI probe never gates.
    expect_classified_archive_rejected gfni_leaks_evex_affine \
        'Leopard2BackendGFNI.cpp.o' \
        'vgf2p8affineqb $0, %ymm0, %ymm0, %ymm0{%k1}'
    # The hardest case, and the reason the 0x62 prefix check cannot be dropped
    # in favour of operand spelling: an allowed mnemonic, a legitimate VEX
    # instruction beside it, YMM-only operands and no mask, yet still an
    # EVEX encoding that requires AVX-512F the GFNI probe never establishes.
    expect_classified_archive_rejected gfni_leaks_evex_encoding_only \
        'Leopard2BackendGFNI.cpp.o' \
        'vpxor %ymm0, %ymm0, %ymm0; {evex} vgf2p8affineqb $0, %ymm0, %ymm0, %ymm0'
    expect_classified_archive_rejected gfni_leaks_zmm \
        'Leopard2BackendGFNI.cpp.o' \
        'vgf2p8affineqb $0, %zmm0, %zmm0, %zmm0'
    expect_classified_archive_rejected gfni_leaks_probe \
        'Leopard2BackendGFNI.cpp.o' 'xgetbv'
    expect_classified_archive_rejected lookalike_ssse3_member \
        'NotLeopard2BackendSSSE3.cpp.o' 'pshufb %xmm0, %xmm0'
    expect_classified_archive_rejected lookalike_avx2_member \
        'NotLeopard2BackendAVX2.cpp.o' 'vpaddd %ymm0, %ymm0, %ymm0'
    expect_classified_archive_rejected lookalike_avx2_t2_k4_prefix \
        'NotLeopard2BackendAVX2T2K4.cpp.o' \
        'vpxor %ymm0, %ymm0, %ymm0'
    expect_classified_archive_rejected lookalike_avx2_t2_k4_suffix \
        'Leopard2BackendAVX2T2K4Extra.cpp.o' \
        'vpxor %ymm0, %ymm0, %ymm0'
    expect_classified_archive_rejected lookalike_avx2_t8_k8_b1024_prefix \
        'NotLeopard2BackendAVX2T8K8B1024.cpp.o' \
        'vpxor %ymm0, %ymm0, %ymm0'
    expect_classified_archive_rejected lookalike_avx2_t8_k8_b1024_suffix \
        'Leopard2BackendAVX2T8K8B1024Extra.cpp.o' \
        'vpxor %ymm0, %ymm0, %ymm0'
    expect_classified_archive_rejected lookalike_avx2_t8_k62_prefix \
        'NotLeopard2BackendAVX2T8K62.cpp.o' \
        'vpxor %ymm0, %ymm0, %ymm0'
    expect_classified_archive_rejected lookalike_avx2_t8_k62_suffix \
        'Leopard2BackendAVX2T8K62Extra.cpp.o' \
        'vpxor %ymm0, %ymm0, %ymm0'
    expect_classified_archive_rejected lookalike_avx2_t16_b64_prefix \
        'NotLeopard2BackendAVX2T16B64.cpp.o' \
        'vpxor %ymm0, %ymm0, %ymm0'
    expect_classified_archive_rejected lookalike_avx2_t16_b64_suffix \
        'Leopard2BackendAVX2T16B64Extra.cpp.o' \
        'vpxor %ymm0, %ymm0, %ymm0'
    expect_classified_archive_rejected lookalike_avx2_t16_k66_prefix \
        'NotLeopard2BackendAVX2T16K66.cpp.o' \
        'vpxor %ymm0, %ymm0, %ymm0'
    expect_classified_archive_rejected lookalike_avx2_t16_k66_suffix \
        'Leopard2BackendAVX2T16K66Extra.cpp.o' \
        'vpxor %ymm0, %ymm0, %ymm0'
    expect_classified_archive_rejected lookalike_avx2_t32_b256_prefix \
        'NotLeopard2BackendAVX2T32B256.cpp.o' \
        'vpxor %ymm0, %ymm0, %ymm0'
    expect_classified_archive_rejected lookalike_avx2_t32_b256_suffix \
        'Leopard2BackendAVX2T32B256Extra.cpp.o' \
        'vpxor %ymm0, %ymm0, %ymm0'
    expect_classified_archive_rejected lookalike_gfni_member \
        'NotLeopard2BackendGFNI.cpp.o' \
        'vgf2p8affineqb $0, %ymm0, %ymm0, %ymm0'
    expect_classified_archive_rejected lookalike_probe_member \
        'NotLeopard2CpuFeatures.cpp.o' 'xgetbv'
    expect_duplicate_archive_rejected duplicate_avx2_member \
        'Leopard2BackendAVX2.cpp.o' \
        'vpxor %ymm0, %ymm0, %ymm0' \
        'vfmadd132ps %ymm0, %ymm0, %ymm0'
    expect_missing_expected_member_rejected
    expect_missing_avx2_xor_member_rejected
    expect_missing_avx2_t2_k4_member_rejected
    expect_missing_avx2_t8_k8_b1024_member_rejected
    expect_missing_avx2_t8_k62_member_rejected
    expect_missing_avx2_t16_b64_member_rejected
    expect_missing_avx2_t16_k66_member_rejected
    expect_missing_avx2_t32_b256_member_rejected
    expect_missing_gfni_member_rejected
    expect_missing_avx512_gfni_t128_member_rejected
    expect_listing_failure_rejected

    expect_metadata_rejected flag_make_ssse3 make '-mssse3'
    expect_metadata_rejected flag_compile_sse41 compile_commands '-msse4.1'
    expect_metadata_rejected flag_march make '-march=native'
    expect_metadata_rejected flag_lto make '-flto=auto'

    avx2_xor_source=/src/Leopard2BackendAVX2Xor.cpp
    avx2_xor_output=CMakeFiles/leopard2_backend_avx2.dir/Leopard2BackendAVX2Xor.cpp.o
    expect_fixture_metadata_accepted exact_avx2_xor_source \
        "c++ -mavx2 -mno-avx512f -o $avx2_xor_output -c $avx2_xor_source" \
        "$avx2_xor_source" "$avx2_xor_output"
    expect_fixture_metadata_rejected avx2_xor_missing_flag \
        "c++ -mno-avx512f -o $avx2_xor_output -c $avx2_xor_source" \
        "$avx2_xor_source" "$avx2_xor_output"
    expect_fixture_metadata_rejected avx2_xor_unrelated_isa \
        "c++ -mavx2 -mno-avx512f -mfma -o $avx2_xor_output -c $avx2_xor_source" \
        "$avx2_xor_source" "$avx2_xor_output"
    expect_fixture_metadata_rejected avx2_xor_lookalike_prefix \
        'c++ -mavx2 -mno-avx512f -o CMakeFiles/leopard2_backend_avx2.dir/NotLeopard2BackendAVX2Xor.cpp.o -c /src/NotLeopard2BackendAVX2Xor.cpp' \
        /src/NotLeopard2BackendAVX2Xor.cpp \
        CMakeFiles/leopard2_backend_avx2.dir/NotLeopard2BackendAVX2Xor.cpp.o
    expect_fixture_metadata_rejected avx2_xor_lookalike_suffix \
        'c++ -mavx2 -mno-avx512f -o CMakeFiles/leopard2_backend_avx2.dir/Leopard2BackendAVX2XorXcpp.o -c /src/Leopard2BackendAVX2XorXcpp' \
        /src/Leopard2BackendAVX2XorXcpp \
        CMakeFiles/leopard2_backend_avx2.dir/Leopard2BackendAVX2XorXcpp.o

    t2_k4_source=/src/Leopard2BackendAVX2T2K4.cpp
    t2_k4_output=CMakeFiles/leopard2_backend_avx2_t2_k4.dir/Leopard2BackendAVX2T2K4.cpp.o
    expect_fixture_metadata_accepted exact_avx2_t2_k4_source \
        "c++ -mavx2 -mno-avx512f -o $t2_k4_output -c $t2_k4_source" \
        "$t2_k4_source" "$t2_k4_output"
    expect_fixture_metadata_rejected avx2_t2_k4_missing_avx2 \
        "c++ -mno-avx512f -o $t2_k4_output -c $t2_k4_source" \
        "$t2_k4_source" "$t2_k4_output"
    expect_fixture_metadata_rejected avx2_t2_k4_missing_no_avx512f \
        "c++ -mavx2 -o $t2_k4_output -c $t2_k4_source" \
        "$t2_k4_source" "$t2_k4_output"
    expect_fixture_metadata_rejected avx2_t2_k4_avx2_flag_lookalike \
        "c++ -mavx2extra -mno-avx512f -o $t2_k4_output -c $t2_k4_source" \
        "$t2_k4_source" "$t2_k4_output"
    expect_fixture_metadata_rejected avx2_t2_k4_no_avx512f_flag_lookalike \
        "c++ -mavx2 -mno-avx512fextra -o $t2_k4_output -c $t2_k4_source" \
        "$t2_k4_source" "$t2_k4_output"
    expect_fixture_metadata_rejected avx2_t2_k4_unrelated_fma \
        "c++ -mavx2 -mno-avx512f -mfma -o $t2_k4_output -c $t2_k4_source" \
        "$t2_k4_source" "$t2_k4_output"
    expect_fixture_metadata_rejected avx2_t2_k4_explicit_avx512 \
        "c++ -mavx2 -mno-avx512f -mavx512f -o $t2_k4_output -c $t2_k4_source" \
        "$t2_k4_source" "$t2_k4_output"
    expect_fixture_metadata_rejected avx2_t2_k4_wrong_target \
        "c++ -mavx2 -mno-avx512f -o CMakeFiles/leopard2_backend_avx2.dir/Leopard2BackendAVX2T2K4.cpp.o -c $t2_k4_source" \
        "$t2_k4_source" \
        CMakeFiles/leopard2_backend_avx2.dir/Leopard2BackendAVX2T2K4.cpp.o
    expect_fixture_metadata_rejected avx2_t2_k4_lookalike_prefix \
        'c++ -mavx2 -mno-avx512f -o CMakeFiles/leopard2_backend_avx2_t2_k4.dir/NotLeopard2BackendAVX2T2K4.cpp.o -c /src/NotLeopard2BackendAVX2T2K4.cpp' \
        /src/NotLeopard2BackendAVX2T2K4.cpp \
        CMakeFiles/leopard2_backend_avx2_t2_k4.dir/NotLeopard2BackendAVX2T2K4.cpp.o
    expect_fixture_metadata_rejected avx2_t2_k4_lookalike_suffix \
        'c++ -mavx2 -mno-avx512f -o CMakeFiles/leopard2_backend_avx2_t2_k4.dir/Leopard2BackendAVX2T2K4Extra.cpp.o -c /src/Leopard2BackendAVX2T2K4Extra.cpp' \
        /src/Leopard2BackendAVX2T2K4Extra.cpp \
        CMakeFiles/leopard2_backend_avx2_t2_k4.dir/Leopard2BackendAVX2T2K4Extra.cpp.o

    t8_k8_b1024_source=/src/Leopard2BackendAVX2T8K8B1024.cpp
    t8_k8_b1024_output=CMakeFiles/leopard2_backend_avx2_t8_k8_b1024.dir/Leopard2BackendAVX2T8K8B1024.cpp.o
    expect_fixture_metadata_accepted exact_avx2_t8_k8_b1024_source \
        "c++ -mavx2 -mno-avx512f -o $t8_k8_b1024_output -c $t8_k8_b1024_source" \
        "$t8_k8_b1024_source" "$t8_k8_b1024_output"
    expect_fixture_metadata_rejected avx2_t8_k8_b1024_missing_avx2 \
        "c++ -mno-avx512f -o $t8_k8_b1024_output -c $t8_k8_b1024_source" \
        "$t8_k8_b1024_source" "$t8_k8_b1024_output"
    expect_fixture_metadata_rejected avx2_t8_k8_b1024_missing_no_avx512f \
        "c++ -mavx2 -o $t8_k8_b1024_output -c $t8_k8_b1024_source" \
        "$t8_k8_b1024_source" "$t8_k8_b1024_output"
    expect_fixture_metadata_rejected avx2_t8_k8_b1024_unrelated_fma \
        "c++ -mavx2 -mno-avx512f -mfma -o $t8_k8_b1024_output -c $t8_k8_b1024_source" \
        "$t8_k8_b1024_source" "$t8_k8_b1024_output"
    expect_fixture_metadata_rejected avx2_t8_k8_b1024_wrong_target \
        "c++ -mavx2 -mno-avx512f -o CMakeFiles/leopard2_backend_avx2.dir/Leopard2BackendAVX2T8K8B1024.cpp.o -c $t8_k8_b1024_source" \
        "$t8_k8_b1024_source" \
        CMakeFiles/leopard2_backend_avx2.dir/Leopard2BackendAVX2T8K8B1024.cpp.o
    expect_fixture_metadata_rejected avx2_t8_k8_b1024_lookalike_prefix \
        'c++ -mavx2 -mno-avx512f -o CMakeFiles/leopard2_backend_avx2_t8_k8_b1024.dir/NotLeopard2BackendAVX2T8K8B1024.cpp.o -c /src/NotLeopard2BackendAVX2T8K8B1024.cpp' \
        /src/NotLeopard2BackendAVX2T8K8B1024.cpp \
        CMakeFiles/leopard2_backend_avx2_t8_k8_b1024.dir/NotLeopard2BackendAVX2T8K8B1024.cpp.o
    expect_fixture_metadata_rejected avx2_t8_k8_b1024_lookalike_suffix \
        'c++ -mavx2 -mno-avx512f -o CMakeFiles/leopard2_backend_avx2_t8_k8_b1024.dir/Leopard2BackendAVX2T8K8B1024Extra.cpp.o -c /src/Leopard2BackendAVX2T8K8B1024Extra.cpp' \
        /src/Leopard2BackendAVX2T8K8B1024Extra.cpp \
        CMakeFiles/leopard2_backend_avx2_t8_k8_b1024.dir/Leopard2BackendAVX2T8K8B1024Extra.cpp.o

    t16_source=/src/Leopard2BackendAVX2T16B64.cpp
    t16_output=CMakeFiles/leopard2_backend_avx2_t16_b64.dir/Leopard2BackendAVX2T16B64.cpp.o
    expect_fixture_metadata_accepted exact_avx2_t16_b64_source \
        "c++ -mavx2 -mno-avx512f -o $t16_output -c $t16_source" \
        "$t16_source" "$t16_output"
    expect_fixture_metadata_rejected avx2_t16_b64_missing_avx2 \
        "c++ -mno-avx512f -o $t16_output -c $t16_source" \
        "$t16_source" "$t16_output"
    expect_fixture_metadata_rejected avx2_t16_b64_missing_no_avx512f \
        "c++ -mavx2 -o $t16_output -c $t16_source" \
        "$t16_source" "$t16_output"
    expect_fixture_metadata_rejected avx2_t16_b64_avx2_flag_lookalike \
        "c++ -mavx2extra -mno-avx512f -o $t16_output -c $t16_source" \
        "$t16_source" "$t16_output"
    expect_fixture_metadata_rejected avx2_t16_b64_no_avx512f_flag_lookalike \
        "c++ -mavx2 -mno-avx512fextra -o $t16_output -c $t16_source" \
        "$t16_source" "$t16_output"
    expect_fixture_metadata_rejected avx2_t16_b64_unrelated_fma \
        "c++ -mavx2 -mno-avx512f -mfma -o $t16_output -c $t16_source" \
        "$t16_source" "$t16_output"
    expect_fixture_metadata_rejected avx2_t16_b64_explicit_avx512 \
        "c++ -mavx2 -mno-avx512f -mavx512f -o $t16_output -c $t16_source" \
        "$t16_source" "$t16_output"
    expect_fixture_metadata_rejected avx2_t16_b64_wrong_target \
        "c++ -mavx2 -mno-avx512f -o CMakeFiles/leopard2_backend_avx2.dir/Leopard2BackendAVX2T16B64.cpp.o -c $t16_source" \
        "$t16_source" \
        CMakeFiles/leopard2_backend_avx2.dir/Leopard2BackendAVX2T16B64.cpp.o
    expect_fixture_metadata_rejected avx2_t16_b64_lookalike_prefix \
        'c++ -mavx2 -mno-avx512f -o CMakeFiles/leopard2_backend_avx2_t16_b64.dir/NotLeopard2BackendAVX2T16B64.cpp.o -c /src/NotLeopard2BackendAVX2T16B64.cpp' \
        /src/NotLeopard2BackendAVX2T16B64.cpp \
        CMakeFiles/leopard2_backend_avx2_t16_b64.dir/NotLeopard2BackendAVX2T16B64.cpp.o
    expect_fixture_metadata_rejected avx2_t16_b64_lookalike_suffix \
        'c++ -mavx2 -mno-avx512f -o CMakeFiles/leopard2_backend_avx2_t16_b64.dir/Leopard2BackendAVX2T16B64Extra.cpp.o -c /src/Leopard2BackendAVX2T16B64Extra.cpp' \
        /src/Leopard2BackendAVX2T16B64Extra.cpp \
        CMakeFiles/leopard2_backend_avx2_t16_b64.dir/Leopard2BackendAVX2T16B64Extra.cpp.o

    t16_k66_source=/src/Leopard2BackendAVX2T16K66.cpp
    t16_k66_output=CMakeFiles/leopard2_backend_avx2_t16_k66.dir/Leopard2BackendAVX2T16K66.cpp.o
    expect_fixture_metadata_accepted exact_avx2_t16_k66_source \
        "c++ -mavx2 -mno-avx512f -o $t16_k66_output -c $t16_k66_source" \
        "$t16_k66_source" "$t16_k66_output"
    expect_fixture_metadata_rejected avx2_t16_k66_missing_avx2 \
        "c++ -mno-avx512f -o $t16_k66_output -c $t16_k66_source" \
        "$t16_k66_source" "$t16_k66_output"
    expect_fixture_metadata_rejected avx2_t16_k66_missing_no_avx512f \
        "c++ -mavx2 -o $t16_k66_output -c $t16_k66_source" \
        "$t16_k66_source" "$t16_k66_output"
    expect_fixture_metadata_rejected avx2_t16_k66_avx2_flag_lookalike \
        "c++ -mavx2extra -mno-avx512f -o $t16_k66_output -c $t16_k66_source" \
        "$t16_k66_source" "$t16_k66_output"
    expect_fixture_metadata_rejected avx2_t16_k66_no_avx512f_flag_lookalike \
        "c++ -mavx2 -mno-avx512fextra -o $t16_k66_output -c $t16_k66_source" \
        "$t16_k66_source" "$t16_k66_output"
    expect_fixture_metadata_rejected avx2_t16_k66_unrelated_fma \
        "c++ -mavx2 -mno-avx512f -mfma -o $t16_k66_output -c $t16_k66_source" \
        "$t16_k66_source" "$t16_k66_output"
    expect_fixture_metadata_rejected avx2_t16_k66_wrong_target \
        "c++ -mavx2 -mno-avx512f -o CMakeFiles/leopard2_backend_avx2.dir/Leopard2BackendAVX2T16K66.cpp.o -c $t16_k66_source" \
        "$t16_k66_source" \
        CMakeFiles/leopard2_backend_avx2.dir/Leopard2BackendAVX2T16K66.cpp.o
    expect_fixture_metadata_rejected avx2_t16_k66_lookalike_prefix \
        'c++ -mavx2 -mno-avx512f -o CMakeFiles/leopard2_backend_avx2_t16_k66.dir/NotLeopard2BackendAVX2T16K66.cpp.o -c /src/NotLeopard2BackendAVX2T16K66.cpp' \
        /src/NotLeopard2BackendAVX2T16K66.cpp \
        CMakeFiles/leopard2_backend_avx2_t16_k66.dir/NotLeopard2BackendAVX2T16K66.cpp.o
    expect_fixture_metadata_rejected avx2_t16_k66_lookalike_suffix \
        'c++ -mavx2 -mno-avx512f -o CMakeFiles/leopard2_backend_avx2_t16_k66.dir/Leopard2BackendAVX2T16K66Extra.cpp.o -c /src/Leopard2BackendAVX2T16K66Extra.cpp' \
        /src/Leopard2BackendAVX2T16K66Extra.cpp \
        CMakeFiles/leopard2_backend_avx2_t16_k66.dir/Leopard2BackendAVX2T16K66Extra.cpp.o

    t32_source=/src/Leopard2BackendAVX2T32B256.cpp
    t32_output=CMakeFiles/leopard2_backend_avx2_t32_b256.dir/Leopard2BackendAVX2T32B256.cpp.o
    expect_fixture_metadata_accepted exact_avx2_t32_b256_source \
        "c++ -mavx2 -mno-avx512f -o $t32_output -c $t32_source" \
        "$t32_source" "$t32_output"
    expect_fixture_metadata_rejected avx2_t32_b256_missing_avx2 \
        "c++ -mno-avx512f -o $t32_output -c $t32_source" \
        "$t32_source" "$t32_output"
    expect_fixture_metadata_rejected avx2_t32_b256_missing_no_avx512f \
        "c++ -mavx2 -o $t32_output -c $t32_source" \
        "$t32_source" "$t32_output"
    expect_fixture_metadata_rejected avx2_t32_b256_avx2_flag_lookalike \
        "c++ -mavx2extra -mno-avx512f -o $t32_output -c $t32_source" \
        "$t32_source" "$t32_output"
    expect_fixture_metadata_rejected avx2_t32_b256_no_avx512f_flag_lookalike \
        "c++ -mavx2 -mno-avx512fextra -o $t32_output -c $t32_source" \
        "$t32_source" "$t32_output"
    expect_fixture_metadata_rejected avx2_t32_b256_unrelated_fma \
        "c++ -mavx2 -mno-avx512f -mfma -o $t32_output -c $t32_source" \
        "$t32_source" "$t32_output"
    expect_fixture_metadata_rejected avx2_t32_b256_explicit_avx512 \
        "c++ -mavx2 -mno-avx512f -mavx512f -o $t32_output -c $t32_source" \
        "$t32_source" "$t32_output"
    expect_fixture_metadata_rejected avx2_t32_b256_wrong_target \
        "c++ -mavx2 -mno-avx512f -o CMakeFiles/leopard2_backend_avx2.dir/Leopard2BackendAVX2T32B256.cpp.o -c $t32_source" \
        "$t32_source" \
        CMakeFiles/leopard2_backend_avx2.dir/Leopard2BackendAVX2T32B256.cpp.o
    expect_fixture_metadata_rejected avx2_t32_b256_lookalike_prefix \
        'c++ -mavx2 -mno-avx512f -o CMakeFiles/leopard2_backend_avx2_t32_b256.dir/NotLeopard2BackendAVX2T32B256.cpp.o -c /src/NotLeopard2BackendAVX2T32B256.cpp' \
        /src/NotLeopard2BackendAVX2T32B256.cpp \
        CMakeFiles/leopard2_backend_avx2_t32_b256.dir/NotLeopard2BackendAVX2T32B256.cpp.o
    expect_fixture_metadata_rejected avx2_t32_b256_lookalike_suffix \
        'c++ -mavx2 -mno-avx512f -o CMakeFiles/leopard2_backend_avx2_t32_b256.dir/Leopard2BackendAVX2T32B256Extra.cpp.o -c /src/Leopard2BackendAVX2T32B256Extra.cpp' \
        /src/Leopard2BackendAVX2T32B256Extra.cpp \
        CMakeFiles/leopard2_backend_avx2_t32_b256.dir/Leopard2BackendAVX2T32B256Extra.cpp.o

    gfni_source=/src/Leopard2BackendGFNI.cpp
    gfni_output=CMakeFiles/leopard2_backend_gfni.dir/Leopard2BackendGFNI.cpp.o
    expect_fixture_metadata_accepted exact_gfni_source \
        "c++ -mavx2 -mgfni -mno-avx512f -o $gfni_output -c $gfni_source" \
        "$gfni_source" "$gfni_output"
    expect_fixture_metadata_rejected gfni_missing_gfni_flag \
        "c++ -mavx2 -mno-avx512f -o $gfni_output -c $gfni_source" \
        "$gfni_source" "$gfni_output"
    expect_fixture_metadata_rejected gfni_missing_avx2_flag \
        "c++ -mgfni -mno-avx512f -o $gfni_output -c $gfni_source" \
        "$gfni_source" "$gfni_output"
    expect_fixture_metadata_rejected gfni_missing_no_avx512f_flag \
        "c++ -mavx2 -mgfni -o $gfni_output -c $gfni_source" \
        "$gfni_source" "$gfni_output"
    expect_fixture_metadata_rejected gfni_unrelated_isa \
        "c++ -mavx2 -mgfni -mno-avx512f -mfma -o $gfni_output -c $gfni_source" \
        "$gfni_source" "$gfni_output"
    expect_fixture_metadata_rejected gfni_lookalike_prefix \
        'c++ -mavx2 -mgfni -mno-avx512f -o CMakeFiles/leopard2_backend_gfni.dir/NotLeopard2BackendGFNI.cpp.o -c /src/NotLeopard2BackendGFNI.cpp' \
        /src/NotLeopard2BackendGFNI.cpp \
        CMakeFiles/leopard2_backend_gfni.dir/NotLeopard2BackendGFNI.cpp.o
    expect_fixture_metadata_rejected gfni_lookalike_suffix \
        'c++ -mavx2 -mgfni -mno-avx512f -o CMakeFiles/leopard2_backend_gfni.dir/Leopard2BackendGFNIXcpp.o -c /src/Leopard2BackendGFNIXcpp' \
        /src/Leopard2BackendGFNIXcpp \
        CMakeFiles/leopard2_backend_gfni.dir/Leopard2BackendGFNIXcpp.o

    t128_source=/src/Leopard2BackendAVX512GFNIT128.cpp
    t128_output=CMakeFiles/leopard2_backend_avx512_gfni_t128.dir/Leopard2BackendAVX512GFNIT128.cpp.o
    t128_flags='-mavx2 -mavx512f -mavx512bw -mavx512vl -mgfni'
    expect_fixture_metadata_accepted exact_avx512_gfni_t128_source \
        "c++ $t128_flags -o $t128_output -c $t128_source" \
        "$t128_source" "$t128_output"
    expect_fixture_metadata_rejected t128_missing_gfni \
        "c++ -mavx2 -mavx512f -mavx512bw -mavx512vl -o $t128_output -c $t128_source" \
        "$t128_source" "$t128_output"
    expect_fixture_metadata_rejected t128_missing_avx512bw \
        "c++ -mavx2 -mavx512f -mavx512vl -mgfni -o $t128_output -c $t128_source" \
        "$t128_source" "$t128_output"
    expect_fixture_metadata_rejected t128_avx512f_lookalike_only \
        "c++ -mavx2 -mavx512fp16 -mavx512bw -mavx512vl -mgfni -o $t128_output -c $t128_source" \
        "$t128_source" "$t128_output"
    expect_fixture_metadata_rejected t128_unrelated_avx512fp16 \
        "c++ $t128_flags -mavx512fp16 -o $t128_output -c $t128_source" \
        "$t128_source" "$t128_output"
    expect_fixture_metadata_rejected t128_unrelated_fma \
        "c++ $t128_flags -mfma -o $t128_output -c $t128_source" \
        "$t128_source" "$t128_output"
    expect_fixture_metadata_rejected t128_lookalike_prefix \
        "c++ $t128_flags -o CMakeFiles/leopard2_backend_avx512_gfni_t128.dir/NotLeopard2BackendAVX512GFNIT128.cpp.o -c /src/NotLeopard2BackendAVX512GFNIT128.cpp" \
        /src/NotLeopard2BackendAVX512GFNIT128.cpp \
        CMakeFiles/leopard2_backend_avx512_gfni_t128.dir/NotLeopard2BackendAVX512GFNIT128.cpp.o

    fixture_prefix='c++ -DLEO2_ENABLE_TEST_HOOKS=1 -DLEO2_PORTABLE_ISA_PRIVILEGED_FIXTURE=1 -mssse3 -mno-avx'
    expect_fixture_metadata_accepted exact_legacy_fixture \
        "$fixture_prefix -o CMakeFiles/leopard_legacy_simd_init_test_lib.dir/LeopardFF8.cpp.o -c /src/LeopardFF8.cpp" \
        /src/LeopardFF8.cpp \
        CMakeFiles/leopard_legacy_simd_init_test_lib.dir/LeopardFF8.cpp.o
    expect_fixture_metadata_accepted exact_legacy_fixture_driver \
        "$fixture_prefix -o CMakeFiles/leopard2_legacy_simd_init_failure_test.dir/tests/leopard2/test_legacy_simd_init_failure.cpp.o -c /src/tests/leopard2/test_legacy_simd_init_failure.cpp" \
        /src/tests/leopard2/test_legacy_simd_init_failure.cpp \
        CMakeFiles/leopard2_legacy_simd_init_failure_test.dir/tests/leopard2/test_legacy_simd_init_failure.cpp.o
    expect_fixture_metadata_rejected lookalike_legacy_fixture \
        "$fixture_prefix -o CMakeFiles/not_leopard_legacy_simd_init_test_lib.dir/LeopardFF8.cpp.o -c /src/LeopardFF8.cpp" \
        /src/LeopardFF8.cpp \
        CMakeFiles/not_leopard_legacy_simd_init_test_lib.dir/LeopardFF8.cpp.o
    expect_fixture_metadata_rejected wrong_fixture_source \
        "$fixture_prefix -o CMakeFiles/leopard_legacy_simd_init_test_lib.dir/LeopardFF8.cpp.o -c /src/LeopardFF16.cpp" \
        /src/LeopardFF16.cpp \
        CMakeFiles/leopard_legacy_simd_init_test_lib.dir/LeopardFF8.cpp.o
    expect_fixture_metadata_rejected unreviewed_fixture_source \
        "$fixture_prefix -o CMakeFiles/leopard_legacy_simd_init_test_lib.dir/NewProductionSource.cpp.o -c /src/NewProductionSource.cpp" \
        /src/NewProductionSource.cpp \
        CMakeFiles/leopard_legacy_simd_init_test_lib.dir/NewProductionSource.cpp.o
    expect_fixture_metadata_rejected production_isa_leak \
        "$fixture_prefix -o CMakeFiles/leopard.dir/LeopardFF8.cpp.o -c /src/LeopardFF8.cpp" \
        /src/LeopardFF8.cpp CMakeFiles/leopard.dir/LeopardFF8.cpp.o
    expect_fixture_metadata_rejected fixture_unrelated_isa \
        "$fixture_prefix -msse4.1 -o CMakeFiles/leopard_legacy_simd_init_test_lib.dir/LeopardFF8.cpp.o -c /src/LeopardFF8.cpp" \
        /src/LeopardFF8.cpp \
        CMakeFiles/leopard_legacy_simd_init_test_lib.dir/LeopardFF8.cpp.o
    expect_fixture_metadata_rejected fixture_missing_marker \
        'c++ -DLEO2_ENABLE_TEST_HOOKS=1 -mssse3 -mno-avx -o CMakeFiles/leopard_legacy_simd_init_test_lib.dir/LeopardFF8.cpp.o -c /src/LeopardFF8.cpp' \
        /src/LeopardFF8.cpp \
        CMakeFiles/leopard_legacy_simd_init_test_lib.dir/LeopardFF8.cpp.o
    expect_required_metadata_rejected missing_compile_commands missing
    expect_required_metadata_rejected empty_compile_commands empty
    expect_required_metadata_rejected unrelated_compile_commands unrelated

    echo "portable ISA checker negative controls: PASS"
}

if [ -n "$cc_bin" ] && [ -n "$ar_bin" ]; then
    run_negative_controls
fi

if [ "$self_test_only" -eq 1 ]; then
    echo "portable ISA checker self-test-only: PASS (archive scan intentionally not requested)"
    exit 0
fi

if ! scan_archive "$archive" "${ar_bin:-}" "$expected_classes"; then
    exit 1
fi

if [ -n "$build_dir" ]; then
    if ! scan_build_metadata "$build_dir" "$metadata_mode"; then
        exit 1
    fi
elif [ "$metadata_mode" = required ]; then
    echo "portable ISA check: required compile metadata needs a build directory" >&2
    exit 1
fi

echo "portable ISA check: PASS (SSE2 baseline; named SSSE3/AVX2/GFNI/AVX-512VL/probe members isolated)"

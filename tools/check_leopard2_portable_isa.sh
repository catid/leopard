#!/bin/sh

set -eu

usage()
{
    echo "usage: $0 OBJDUMP STATIC_ARCHIVE [BUILD_DIRECTORY [EXPECTED_CLASSES] | [BUILD_DIRECTORY CC AR [EXPECTED_CLASSES [optional|required]]]]" >&2
    exit 2
}

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
case "$metadata_mode" in
    optional|required) ;;
    *) usage ;;
esac

if [ ! -f "$archive" ]; then
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
allowed_avx2_vex_mnemonics='^(vbroadcastf128|vbroadcasti128|vinserti128|vmovd|vmovdqa|vmovdqu|vmovq|vmovups|vpand|vpbroadcastb|vpbroadcastq|vpinsrb|vpshufb|vpsrlq|vpunpckldq|vpunpcklqdq|vpunpcklwd|vpxor|vxorps|vzeroupper)$'

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
            Leopard2BackendAVX2.cpp.o|Leopard2BackendAVX2.cpp.obj)
                object_class=avx2 ;;
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
        violating_lines="$scratch_root/metadata-violations"
        LC_ALL=C grep -Ein -- "$forbidden_metadata" "$compile_commands" |
            grep -Ev 'Leopard2Backend(SSSE3|AVX2)[.]cpp' > "$violating_lines" || true
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
        if grep -q 'Leopard2BackendAVX2[.]cpp' "$compile_commands" &&
           ! grep 'Leopard2BackendAVX2[.]cpp' "$compile_commands" |
               grep -q -- '-mavx2'
        then
            echo "portable ISA check: AVX2 object lacks its ISA flag" >&2
            return 1
        fi
        avx2_command="$scratch_root/avx2-command"
        grep 'Leopard2BackendAVX2[.]cpp' "$compile_commands" |
            sed 's/-mavx2//g' > "$avx2_command" || true
        if [ -s "$avx2_command" ] &&
           grep -Ein -- "$forbidden_metadata" "$avx2_command"
        then
            echo "portable ISA check: AVX2 object has an unrelated ISA/LTO flag" >&2
            return 1
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

    expect_classified_archive_accepted good_probe \
        'Leopard2CpuFeatures.cpp.o' 'xgetbv'
    expect_classified_archive_accepted good_ssse3 \
        'Leopard2BackendSSSE3.cpp.o' 'pshufb %xmm0, %xmm0'
    expect_classified_archive_accepted good_avx2 \
        'Leopard2BackendAVX2.cpp.o' 'vpxor %ymm0, %ymm0, %ymm0'
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
    expect_classified_archive_rejected avx2_leaks_fma \
        'Leopard2BackendAVX2.cpp.o' \
        'vfmadd132ps %ymm0, %ymm0, %ymm0'
    expect_classified_archive_rejected avx2_leaks_f16c \
        'Leopard2BackendAVX2.cpp.o' 'vcvtph2ps %xmm0, %ymm0'
    expect_classified_archive_rejected avx2_leaks_evex \
        'Leopard2BackendAVX2.cpp.o' \
        'vpternlogd $0, %ymm0, %ymm0, %ymm0'
    expect_classified_archive_rejected probe_leaks_avx2 \
        'Leopard2CpuFeatures.cpp.o' 'vpaddd %ymm0, %ymm0, %ymm0'
    expect_classified_archive_rejected ssse3_leaks_avx2 \
        'Leopard2BackendSSSE3.cpp.o' 'vpaddd %ymm0, %ymm0, %ymm0'
    expect_classified_archive_rejected avx2_leaks_probe \
        'Leopard2BackendAVX2.cpp.o' 'xgetbv'
    expect_classified_archive_rejected avx2_leaks_avx512 \
        'Leopard2BackendAVX2.cpp.o' 'vpxord %zmm0, %zmm0, %zmm0'
    expect_classified_archive_rejected lookalike_ssse3_member \
        'NotLeopard2BackendSSSE3.cpp.o' 'pshufb %xmm0, %xmm0'
    expect_classified_archive_rejected lookalike_avx2_member \
        'NotLeopard2BackendAVX2.cpp.o' 'vpaddd %ymm0, %ymm0, %ymm0'
    expect_classified_archive_rejected lookalike_probe_member \
        'NotLeopard2CpuFeatures.cpp.o' 'xgetbv'
    expect_duplicate_archive_rejected duplicate_avx2_member \
        'Leopard2BackendAVX2.cpp.o' \
        'vpxor %ymm0, %ymm0, %ymm0' \
        'vfmadd132ps %ymm0, %ymm0, %ymm0'
    expect_missing_expected_member_rejected
    expect_listing_failure_rejected

    expect_metadata_rejected flag_make_ssse3 make '-mssse3'
    expect_metadata_rejected flag_compile_sse41 compile_commands '-msse4.1'
    expect_metadata_rejected flag_march make '-march=native'
    expect_metadata_rejected flag_lto make '-flto=auto'
    expect_required_metadata_rejected missing_compile_commands missing
    expect_required_metadata_rejected empty_compile_commands empty
    expect_required_metadata_rejected unrelated_compile_commands unrelated

    echo "portable ISA checker negative controls: PASS"
}

if [ -n "$cc_bin" ] && [ -n "$ar_bin" ]; then
    run_negative_controls
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

echo "portable ISA check: PASS (SSE2 baseline; named SSSE3/AVX2/probe members isolated)"

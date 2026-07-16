#!/bin/sh

set -eu

usage()
{
    echo "usage: $0 OBJDUMP STATIC_ARCHIVE [BUILD_DIRECTORY [CC AR]]" >&2
    exit 2
}

case "$#" in
    2|3|5) ;;
    *) usage ;;
esac

objdump_bin=$1
archive=$2
build_dir=${3:-}
cc_bin=${4:-}
ar_bin=${5:-}

if [ ! -f "$archive" ]; then
    echo "portable ISA check: archive not found: $archive" >&2
    exit 2
fi

# This is intentionally an x86-64 post-baseline denylist boundary: SSE2
# instructions are allowed and known compiler-generated extension families are
# rejected. Matching is performed only on objdump's mnemonic field, so symbol
# names and operands cannot trigger a false positive.
forbidden_mnemonics='^(addsubp[ds]|haddp[ds]|hsubp[ds]|lddqu|movddup|movshdup|movsldup|fisttp[[:alnum:]]*|monitor|mwait|monitorx|mwaitx|prefetchw|prefetchwt1|lahf|sahf|pabs[bdw]|palignr|phadd(d|sw|w)|phsub(d|sw|w)|pmaddubsw|pmulhrsw|pshufb|psign[bdw]|blendp[ds]|blendvp[ds]|dpp[ds]|extractps|insertps|movntdqa|mpsadbw|packusdw|pblendvb|pblendw|pcmpeqq|pextr[bdq]|phminposuw|pinsr[bdq]|pmax(sb|sd|uw|ud)|pmin(sb|sd|uw|ud)|pmov(sx|zx)(bd|bq|bw|dq|dw|wd|wq)|pmuldq|pmulld|ptest|round(p[ds]|s[ds])|crc32[[:alnum:]]*|pcmp(e|i)str[im]|pcmpgtq|popcnt|extrq|insertq|movntsd|movntss|v[[:alnum:]_.]*|cmpxchg16b|lzcnt|tzcnt|andn|bextr|blcfill|blci|blcic|blcmsk|blcs|blsfill|blsic|t1mskc|tzmsk|blsi|blsmsk|blsr|bzhi|pdep|pext|mulx|rorx|sarx|shlx|shrx|adcx|adox|movbe|rdrand|rdseed|rdpid|rdtscp|clflushopt|clwb|clzero|cldemote|wbnoinvd|serialize|umonitor|umwait|tpause|movdiri|movdir64b|enqcmd|enqcmds|xgetbv|xsetbv|xsave[[:alnum:]]*|xrstor[[:alnum:]]*|aes(enc|enclast|dec|declast|imc|keygenassist)|pclmul[[:alnum:]]*|gf2p8[[:alnum:]]*|sha1(msg1|msg2|nexte|rnds4)|sha256(msg1|msg2|rnds2)|k(add|and|andn|mov|not|or|ortest|shiftl|shiftr|test|unpck|xnor|xor)[bdqw]*)$'

# Reject target-raising options in Make, Ninja, or compilation-database
# metadata.  -mno-* and the x86-64 SSE2 baseline remain allowed.  All -march
# and -mcpu values are conservatively rejected because proving that an
# arbitrary spelling retains the baseline is toolchain-specific.
forbidden_flags='(^|[[:space:]"=,:])(-march=[^[:space:]",]+|-mcpu=[^[:space:]",]+|-m(sse3|ssse3|sse4|sse4a|sse4[.][12]|avx[0-9.]*|avx512[^[:space:]",]*|fma|fma4|f16c|bmi|bmi2|lzcnt|popcnt|cx16|lahf-sahf|prfchw|mwaitx|clzero|wbnoinvd|aes|pclmul|vpclmulqdq|gfni|vaes|sha|adx|movbe|xsave[^[:space:]",]*|rdrnd|rdseed|xop|tbm|3dnow)([[:space:]",]|$)|/arch:(avx|avx2|avx512[^[:space:]",]*|sse3|ssse3|sse4([.][12])?)([[:space:]",]|$)|[+](sse3|ssse3|sse4a?|sse4[.]?[12]|avx[0-9.]*|avx512[^[:space:]",]*|fma|f16c|bmi2?|lzcnt|popcnt|cx16|lahf-sahf|prfchw|mwaitx|clzero|wbnoinvd|aes|pclmul|gfni|vaes|sha|adx|movbe)([[:space:]",]|$))'

scratch_root=$(mktemp -d "${TMPDIR:-/tmp}/leopard2-portable-isa.XXXXXX")
trap 'rm -rf "$scratch_root"' EXIT HUP INT TERM
scan_index=0

scan_archive()
{
    scan_index=$((scan_index + 1))
    disassembly_file="$scratch_root/disassembly.$scan_index"
    mnemonics_file="$scratch_root/mnemonics.$scan_index"

    LC_ALL=C "$objdump_bin" -d --no-show-raw-insn "$1" > "$disassembly_file"
    sed -nE 's/^[[:space:]]*[0-9a-f]+:[[:space:]]+([[:alnum:]_.]+).*/\1/p' \
        "$disassembly_file" | sort -u > "$mnemonics_file"

    if [ ! -s "$mnemonics_file" ]; then
        echo "portable ISA check: objdump produced no instruction mnemonics" >&2
        return 1
    fi
    if grep -E "$forbidden_mnemonics" "$mnemonics_file"; then
        echo "portable ISA check: archive contains a post-SSE2 instruction" >&2
        return 1
    fi
    return 0
}

scan_build_metadata()
{
    metadata_root=$1
    for metadata_file in \
        "$metadata_root/CMakeFiles/libleopard.dir/flags.make" \
        "$metadata_root/compile_commands.json" \
        "$metadata_root/build.ninja"
    do
        if [ -f "$metadata_file" ] &&
           LC_ALL=C grep -Ein -- "$forbidden_flags" "$metadata_file"
        then
            echo "portable ISA check: compiler metadata raises the ISA floor: $metadata_file" >&2
            return 1
        fi
    done
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

expect_metadata_rejected()
{
    fixture_name=$1
    metadata_kind=$2
    fixture_flag=$3
    fixture_dir="$scratch_root/metadata-$fixture_name"
    mkdir -p "$fixture_dir/CMakeFiles/libleopard.dir"

    case "$metadata_kind" in
        make)
            printf 'CXX_FLAGS = -O2 %s\n' "$fixture_flag" > \
                "$fixture_dir/CMakeFiles/libleopard.dir/flags.make"
            ;;
        compile_commands)
            printf '[{"command":"c++ -O2 %s -c x.cpp","file":"x.cpp"}]\n' \
                "$fixture_flag" > "$fixture_dir/compile_commands.json"
            ;;
        ninja)
            printf 'FLAGS = -O2 %s\n' "$fixture_flag" > \
                "$fixture_dir/build.ninja"
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

    expect_metadata_rejected flag_make_ssse3 make '-mssse3'
    expect_metadata_rejected flag_compile_sse41 compile_commands '-msse4.1'
    expect_metadata_rejected flag_ninja_avx2 ninja '-mavx2'
    expect_metadata_rejected flag_march make '-march=native'

    echo "portable ISA checker negative controls: PASS"
}

if [ -n "$cc_bin" ] && [ -n "$ar_bin" ]; then
    run_negative_controls
fi

if ! scan_archive "$archive"; then
    exit 1
fi

if [ -n "$build_dir" ] && ! scan_build_metadata "$build_dir"; then
    exit 1
fi

echo "portable ISA check: PASS (x86-64 SSE2 baseline)"

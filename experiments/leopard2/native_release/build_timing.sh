#!/usr/bin/env bash
# Run under the canonical lock and a 512-MiB/no-swap, single-job build scope.
set -euxo pipefail
[[ $# == 2 ]] || { echo 'usage: build_timing.sh qualified_lane fresh_output' >&2; exit 1; }
qualified=$(realpath -e -- "$1")
lane=$(realpath -e -- "$2")
source_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
[[ $(git -C "$qualified/candidate-source" rev-parse HEAD) == e35b1f0c3a5ef9f241ad3b174b0fbfdb61bbccab ]]
[[ $(git -C "$qualified/baseline-source" rev-parse HEAD) == 6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198 ]]
[[ -z $(git -C "$qualified/candidate-source" status --porcelain) ]]
[[ -z $(git -C "$qualified/baseline-source" status --porcelain) ]]
printf '%s  %s\n' d78cf8ef942b3ee93063b69c186e1512d8fd224e31aa2ae60ab72da8b6374847 "$qualified/frozen/current.a" b183093f2d6dd935cc1119a33fac5e08dec8c654836402b5152a1dc85ae2b285 "$qualified/frozen/main.a" | sha256sum --check -
mkdir "$lane/source" "$lane/build" "$lane/frozen"
cp "$source_dir/encode_probe.cpp" "$source_dir/encode_timing.cpp" "$source_dir/EncodeTiming.h" "$source_dir/timing_clock.cpp" "$source_dir/timing_witness.cpp" "$source_dir/test_encode_timing.cpp" "$lane/source/"
cp "$qualified/frozen/main.a" "$qualified/frozen/current.a" "$lane/frozen/"
common=(-std=c++11 -Wall -Wextra -Werror -O3 -fopenmp)
/usr/bin/c++ "${common[@]}" -DNDEBUG -I "$qualified/candidate-source" '-DLEO_NATIVE_RELEASE_CODEC_COMMIT="e35b1f0c3a5ef9f241ad3b174b0fbfdb61bbccab"' -c "$lane/source/encode_timing.cpp" -o "$lane/build/current.o"
/usr/bin/c++ "${common[@]}" -march=native -g -O0 -O3 -I "$qualified/baseline-source" '-DLEO_NATIVE_RELEASE_CODEC_COMMIT="6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198"' -DLEO_NATIVE_RELEASE_BASELINE=1 -c "$lane/source/encode_timing.cpp" -o "$lane/build/main.o"
for kind in steady synthetic abort; do
    clock_flags=()
    if [[ $kind == synthetic ]]; then clock_flags=(-DLEO_NATIVE_CLOCK_SYNTHETIC=1); fi
    if [[ $kind == abort ]]; then clock_flags=(-DLEO_NATIVE_CLOCK_ABORT=1); fi
    /usr/bin/c++ "${common[@]}" "${clock_flags[@]}" -c "$lane/source/timing_clock.cpp" -o "$lane/build/clock-$kind.o"
done
/usr/bin/c++ "${common[@]}" -I "$qualified/candidate-source" -c "$lane/source/timing_witness.cpp" -o "$lane/build/witness-current.o"
/usr/bin/c++ "${common[@]}" -I "$qualified/baseline-source" -DLEO_NATIVE_RELEASE_BASELINE=1 -c "$lane/source/timing_witness.cpp" -o "$lane/build/witness-main.o"
for implementation in main current; do
    for kind in steady synthetic abort; do
        extra=()
        if [[ $kind == synthetic ]]; then
            symbol=leo2_encode
            if [[ $implementation == main ]]; then symbol=leo_encode; fi
            extra=("$lane/build/witness-$implementation.o" "-Wl,--wrap=$symbol")
        fi
        /usr/bin/c++ "$lane/build/$implementation.o" "${extra[@]}" "$lane/frozen/$implementation.a" "$lane/build/clock-$kind.o" -fopenmp -o "$lane/frozen/$implementation-$kind"
    done
done
/usr/bin/c++ "${common[@]}" "$lane/source/test_encode_timing.cpp" -o "$lane/frozen/timing-unit"
/usr/bin/c++ "${common[@]}" -O1 -g1 -fsanitize=address,undefined -fno-omit-frame-pointer -fno-pie -no-pie "$lane/source/test_encode_timing.cpp" -o "$lane/frozen/timing-unit-sanitized"
cd "$lane/frozen"
chmod 555 main-steady main-synthetic main-abort current-steady current-synthetic current-abort timing-unit timing-unit-sanitized
chmod 444 main.a current.a
sha256sum main-steady main-synthetic main-abort current-steady current-synthetic current-abort timing-unit timing-unit-sanitized main.a current.a > SHA256SUMS
chmod 444 SHA256SUMS
chmod 555 "$lane/frozen"

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

#include "direct_oracle.h"

#include <algorithm>
#include <stdint.h>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace {

using leopard2_test::BinaryField;
using leopard2_test::Element;
using leopard2_test::Matrix;
using leopard2_test::ProfileLayout;

typedef std::vector<uint8_t> Bytes;
typedef std::vector<Bytes> Shards;

void require(bool condition, const std::string& message)
{
    if (!condition)
        throw std::runtime_error(message);
}

size_t rounded_bytes(size_t bytes)
{
    return (bytes + 63u) & ~static_cast<size_t>(63u);
}

Bytes compact_pack(const Bytes& input)
{
    require(!input.empty() && (input.size() & 1u) == 0,
        "compact GF16 packing requires a positive even byte count");
    Bytes output(rounded_bytes(input.size()), 0);
    const size_t complete = input.size() & ~static_cast<size_t>(63u);
    std::copy(input.begin(), input.begin() + complete, output.begin());
    const size_t residual = input.size() - complete;
    const size_t symbols = residual / 2;
    for (size_t i = 0; i < symbols; ++i)
    {
        output[complete + i] = input[complete + i];
        output[complete + 32 + i] = input[complete + symbols + i];
    }
    return output;
}

Bytes compact_gather(const Bytes& input, size_t bytes)
{
    require(bytes != 0 && (bytes & 1u) == 0 && input.size() == rounded_bytes(bytes),
        "invalid compact GF16 gather");
    Bytes output(bytes, 0);
    const size_t complete = bytes & ~static_cast<size_t>(63u);
    std::copy(input.begin(), input.begin() + complete, output.begin());
    const size_t residual = bytes - complete;
    const size_t symbols = residual / 2;
    for (size_t i = 0; i < symbols; ++i)
    {
        output[complete + i] = input[complete + i];
        output[complete + symbols + i] = input[complete + 32 + i];
    }
    return output;
}

Element load_altmap(const Bytes& data, size_t tile, size_t lane)
{
    return static_cast<Element>(data[tile + lane] |
        (static_cast<unsigned>(data[tile + 32 + lane]) << 8));
}

void store_altmap(Bytes* data, size_t tile, size_t lane, Element value)
{
    (*data)[tile + lane] = static_cast<uint8_t>(value);
    (*data)[tile + 32 + lane] = static_cast<uint8_t>(value >> 8);
}

unsigned binary_rank(std::vector<std::vector<uint8_t> > matrix)
{
    if (matrix.empty())
        return 0;
    const unsigned rows = static_cast<unsigned>(matrix.size());
    const unsigned columns = static_cast<unsigned>(matrix[0].size());
    unsigned rank = 0;
    for (unsigned column = 0; column < columns && rank < rows; ++column)
    {
        unsigned pivot = rank;
        while (pivot < rows && matrix[pivot][column] == 0)
            ++pivot;
        if (pivot == rows)
            continue;
        std::swap(matrix[pivot], matrix[rank]);
        for (unsigned row = 0; row < rows; ++row)
            if (row != rank && matrix[row][column])
                for (unsigned j = column; j < columns; ++j)
                    matrix[row][j] ^= matrix[rank][j];
        ++rank;
    }
    return rank;
}

std::vector<std::vector<uint8_t> > projected_binary_matrix(
    const BinaryField& field,
    const Matrix& generator,
    const std::vector<unsigned>& selected_rows)
{
    const unsigned k = static_cast<unsigned>(generator[0].size());
    std::vector<std::vector<uint8_t> > result(
        selected_rows.size() * 8, std::vector<uint8_t>(k * 8, 0));
    for (unsigned selected = 0; selected < selected_rows.size(); ++selected)
    {
        const std::vector<Element>& row = generator[selected_rows[selected]];
        for (unsigned original = 0; original < k; ++original)
            for (unsigned input_bit = 0; input_bit < 8; ++input_bit)
            {
                const Element product = field.multiply(
                    row[original], static_cast<Element>(1u << input_bit));
                for (unsigned output_bit = 0; output_bit < 8; ++output_bit)
                    result[selected * 8 + output_bit][original * 8 + input_bit] =
                        static_cast<uint8_t>((product >> output_bit) & 1u);
            }
    }
    return result;
}

std::vector<uint8_t> nonzero_null_vector(
    std::vector<std::vector<uint8_t> > matrix)
{
    const unsigned rows = static_cast<unsigned>(matrix.size());
    const unsigned columns = static_cast<unsigned>(matrix[0].size());
    std::vector<int> pivot_column(rows, -1);
    std::vector<uint8_t> is_pivot(columns, 0);
    unsigned rank = 0;
    for (unsigned column = 0; column < columns && rank < rows; ++column)
    {
        unsigned pivot = rank;
        while (pivot < rows && matrix[pivot][column] == 0)
            ++pivot;
        if (pivot == rows)
            continue;
        std::swap(matrix[pivot], matrix[rank]);
        for (unsigned row = 0; row < rows; ++row)
            if (row != rank && matrix[row][column])
                for (unsigned j = column; j < columns; ++j)
                    matrix[row][j] ^= matrix[rank][j];
        pivot_column[rank] = static_cast<int>(column);
        is_pivot[column] = 1;
        ++rank;
    }
    require(rank < columns, "null-vector request was made for a full-rank matrix");
    unsigned free_column = 0;
    while (is_pivot[free_column])
        ++free_column;
    std::vector<uint8_t> vector(columns, 0);
    vector[free_column] = 1;
    for (unsigned row = rank; row-- > 0;)
    {
        const unsigned pivot = static_cast<unsigned>(pivot_column[row]);
        uint8_t sum = 0;
        for (unsigned column = pivot + 1; column < columns; ++column)
            sum ^= static_cast<uint8_t>(matrix[row][column] & vector[column]);
        vector[pivot] = sum;
    }
    bool nonzero = false;
    for (unsigned column = 0; column < columns; ++column)
        nonzero = nonzero || vector[column] != 0;
    require(nonzero, "constructed binary null vector is zero");
    for (unsigned row = 0; row < rows; ++row)
    {
        uint8_t sum = 0;
        for (unsigned column = 0; column < columns; ++column)
            sum ^= static_cast<uint8_t>(matrix[row][column] & vector[column]);
        require(sum == 0, "constructed vector is not in the binary nullspace");
    }
    return vector;
}

std::vector<std::vector<unsigned> > coordinate_subsets(unsigned n, unsigned k)
{
    std::vector<std::vector<unsigned> > result;
    std::vector<unsigned> subset(k);
    for (unsigned i = 0; i < k; ++i)
        subset[i] = i;
    for (;;)
    {
        result.push_back(subset);
        int i = static_cast<int>(k) - 1;
        while (i >= 0 && subset[static_cast<unsigned>(i)] == n - k + static_cast<unsigned>(i))
            --i;
        if (i < 0)
            break;
        ++subset[static_cast<unsigned>(i)];
        for (unsigned j = static_cast<unsigned>(i) + 1; j < k; ++j)
            subset[j] = subset[j - 1] + 1;
    }
    return result;
}

Shards deterministic_shards(unsigned count, size_t bytes)
{
    Shards shards(count, Bytes(bytes, 0));
    uint64_t state = UINT64_C(0x9e3779b97f4a7c15) ^ bytes;
    for (unsigned shard = 0; shard < count; ++shard)
        for (size_t i = 0; i < bytes; ++i)
        {
            state ^= state >> 12;
            state ^= state << 25;
            state ^= state >> 27;
            shards[shard][i] = static_cast<uint8_t>(
                (state * UINT64_C(2685821657736338717)) >> 56);
        }
    return shards;
}

uint64_t test_compact_even_tails(const BinaryField& field)
{
    const unsigned k = 3;
    const unsigned r = 5;
    const ProfileLayout layout = leopard2_test::make_profile_layout(
        leopard2_test::kLow, k, r);
    const Matrix generator = leopard2_test::direct_systematic_generator(field, layout);
    const std::vector<std::vector<unsigned> > subsets = coordinate_subsets(k + r, k);
    std::vector<Matrix> inverses;
    for (unsigned subset_i = 0; subset_i < subsets.size(); ++subset_i)
    {
        Matrix selected(k, std::vector<Element>(k, 0));
        for (unsigned row = 0; row < k; ++row)
            selected[row] = generator[subsets[subset_i][row]];
        Matrix inverse;
        require(leopard2_test::invert_matrix(field, selected, &inverse),
            "direct GF16 generator unexpectedly has a singular K-subset");
        inverses.push_back(inverse);
    }

    std::vector<size_t> byte_counts;
    for (size_t bytes = 2; bytes <= 64; bytes += 2)
        byte_counts.push_back(bytes);
    byte_counts.push_back(66);
    byte_counts.push_back(1024);
    byte_counts.push_back(1026);

    uint64_t recovered_symbols = 0;
    for (unsigned count_i = 0; count_i < byte_counts.size(); ++count_i)
    {
        const size_t bytes = byte_counts[count_i];
        const Shards originals = deterministic_shards(k, bytes);
        Shards packed(k);
        for (unsigned i = 0; i < k; ++i)
        {
            packed[i] = compact_pack(originals[i]);
            require(compact_gather(packed[i], bytes) == originals[i],
                "compact GF16 pack/gather did not round trip");
        }
        Shards codeword(k + r, Bytes(rounded_bytes(bytes), 0));
        for (unsigned i = 0; i < k; ++i)
            codeword[i] = packed[i];

        for (size_t tile = 0; tile < rounded_bytes(bytes); tile += 64)
            for (size_t lane = 0; lane < 32; ++lane)
            {
                std::vector<Element> message(k, 0);
                for (unsigned i = 0; i < k; ++i)
                    message[i] = load_altmap(packed[i], tile, lane);
                const std::vector<Element> encoded =
                    leopard2_test::matrix_vector_multiply(field, generator, message);
                for (unsigned row = k; row < k + r; ++row)
                    store_altmap(&codeword[row], tile, lane, encoded[row]);

                for (unsigned subset_i = 0; subset_i < subsets.size(); ++subset_i)
                {
                    std::vector<Element> received(k, 0);
                    for (unsigned row = 0; row < k; ++row)
                        received[row] = encoded[subsets[subset_i][row]];
                    const std::vector<Element> recovered =
                        leopard2_test::matrix_vector_multiply(
                            field, inverses[subset_i], received);
                    require(recovered == message,
                        "compact GF16 tail failed direct K-subset recovery");
                    recovered_symbols += k;
                }
            }

        for (unsigned row = 0; row < k + r; ++row)
        {
            const Bytes compact = compact_gather(codeword[row], bytes);
            require(compact_pack(compact) == codeword[row],
                "compact parity did not preserve every active GF16 symbol");
        }
        if ((bytes & 63u) == 0)
            for (unsigned i = 0; i < k; ++i)
                require(packed[i] == originals[i],
                    "complete ALTMAP tile changed under compact packing");
    }
    return recovered_symbols;
}

std::string test_one_byte_projection_is_not_mds(const BinaryField& field)
{
    const ProfileLayout layout = leopard2_test::make_profile_layout(
        leopard2_test::kLow, 2, 256);
    require(layout.parent_size == 512,
        "one-byte impossibility witness did not select the intended GF16 parent");
    const Matrix generator = leopard2_test::direct_systematic_generator(field, layout);
    for (unsigned first = 0; first < generator.size(); ++first)
        for (unsigned second = first + 1; second < generator.size(); ++second)
        {
            const std::vector<unsigned> selected = { first, second };
            const std::vector<std::vector<uint8_t> > projected =
                projected_binary_matrix(field, generator, selected);
            const unsigned rank = binary_rank(projected);
            if (rank == 16)
                continue;
            const std::vector<uint8_t> collision = nonzero_null_vector(projected);
            unsigned message0 = 0;
            unsigned message1 = 0;
            for (unsigned bit = 0; bit < 8; ++bit)
            {
                message0 |= static_cast<unsigned>(collision[bit]) << bit;
                message1 |= static_cast<unsigned>(collision[8 + bit]) << bit;
            }
            require(message0 != 0 || message1 != 0,
                "projected collision did not change the source message");
            Element full_output0 = 0;
            Element full_output1 = 0;
            const Element messages[2] = {
                static_cast<Element>(message0), static_cast<Element>(message1)
            };
            for (unsigned original = 0; original < 2; ++original)
            {
                full_output0 ^= field.multiply(generator[first][original], messages[original]);
                full_output1 ^= field.multiply(generator[second][original], messages[original]);
            }
            require((full_output0 & 0xffu) == 0 && (full_output1 & 0xffu) == 0,
                "binary null vector did not collide after low-byte projection");
            require(full_output0 != 0 || full_output1 != 0,
                "projection witness did not discard any GF16 information");
            std::ostringstream stream;
            stream << "rows=" << first << ',' << second
                   << " rank=" << rank
                   << " colliding_message=0x" << std::hex << message0
                   << ",0x" << message1
                   << " full_outputs=0x" << full_output0
                   << ",0x" << full_output1;
            return stream.str();
        }
    throw std::runtime_error(
        "projected low-byte [258,2] code unexpectedly passed every pair");
}

} // namespace

int main()
{
    try
    {
        const BinaryField field = leopard2_test::make_legacy_gf16();
        const uint64_t recovered = test_compact_even_tails(field);
        const std::string witness = test_one_byte_projection_is_not_mds(field);
        std::cout << "GF16 tail oracle passed: recovered_symbols=" << recovered
                  << " unsafe_projection_witness=" << witness << std::endl;
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "GF16 tail oracle failed: " << error.what() << std::endl;
        return 1;
    }
}

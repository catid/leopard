/*
    Copyright (c) 2017 Christopher A. Taylor.  All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the conditions in the Leopard-RS
    BSD license are met.  See the repository LICENSE file.
*/

// C2 research checkpoint: A compact parent-preserving truncated LCH schedule.
//
// This standalone program is deliberately not part of the default build.  It
// calls the production full-subtree kernels through their test-only hooks, but
// executes ragged boundary butterflies with an independent scalar field.  It
// compares every requested byte with the complete padded parent transform.

#include "LeopardFF8.h"
#include "LeopardFF16.h"
#include "leopard2.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#if defined(__linux__)
#include <sched.h>
#include <sys/utsname.h>
#include <unistd.h>
#endif

#if defined(_OPENMP)
#include <omp.h>
#endif

#if !defined(LEO2_C2_SOURCE_SHA256)
#define LEO2_C2_SOURCE_SHA256 "unbound"
#endif

#if !defined(LEO2_C2_CORE_GIT_SHA)
#define LEO2_C2_CORE_GIT_SHA "unbound"
#endif

#if !defined(LEO2_C2_LIBRARY_SHA256)
#define LEO2_C2_LIBRARY_SHA256 "unbound"
#endif

#if !defined(LEO2_C2_CORE_MATRIX_SHA256)
#define LEO2_C2_CORE_MATRIX_SHA256 "unbound"
#endif

#if !defined(LEO2_C2_SANITIZER_MODE)
#define LEO2_C2_SANITIZER_MODE "unbound"
#endif

namespace {

typedef uint16_t Element;
typedef std::array<uint8_t, 64> Shard;
typedef std::vector<uint8_t> Mask;

static const int kMissing = -1;

void require(bool condition, const std::string& message)
{
    if (!condition)
        throw std::runtime_error(message);
}

class BinaryField
{
public:
    BinaryField(
        unsigned bits,
        uint32_t polynomial,
        const Element* basis)
        : bits_(bits)
        , order_(1u << bits)
        , polynomial_(polynomial)
        , coordinate_to_polynomial_(order_, 0)
        , polynomial_to_coordinate_(order_, UINT32_MAX)
    {
        require(bits == 8 || bits == 16, "C2 supports only legacy GF8/GF16");
        require((polynomial & (1u << bits)) != 0,
            "field polynomial has the wrong degree");
        for (unsigned coordinate = 0; coordinate < order_; ++coordinate)
        {
            Element polynomial_value = 0;
            for (unsigned bit = 0; bit < bits_; ++bit)
                if ((coordinate & (1u << bit)) != 0)
                    polynomial_value ^= basis[bit];
            require(polynomial_value < order_, "basis value outside field");
            require(polynomial_to_coordinate_[polynomial_value] == UINT32_MAX,
                "coordinate basis is not independent");
            coordinate_to_polynomial_[coordinate] = polynomial_value;
            polynomial_to_coordinate_[polynomial_value] = coordinate;
        }
    }

    unsigned bits() const { return bits_; }
    unsigned order() const { return order_; }
    unsigned lanes() const { return bits_ == 8 ? 64u : 32u; }

    Element multiply(Element left, Element right) const
    {
        uint32_t a = coordinate_to_polynomial_[left];
        uint32_t b = coordinate_to_polynomial_[right];
        uint32_t product = 0;
        while (b != 0)
        {
            if ((b & 1u) != 0)
                product ^= a;
            b >>= 1;
            a <<= 1;
        }
        for (int bit = static_cast<int>(2 * bits_ - 2);
             bit >= static_cast<int>(bits_); --bit)
        {
            if ((product & (1u << bit)) != 0)
                product ^= polynomial_ << (bit - bits_);
        }
        if (product >= order_ ||
            polynomial_to_coordinate_[product] == UINT32_MAX)
            throw std::logic_error("field reduction escaped the coordinate map");
        return static_cast<Element>(polynomial_to_coordinate_[product]);
    }

    Element power(Element value, uint32_t exponent) const
    {
        Element result = 1;
        while (exponent != 0)
        {
            if ((exponent & 1u) != 0)
                result = multiply(result, value);
            exponent >>= 1;
            if (exponent != 0)
                value = multiply(value, value);
        }
        return result;
    }

    Element inverse(Element value) const
    {
        require(value != 0, "zero has no inverse");
        return power(value, order_ - 2);
    }

    Element get(const Shard& shard, unsigned lane) const
    {
        if (lane >= lanes())
            throw std::out_of_range("lane outside shard");
        if (bits_ == 8)
            return shard[lane];
        return static_cast<Element>(shard[lane] |
            (static_cast<unsigned>(shard[32 + lane]) << 8));
    }

    void set(Shard* shard, unsigned lane, Element value) const
    {
        if (lane >= lanes())
            throw std::out_of_range("lane outside shard");
        if (bits_ == 8)
        {
            (*shard)[lane] = static_cast<uint8_t>(value);
            return;
        }
        (*shard)[lane] = static_cast<uint8_t>(value);
        (*shard)[32 + lane] = static_cast<uint8_t>(value >> 8);
    }

private:
    unsigned bits_;
    unsigned order_;
    uint32_t polynomial_;
    std::vector<Element> coordinate_to_polynomial_;
    std::vector<uint32_t> polynomial_to_coordinate_;
};

BinaryField make_gf8()
{
    static const Element basis[] = {
        1, 214, 152, 146, 86, 200, 88, 230
    };
    return BinaryField(8, 0x11d, basis);
}

BinaryField make_gf16()
{
    static const Element basis[] = {
        0x0001, 0xacca, 0x3c0e, 0x163e,
        0xc582, 0xed2e, 0x914c, 0x4012,
        0x6c98, 0x10d8, 0x6a72, 0xb900,
        0xfdb8, 0xfb34, 0xff38, 0x991e
    };
    return BinaryField(16, 0x1002d, basis);
}

std::vector<Element> make_independent_skew(const BinaryField& field)
{
    std::vector<Element> skew(field.order() - 1, 0);
    std::vector<Element> temp;
    for (unsigned bit = 1; bit < field.bits(); ++bit)
        temp.push_back(static_cast<Element>(1u << bit));
    for (unsigned layer = 0; layer + 1 < field.bits(); ++layer)
    {
        const unsigned step = 1u << (layer + 1);
        const unsigned first = (1u << layer) - 1;
        skew[first] = 0;
        for (unsigned index = layer; index + 1 < field.bits(); ++index)
        {
            const unsigned offset = 1u << (index + 1);
            for (unsigned position = first; position < offset; position += step)
                skew[position + offset] = skew[position] ^ temp[index];
        }
        const Element product = field.multiply(
            temp[layer], static_cast<Element>(temp[layer] ^ 1));
        const Element scale = field.inverse(product);
        for (unsigned index = layer + 1; index + 1 < field.bits(); ++index)
        {
            temp[index] = field.multiply(
                field.multiply(temp[index],
                    static_cast<Element>(temp[index] ^ 1)), scale);
        }
    }
    return skew;
}

Element production_multiplier(const BinaryField& field, unsigned index)
{
    return field.bits() == 8
        ? static_cast<Element>(leopard::ff8::TestOnlyFFTMultiplier(index))
        : static_cast<Element>(leopard::ff16::TestOnlyFFTMultiplier(index));
}

enum Direction
{
    kForward,
    kInverse
};

enum NodeKind
{
    kEmpty,
    kLeaf,
    kComplete,
    kBoundary
};

bool any(const Mask& mask)
{
    return std::find(mask.begin(), mask.end(), 1) != mask.end();
}

bool all(const Mask& mask)
{
    return std::find(mask.begin(), mask.end(), 0) == mask.end();
}

Mask slice(const Mask& mask, unsigned begin, unsigned end)
{
    return Mask(mask.begin() + begin, mask.begin() + end);
}

Mask concatenate(const Mask& left, const Mask& right)
{
    Mask result(left);
    result.insert(result.end(), right.begin(), right.end());
    return result;
}

struct Node
{
    unsigned length;
    unsigned shift;
    Direction direction;
    Mask input_mask;
    Mask output_mask;
    Mask required_inputs;
    Mask live_outputs;
    NodeKind kind;
    Element multiplier;
    std::unique_ptr<Node> left;
    std::unique_ptr<Node> right;
    Mask left_inputs;
    Mask right_inputs;
};

std::unique_ptr<Node> make_node(
    unsigned length,
    unsigned shift,
    Direction direction,
    const Mask& input_mask,
    const Mask& output_mask,
    NodeKind kind)
{
    std::unique_ptr<Node> node(new Node());
    node->length = length;
    node->shift = shift;
    node->direction = direction;
    node->input_mask = input_mask;
    node->output_mask = output_mask;
    node->required_inputs.assign(length, 0);
    node->live_outputs.assign(length, 0);
    node->kind = kind;
    node->multiplier = 0;
    return node;
}

std::unique_ptr<Node> compile_forward_node(
    const std::vector<Element>& skew,
    unsigned length,
    unsigned shift,
    const Mask& input_mask,
    const Mask& output_mask)
{
    if (!any(input_mask) || !any(output_mask))
        return make_node(length, shift, kForward, input_mask, output_mask, kEmpty);
    if (length == 1)
    {
        std::unique_ptr<Node> node = make_node(
            length, shift, kForward, input_mask, output_mask, kLeaf);
        node->required_inputs[0] = input_mask[0] && output_mask[0];
        node->live_outputs = node->required_inputs;
        return node;
    }
    if (all(input_mask) && all(output_mask))
    {
        std::unique_ptr<Node> node = make_node(
            length, shift, kForward, input_mask, output_mask, kComplete);
        node->required_inputs = input_mask;
        node->live_outputs = output_mask;
        return node;
    }

    const unsigned half = length / 2;
    const Element multiplier = skew[shift + half - 1];
    const Mask a_live = slice(input_mask, 0, half);
    const Mask b_live = slice(input_mask, half, length);
    Mask u_live(half, 0), v_live(half, 0);
    for (unsigned i = 0; i < half; ++i)
    {
        u_live[i] = a_live[i] || (b_live[i] && multiplier != 0);
        v_live[i] = a_live[i] ||
            (b_live[i] && static_cast<Element>(multiplier ^ 1) != 0);
    }
    std::unique_ptr<Node> left = compile_forward_node(
        skew, half, shift, u_live, slice(output_mask, 0, half));
    std::unique_ptr<Node> right = compile_forward_node(
        skew, half, shift + half, v_live, slice(output_mask, half, length));
    Mask required_a(half, 0), required_b(half, 0);
    for (unsigned i = 0; i < half; ++i)
    {
        const bool need_u = left->required_inputs[i] != 0;
        const bool need_v = right->required_inputs[i] != 0;
        required_a[i] = a_live[i] && (need_u || need_v);
        required_b[i] = b_live[i] &&
            ((need_u && multiplier != 0) ||
             (need_v && static_cast<Element>(multiplier ^ 1) != 0));
    }
    std::unique_ptr<Node> node = make_node(
        length, shift, kForward, input_mask, output_mask, kBoundary);
    node->multiplier = multiplier;
    node->required_inputs = concatenate(required_a, required_b);
    node->left_inputs = left->required_inputs;
    node->right_inputs = right->required_inputs;
    node->live_outputs = concatenate(left->live_outputs, right->live_outputs);
    node->left = std::move(left);
    node->right = std::move(right);
    if (!any(node->required_inputs))
    {
        node->kind = kEmpty;
        node->live_outputs.assign(length, 0);
        node->left.reset();
        node->right.reset();
    }
    return node;
}

std::unique_ptr<Node> compile_inverse_node(
    const std::vector<Element>& skew,
    unsigned length,
    unsigned shift,
    const Mask& input_mask,
    const Mask& output_mask)
{
    if (!any(input_mask) || !any(output_mask))
        return make_node(length, shift, kInverse, input_mask, output_mask, kEmpty);
    if (length == 1)
    {
        std::unique_ptr<Node> node = make_node(
            length, shift, kInverse, input_mask, output_mask, kLeaf);
        node->required_inputs[0] = input_mask[0] && output_mask[0];
        node->live_outputs = node->required_inputs;
        return node;
    }
    if (all(input_mask) && all(output_mask))
    {
        std::unique_ptr<Node> node = make_node(
            length, shift, kInverse, input_mask, output_mask, kComplete);
        node->required_inputs = input_mask;
        node->live_outputs = output_mask;
        return node;
    }

    const unsigned half = length / 2;
    const Element multiplier = skew[shift + half - 1];
    const Mask need_x = slice(output_mask, 0, half);
    const Mask need_y = slice(output_mask, half, length);
    Mask need_u(half, 0), need_v(half, 0);
    for (unsigned i = 0; i < half; ++i)
    {
        need_u[i] = need_y[i] ||
            (need_x[i] && static_cast<Element>(multiplier ^ 1) != 0);
        need_v[i] = need_y[i] || (need_x[i] && multiplier != 0);
    }
    std::unique_ptr<Node> left = compile_inverse_node(
        skew, half, shift, slice(input_mask, 0, half), need_u);
    std::unique_ptr<Node> right = compile_inverse_node(
        skew, half, shift + half, slice(input_mask, half, length), need_v);
    Mask live_x(half, 0), live_y(half, 0);
    for (unsigned i = 0; i < half; ++i)
    {
        live_x[i] = need_x[i] &&
            ((left->live_outputs[i] &&
                static_cast<Element>(multiplier ^ 1) != 0) ||
             (right->live_outputs[i] && multiplier != 0));
        live_y[i] = need_y[i] &&
            (left->live_outputs[i] || right->live_outputs[i]);
    }
    std::unique_ptr<Node> node = make_node(
        length, shift, kInverse, input_mask, output_mask, kBoundary);
    node->multiplier = multiplier;
    node->required_inputs = concatenate(
        left->required_inputs, right->required_inputs);
    node->left_inputs = left->live_outputs;
    node->right_inputs = right->live_outputs;
    node->live_outputs = concatenate(live_x, live_y);
    node->left = std::move(left);
    node->right = std::move(right);
    if (!any(node->required_inputs))
    {
        node->kind = kEmpty;
        node->live_outputs.assign(length, 0);
        node->left.reset();
        node->right.reset();
    }
    return node;
}

enum OperationKind
{
    kPairOperation,
    kCompleteOperation
};

struct VirtualOperation
{
    OperationKind kind;
    Direction direction;
    Element multiplier;
    int input_left;
    int input_right;
    int output_left;
    int output_right;
    unsigned length;
    unsigned shift;
    std::vector<int> values;
};

struct VirtualProgram
{
    int next_value;
    std::vector<int> compact_inputs;
    std::vector<int> compact_outputs;
    std::vector<VirtualOperation> operations;

    VirtualProgram() : next_value(0) {}
    int value() { return next_value++; }
};

std::vector<int> emit_virtual(
    const Node& node,
    const std::vector<int>& inputs,
    VirtualProgram* program)
{
    require(inputs.size() == node.length, "virtual input map length mismatch");
    if (node.kind == kEmpty)
        return std::vector<int>(node.length, kMissing);
    if (node.kind == kLeaf)
        return std::vector<int>(1,
            node.live_outputs[0] ? inputs[0] : kMissing);
    if (node.kind == kComplete)
    {
        VirtualOperation op;
        op.kind = kCompleteOperation;
        op.direction = node.direction;
        op.multiplier = 0;
        op.input_left = op.input_right = kMissing;
        op.output_left = op.output_right = kMissing;
        op.length = node.length;
        op.shift = node.shift;
        op.values = inputs;
        require(std::find(inputs.begin(), inputs.end(), kMissing) == inputs.end(),
            "complete subtree has a missing value");
        program->operations.push_back(op);
        return inputs;
    }

    require(node.left.get() && node.right.get(), "boundary lacks children");
    const unsigned half = node.length / 2;
    if (node.direction == kForward)
    {
        std::vector<int> left_inputs(half, kMissing);
        std::vector<int> right_inputs(half, kMissing);
        for (unsigned i = 0; i < half; ++i)
        {
            const bool need_left = node.left_inputs[i] != 0;
            const bool need_right = node.right_inputs[i] != 0;
            if (!need_left && !need_right)
                continue;
            VirtualOperation op;
            op.kind = kPairOperation;
            op.direction = kForward;
            op.multiplier = node.multiplier;
            op.input_left = node.required_inputs[i] ? inputs[i] : kMissing;
            op.input_right = node.required_inputs[half + i]
                ? inputs[half + i] : kMissing;
            op.output_left = need_left ? program->value() : kMissing;
            op.output_right = need_right ? program->value() : kMissing;
            op.length = 0;
            op.shift = 0;
            require(op.input_left != kMissing || op.input_right != kMissing,
                "forward pair has no live input");
            program->operations.push_back(op);
            left_inputs[i] = op.output_left;
            right_inputs[i] = op.output_right;
        }
        const std::vector<int> left_output = emit_virtual(
            *node.left, left_inputs, program);
        const std::vector<int> right_output = emit_virtual(
            *node.right, right_inputs, program);
        std::vector<int> result(left_output);
        result.insert(result.end(), right_output.begin(), right_output.end());
        return result;
    }

    const std::vector<int> left_input(inputs.begin(), inputs.begin() + half);
    const std::vector<int> right_input(inputs.begin() + half, inputs.end());
    const std::vector<int> left_output = emit_virtual(
        *node.left, left_input, program);
    const std::vector<int> right_output = emit_virtual(
        *node.right, right_input, program);
    std::vector<int> result(node.length, kMissing);
    for (unsigned i = 0; i < half; ++i)
    {
        const bool need_left = node.live_outputs[i] != 0;
        const bool need_right = node.live_outputs[half + i] != 0;
        if (!need_left && !need_right)
            continue;
        VirtualOperation op;
        op.kind = kPairOperation;
        op.direction = kInverse;
        op.multiplier = node.multiplier;
        op.input_left = node.left_inputs[i] ? left_output[i] : kMissing;
        op.input_right = node.right_inputs[i] ? right_output[i] : kMissing;
        op.output_left = need_left ? program->value() : kMissing;
        op.output_right = need_right ? program->value() : kMissing;
        op.length = 0;
        op.shift = 0;
        require(op.input_left != kMissing || op.input_right != kMissing,
            "inverse pair has no live input");
        program->operations.push_back(op);
        result[i] = op.output_left;
        result[half + i] = op.output_right;
    }
    return result;
}

VirtualProgram make_virtual_program(const Node& root)
{
    VirtualProgram program;
    std::vector<int> inputs(root.length, kMissing);
    for (unsigned i = 0; i < root.length; ++i)
    {
        if (!root.input_mask[i])
            continue;
        const int id = root.required_inputs[i] ? program.value() : kMissing;
        program.compact_inputs.push_back(id);
        inputs[i] = id;
    }
    const std::vector<int> outputs = emit_virtual(root, inputs, &program);
    for (unsigned i = 0; i < root.length; ++i)
        if (root.output_mask[i])
            program.compact_outputs.push_back(outputs[i]);
    return program;
}

struct PhysicalOperation
{
    OperationKind kind;
    Direction direction;
    Element multiplier;
    int input_left;
    int input_right;
    int output_left;
    int output_right;
    unsigned length;
    unsigned shift;
    std::vector<int> slots;
};

struct Program
{
    unsigned parent_size;
    unsigned shift;
    Direction direction;
    unsigned active_count;
    unsigned requested_count;
    std::vector<int> compact_input_slots;
    std::vector<int> compact_output_slots;
    std::vector<PhysicalOperation> operations;
    unsigned slot_count;
    unsigned peak_live_slots;
    unsigned maximum_complete_block;
    uint64_t serialized_schedule_bytes;
    uint64_t resident_plan_bytes_excluding_allocator_headers;

    uint64_t execution_scratch_bytes() const
    {
        return static_cast<uint64_t>(slot_count) * sizeof(Shard) +
            4 * sizeof(Shard) +
            static_cast<uint64_t>(maximum_complete_block) * sizeof(void*);
    }
};

int take_slot(std::vector<int>* free_slots, unsigned* slot_count)
{
    if (!free_slots->empty())
    {
        const int result = free_slots->back();
        free_slots->pop_back();
        return result;
    }
    return static_cast<int>((*slot_count)++);
}

Program allocate_program(
    const VirtualProgram& source,
    unsigned parent_size,
    unsigned shift,
    Direction direction)
{
    Program result;
    result.parent_size = parent_size;
    result.shift = shift;
    result.direction = direction;
    result.active_count = source.compact_inputs.size();
    result.requested_count = source.compact_outputs.size();
    result.slot_count = 0;
    result.peak_live_slots = 0;
    result.maximum_complete_block = 0;

    std::vector<unsigned> uses(source.next_value, 0);
    for (size_t op_i = 0; op_i < source.operations.size(); ++op_i)
    {
        const VirtualOperation& op = source.operations[op_i];
        if (op.kind != kPairOperation)
            continue;
        if (op.input_left != kMissing)
            ++uses[op.input_left];
        if (op.input_right != kMissing)
            ++uses[op.input_right];
    }
    for (size_t i = 0; i < source.compact_outputs.size(); ++i)
        if (source.compact_outputs[i] != kMissing)
            ++uses[source.compact_outputs[i]];

    std::vector<int> value_slots(source.next_value, kMissing);
    unsigned live = 0;
    for (size_t i = 0; i < source.compact_inputs.size(); ++i)
    {
        const int id = source.compact_inputs[i];
        if (id == kMissing)
        {
            result.compact_input_slots.push_back(kMissing);
            continue;
        }
        require(uses[id] != 0, "required input has no consumer");
        const int slot = static_cast<int>(result.slot_count++);
        value_slots[id] = slot;
        result.compact_input_slots.push_back(slot);
        ++live;
    }
    result.peak_live_slots = live;
    std::vector<int> free_slots;

    for (size_t op_i = 0; op_i < source.operations.size(); ++op_i)
    {
        const VirtualOperation& op = source.operations[op_i];
        PhysicalOperation physical;
        physical.kind = op.kind;
        physical.direction = op.direction;
        physical.multiplier = op.multiplier;
        physical.input_left = physical.input_right = kMissing;
        physical.output_left = physical.output_right = kMissing;
        physical.length = op.length;
        physical.shift = op.shift;
        if (op.kind == kCompleteOperation)
        {
            for (size_t i = 0; i < op.values.size(); ++i)
            {
                require(value_slots[op.values[i]] != kMissing,
                    "complete operation reads an unallocated value");
                physical.slots.push_back(value_slots[op.values[i]]);
            }
            result.maximum_complete_block = std::max(
                result.maximum_complete_block, op.length);
            result.operations.push_back(physical);
            continue;
        }

        if (op.input_left != kMissing)
            physical.input_left = value_slots[op.input_left];
        if (op.input_right != kMissing)
            physical.input_right = value_slots[op.input_right];
        const int input_ids[] = { op.input_left, op.input_right };
        for (unsigned side = 0; side < 2; ++side)
        {
            const int id = input_ids[side];
            if (id == kMissing)
                continue;
            require(uses[id] != 0, "value use underflow");
            if (--uses[id] == 0)
            {
                free_slots.push_back(value_slots[id]);
                value_slots[id] = kMissing;
                --live;
            }
        }
        if (op.output_left != kMissing)
        {
            physical.output_left = take_slot(&free_slots, &result.slot_count);
            value_slots[op.output_left] = physical.output_left;
            ++live;
        }
        if (op.output_right != kMissing)
        {
            physical.output_right = take_slot(&free_slots, &result.slot_count);
            value_slots[op.output_right] = physical.output_right;
            ++live;
        }
        result.peak_live_slots = std::max(result.peak_live_slots, live);
        result.operations.push_back(physical);
    }

    for (size_t i = 0; i < source.compact_outputs.size(); ++i)
    {
        const int id = source.compact_outputs[i];
        result.compact_output_slots.push_back(
            id == kMissing ? kMissing : value_slots[id]);
        if (id != kMissing)
            require(value_slots[id] != kMissing, "output value was released");
    }

    result.serialized_schedule_bytes = 32 +
        4 * (result.compact_input_slots.size() +
             result.compact_output_slots.size());
    for (size_t i = 0; i < result.operations.size(); ++i)
    {
        result.serialized_schedule_bytes += 16;
        result.serialized_schedule_bytes +=
            result.operations[i].kind == kPairOperation
                ? 6 * 4
                : 4 * result.operations[i].slots.size();
    }
    result.resident_plan_bytes_excluding_allocator_headers = sizeof(Program) +
        result.compact_input_slots.capacity() * sizeof(int) +
        result.compact_output_slots.capacity() * sizeof(int) +
        result.operations.capacity() * sizeof(PhysicalOperation);
    for (size_t i = 0; i < result.operations.size(); ++i)
        result.resident_plan_bytes_excluding_allocator_headers +=
            result.operations[i].slots.capacity() * sizeof(int);
    return result;
}

Program compile_program(
    const std::vector<Element>& skew,
    unsigned parent_size,
    unsigned shift,
    Direction direction,
    const Mask& input_mask,
    const Mask& output_mask)
{
    require(parent_size != 0 && (parent_size & (parent_size - 1)) == 0,
        "parent must be dyadic");
    require(input_mask.size() == parent_size &&
            output_mask.size() == parent_size,
        "mask size differs from parent");
    require((shift & (parent_size - 1)) == 0 &&
            shift + parent_size <= skew.size() + 1,
        "shift is not an aligned in-field coset");
    std::unique_ptr<Node> root = direction == kForward
        ? compile_forward_node(skew, parent_size, shift,
            input_mask, output_mask)
        : compile_inverse_node(skew, parent_size, shift,
            input_mask, output_mask);
    const VirtualProgram virtual_program = make_virtual_program(*root);
    return allocate_program(
        virtual_program, parent_size, shift, direction);
}

void production_transform(
    const BinaryField& field,
    Direction direction,
    unsigned length,
    unsigned shift,
    std::vector<void*>* pointers)
{
    if (pointers->size() < length)
        throw std::invalid_argument("production pointer workspace too small");
    if (field.bits() == 8)
    {
        if (direction == kForward)
            leopard::ff8::TestOnlyLchForward(
                sizeof(Shard), length, shift, length, &(*pointers)[0]);
        else
            leopard::ff8::TestOnlyLchInverse(
                sizeof(Shard), length, shift, length, &(*pointers)[0]);
    }
    else
    {
        if (direction == kForward)
            leopard::ff16::TestOnlyLchForward(
                sizeof(Shard), length, shift, length, &(*pointers)[0]);
        else
            leopard::ff16::TestOnlyLchInverse(
                sizeof(Shard), length, shift, length, &(*pointers)[0]);
    }
}

class Executor
{
public:
    Executor(const BinaryField& field, const Program& program)
        : field_(field)
        , program_(program)
        , arena_(program.slot_count)
        , pointers_(program.maximum_complete_block, NULL)
    {}

    void execute(
        const std::vector<Shard>& compact_input,
        std::vector<Shard>* compact_output)
    {
        if (compact_input.size() != program_.compact_input_slots.size())
            throw std::invalid_argument("compact input count differs from plan");
        if (compact_output->size() != program_.compact_output_slots.size())
            throw std::invalid_argument("compact output count differs from plan");
        for (size_t i = 0; i < compact_input.size(); ++i)
        {
            const int slot = program_.compact_input_slots[i];
            if (slot != kMissing)
                arena_[slot] = compact_input[i];
        }
        for (size_t op_i = 0; op_i < program_.operations.size(); ++op_i)
        {
            const PhysicalOperation& op = program_.operations[op_i];
            if (op.kind == kCompleteOperation)
            {
                if (op.slots.size() != op.length)
                    throw std::logic_error("complete operation slot count mismatch");
                for (unsigned i = 0; i < op.length; ++i)
                    pointers_[i] = arena_[op.slots[i]].data();
                production_transform(
                    field_, op.direction, op.length, op.shift, &pointers_);
                continue;
            }
            execute_pair(op);
        }
        const Shard zero = {{ 0 }};
        for (size_t i = 0; i < compact_output->size(); ++i)
        {
            const int slot = program_.compact_output_slots[i];
            (*compact_output)[i] = slot == kMissing ? zero : arena_[slot];
        }
    }

private:
    void execute_pair(const PhysicalOperation& op)
    {
        Shard left = {{ 0 }};
        Shard right = {{ 0 }};
        if (op.input_left != kMissing)
            left = arena_[op.input_left];
        if (op.input_right != kMissing)
            right = arena_[op.input_right];
        Shard output_left = {{ 0 }};
        Shard output_right = {{ 0 }};
        for (unsigned lane = 0; lane < field_.lanes(); ++lane)
        {
            const Element a = field_.get(left, lane);
            const Element b = field_.get(right, lane);
            Element x = 0, y = 0;
            if (op.direction == kForward)
            {
                x = a ^ multiply_constant(op.multiplier, b);
                y = a ^ multiply_constant(
                    static_cast<Element>(op.multiplier ^ 1), b);
            }
            else
            {
                x = multiply_constant(
                        static_cast<Element>(op.multiplier ^ 1), a) ^
                    multiply_constant(op.multiplier, b);
                y = a ^ b;
            }
            if (op.output_left != kMissing)
                field_.set(&output_left, lane, x);
            if (op.output_right != kMissing)
                field_.set(&output_right, lane, y);
        }
        if (op.output_left != kMissing)
            arena_[op.output_left] = output_left;
        if (op.output_right != kMissing)
            arena_[op.output_right] = output_right;
    }

    Element multiply_constant(Element coefficient, Element value) const
    {
        if (coefficient == 0 || value == 0)
            return 0;
        if (coefficient == 1)
            return value;
        return field_.multiply(coefficient, value);
    }

    const BinaryField& field_;
    const Program& program_;
    std::vector<Shard> arena_;
    std::vector<void*> pointers_;
};

Mask prefix_mask(unsigned length, unsigned count)
{
    require(count <= length, "prefix exceeds parent");
    Mask result(length, 0);
    std::fill(result.begin(), result.begin() + count, 1);
    return result;
}

Mask irregular_mask(unsigned length, unsigned multiplier, unsigned add, unsigned mod)
{
    Mask result(length, 0);
    for (unsigned i = 0; i < length; ++i)
        result[i] = ((i * multiplier + add) % mod) < (mod / 2);
    return result;
}

Shard deterministic_shard(
    const BinaryField& field,
    uint64_t seed,
    unsigned coordinate)
{
    Shard result = {{ 0 }};
    for (unsigned lane = 0; lane < field.lanes(); ++lane)
    {
        uint64_t value = seed ^
            (static_cast<uint64_t>(coordinate + 1) * UINT64_C(0x9e3779b97f4a7c15)) ^
            (static_cast<uint64_t>(lane + 1) * UINT64_C(0xbf58476d1ce4e5b9));
        value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
        value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
        value ^= value >> 31;
        Element symbol = static_cast<Element>(value & (field.order() - 1));
        if (field.bits() == 16)
            symbol ^= static_cast<Element>(((lane * 37u + coordinate * 13u + 1u)
                & 0xffu) << 8);
        field.set(&result, lane, symbol);
    }
    return result;
}

std::vector<unsigned> indices(const Mask& mask)
{
    std::vector<unsigned> result;
    for (unsigned i = 0; i < mask.size(); ++i)
        if (mask[i])
            result.push_back(i);
    return result;
}

struct Case
{
    std::string name;
    const BinaryField* field;
    unsigned parent_size;
    unsigned shift;
    Direction direction;
    Mask input_mask;
    Mask output_mask;
    bool benchmark;
};

struct CaseResult
{
    std::string name;
    unsigned field_bits;
    unsigned parent_size;
    unsigned shift;
    std::string direction;
    unsigned active;
    unsigned requested;
    uint64_t compared_bytes;
    unsigned pair_operations;
    unsigned complete_operations;
    unsigned maximum_complete_block;
    unsigned scratch_slots;
    unsigned peak_live_slots;
    uint64_t execution_scratch_bytes;
    uint64_t padded_scratch_bytes;
    uint64_t schedule_bytes;
    uint64_t resident_plan_bytes;
};

uint64_t fnv_update(uint64_t digest, const void* data, size_t bytes)
{
    const uint8_t* source = static_cast<const uint8_t*>(data);
    for (size_t i = 0; i < bytes; ++i)
    {
        digest ^= source[i];
        digest *= UINT64_C(1099511628211);
    }
    return digest;
}

void full_padded_execute(
    const Case& test_case,
    const std::vector<Shard>& compact_input,
    const std::vector<unsigned>& inputs,
    const std::vector<unsigned>& outputs,
    std::vector<Shard>* full_work,
    std::vector<void*>* pointers,
    std::vector<Shard>* compact_output)
{
    if (inputs.size() != compact_input.size())
        throw std::invalid_argument("full input count mismatch");
    if (outputs.size() != compact_output->size())
        throw std::invalid_argument("full output count mismatch");
    const Shard zero = {{ 0 }};
    std::fill(full_work->begin(), full_work->end(), zero);
    for (size_t i = 0; i < inputs.size(); ++i)
        (*full_work)[inputs[i]] = compact_input[i];
    for (unsigned i = 0; i < test_case.parent_size; ++i)
        (*pointers)[i] = (*full_work)[i].data();
    production_transform(*test_case.field, test_case.direction,
        test_case.parent_size, test_case.shift, pointers);
    for (size_t i = 0; i < outputs.size(); ++i)
        (*compact_output)[i] = (*full_work)[outputs[i]];
}

CaseResult verify_case(const Case& test_case, uint64_t* digest)
{
    const std::vector<Element> skew = make_independent_skew(*test_case.field);
    const Program program = compile_program(skew, test_case.parent_size,
        test_case.shift, test_case.direction,
        test_case.input_mask, test_case.output_mask);
    const std::vector<unsigned> input_indices = indices(test_case.input_mask);
    const std::vector<unsigned> output_indices = indices(test_case.output_mask);
    std::vector<Shard> input;
    for (size_t i = 0; i < input_indices.size(); ++i)
        input.push_back(deterministic_shard(*test_case.field,
            UINT64_C(0xc2c2f00d12345678) ^ test_case.shift,
            input_indices[i]));
    std::vector<Shard> actual(output_indices.size());
    std::vector<Shard> expected(output_indices.size());
    Executor executor(*test_case.field, program);
    executor.execute(input, &actual);
    std::vector<Shard> full_work(test_case.parent_size);
    std::vector<void*> pointers(test_case.parent_size, NULL);
    full_padded_execute(test_case, input, input_indices, output_indices,
        &full_work, &pointers, &expected);
    require(actual == expected, "C2 output differs from padded production: " +
        test_case.name);
    for (size_t i = 0; i < actual.size(); ++i)
        *digest = fnv_update(*digest, actual[i].data(), actual[i].size());

    CaseResult result;
    result.name = test_case.name;
    result.field_bits = test_case.field->bits();
    result.parent_size = test_case.parent_size;
    result.shift = test_case.shift;
    result.direction = test_case.direction == kForward ? "forward" : "inverse";
    result.active = input_indices.size();
    result.requested = output_indices.size();
    result.compared_bytes = static_cast<uint64_t>(actual.size()) * sizeof(Shard);
    result.pair_operations = 0;
    result.complete_operations = 0;
    for (size_t i = 0; i < program.operations.size(); ++i)
    {
        if (program.operations[i].kind == kPairOperation)
            ++result.pair_operations;
        else
            ++result.complete_operations;
    }
    result.maximum_complete_block = program.maximum_complete_block;
    result.scratch_slots = program.slot_count;
    result.peak_live_slots = program.peak_live_slots;
    result.execution_scratch_bytes = program.execution_scratch_bytes();
    result.padded_scratch_bytes =
        static_cast<uint64_t>(test_case.parent_size) * sizeof(Shard) +
        static_cast<uint64_t>(test_case.parent_size) * sizeof(void*);
    result.schedule_bytes = program.serialized_schedule_bytes;
    result.resident_plan_bytes =
        program.resident_plan_bytes_excluding_allocator_headers;
    return result;
}

double median(std::vector<double> values)
{
    require(!values.empty(), "median of empty sample");
    std::sort(values.begin(), values.end());
    const size_t middle = values.size() / 2;
    return values.size() % 2 != 0
        ? values[middle]
        : (values[middle - 1] + values[middle]) / 2.;
}

double mad(const std::vector<double>& values, double center)
{
    std::vector<double> deviations;
    for (size_t i = 0; i < values.size(); ++i)
        deviations.push_back(std::fabs(values[i] - center));
    return median(deviations);
}

struct BenchmarkResult
{
    std::string name;
    unsigned field_bits;
    unsigned parent_size;
    std::string direction;
    unsigned repetitions;
    double setup_median_us;
    double setup_mad_us;
    double candidate_median_us;
    double candidate_mad_us;
    double padded_median_us;
    double padded_mad_us;
    double speedup;
    uint64_t candidate_scratch_bytes;
    uint64_t padded_scratch_bytes;
};

typedef std::chrono::steady_clock Clock;

double elapsed_us(Clock::time_point begin, Clock::time_point end)
{
    return std::chrono::duration_cast<
        std::chrono::duration<double, std::micro> >(end - begin).count();
}

BenchmarkResult benchmark_case(const Case& test_case, unsigned scale)
{
    const std::vector<Element> skew = make_independent_skew(*test_case.field);
    unsigned repetitions = test_case.parent_size <= 256 ? 41 :
        (test_case.parent_size <= 1024 ? 23 : 11);
    repetitions = std::max(3u, repetitions * std::max(1u, scale));
    std::vector<double> setup_samples;
    Program reference;
    for (unsigned i = 0; i < repetitions; ++i)
    {
        const Clock::time_point begin = Clock::now();
        Program program = compile_program(skew, test_case.parent_size,
            test_case.shift, test_case.direction,
            test_case.input_mask, test_case.output_mask);
        const Clock::time_point end = Clock::now();
        setup_samples.push_back(elapsed_us(begin, end));
        if (i == 0)
            reference = std::move(program);
    }

    const std::vector<unsigned> input_indices = indices(test_case.input_mask);
    const std::vector<unsigned> output_indices = indices(test_case.output_mask);
    std::vector<Shard> input;
    for (size_t i = 0; i < input_indices.size(); ++i)
        input.push_back(deterministic_shard(*test_case.field,
            UINT64_C(0xbadc2ffe55aa7711), input_indices[i]));
    std::vector<Shard> candidate_output(output_indices.size());
    std::vector<Shard> padded_output(output_indices.size());
    Executor executor(*test_case.field, reference);
    std::vector<Shard> full_work(test_case.parent_size);
    std::vector<void*> pointers(test_case.parent_size, NULL);

    for (unsigned i = 0; i < 3; ++i)
    {
        executor.execute(input, &candidate_output);
        full_padded_execute(test_case, input, input_indices, output_indices,
            &full_work, &pointers, &padded_output);
    }
    require(candidate_output == padded_output,
        "benchmark warmup differs from padded transform");

    std::vector<double> candidate_samples, padded_samples;
    volatile uint64_t sink = 0;
    for (unsigned i = 0; i < repetitions; ++i)
    {
        Clock::time_point begin = Clock::now();
        executor.execute(input, &candidate_output);
        Clock::time_point end = Clock::now();
        candidate_samples.push_back(elapsed_us(begin, end));
        sink ^= candidate_output.empty() ? 0 : candidate_output[i % candidate_output.size()][i % 64];

        begin = Clock::now();
        full_padded_execute(test_case, input, input_indices, output_indices,
            &full_work, &pointers, &padded_output);
        end = Clock::now();
        padded_samples.push_back(elapsed_us(begin, end));
        sink ^= padded_output.empty() ? 0 : padded_output[i % padded_output.size()][i % 64];
    }
    (void)sink;
    require(candidate_output == padded_output,
        "benchmark result differs from padded transform");

    BenchmarkResult result;
    result.name = test_case.name;
    result.field_bits = test_case.field->bits();
    result.parent_size = test_case.parent_size;
    result.direction = test_case.direction == kForward ? "forward" : "inverse";
    result.repetitions = repetitions;
    result.setup_median_us = median(setup_samples);
    result.setup_mad_us = mad(setup_samples, result.setup_median_us);
    result.candidate_median_us = median(candidate_samples);
    result.candidate_mad_us = mad(candidate_samples, result.candidate_median_us);
    result.padded_median_us = median(padded_samples);
    result.padded_mad_us = mad(padded_samples, result.padded_median_us);
    result.speedup = result.padded_median_us / result.candidate_median_us;
    result.candidate_scratch_bytes = reference.execution_scratch_bytes();
    result.padded_scratch_bytes =
        static_cast<uint64_t>(test_case.parent_size) * sizeof(Shard) +
        static_cast<uint64_t>(test_case.parent_size) * sizeof(void*);
    return result;
}

const char* backend_name(leo2_backend backend)
{
    switch (backend)
    {
    case LEO2_BACKEND_SCALAR: return "scalar";
    case LEO2_BACKEND_SSSE3: return "ssse3";
    case LEO2_BACKEND_AVX2: return "avx2";
    case LEO2_BACKEND_NEON: return "neon";
    default: return "auto";
    }
}

std::string json_escape(const std::string& value)
{
    std::ostringstream output;
    for (size_t i = 0; i < value.size(); ++i)
    {
        const unsigned char ch = static_cast<unsigned char>(value[i]);
        if (ch == '\\' || ch == '"')
            output << '\\' << ch;
        else if (ch >= 0x20)
            output << ch;
    }
    return output.str();
}

struct RuntimeEnvironment
{
    std::string affinity_api;
    std::vector<unsigned> affinity_cpus;
    std::string omp_num_threads;
    unsigned openmp_max_threads;
    std::string hostname;
    std::string machine;
    std::string cpu_model;
};

std::string trim(const std::string& value)
{
    const size_t first = value.find_first_not_of(" \t\r\n");
    if (first == std::string::npos)
        return "";
    const size_t last = value.find_last_not_of(" \t\r\n");
    return value.substr(first, last - first + 1);
}

RuntimeEnvironment read_runtime_environment()
{
    RuntimeEnvironment result;
    const char* omp_threads = std::getenv("OMP_NUM_THREADS");
    result.omp_num_threads = omp_threads ? omp_threads : "unset";
#if defined(_OPENMP)
    result.openmp_max_threads = static_cast<unsigned>(omp_get_max_threads());
#else
    result.openmp_max_threads = 0;
#endif
#if defined(__linux__)
    cpu_set_t allowed;
    CPU_ZERO(&allowed);
    if (sched_getaffinity(0, sizeof(allowed), &allowed) == 0)
    {
        result.affinity_api = "sched_getaffinity";
        for (unsigned cpu = 0; cpu < CPU_SETSIZE; ++cpu)
            if (CPU_ISSET(cpu, &allowed))
                result.affinity_cpus.push_back(cpu);
    }
    else
    {
        result.affinity_api = "unavailable";
    }
    std::array<char, 256> hostname = {{ 0 }};
    result.hostname = gethostname(hostname.data(), hostname.size()) == 0
        ? hostname.data() : "unavailable";
    struct utsname name;
    if (uname(&name) == 0)
    {
        result.machine = std::string(name.sysname) + " " + name.release +
            " " + name.machine;
    }
    else
    {
        result.machine = "unavailable";
    }
    std::ifstream cpuinfo("/proc/cpuinfo");
    std::string line;
    while (std::getline(cpuinfo, line))
    {
        const size_t colon = line.find(':');
        if (colon != std::string::npos &&
            trim(line.substr(0, colon)) == "model name")
        {
            result.cpu_model = trim(line.substr(colon + 1));
            break;
        }
    }
#else
    result.affinity_api = "unavailable";
    result.hostname = "unavailable";
    result.machine = "unavailable";
#endif
    if (result.cpu_model.empty())
        result.cpu_model = "unavailable";
    return result;
}

void write_results(
    std::ostream& output,
    const std::string& requested_backend,
    const std::string& mode,
    leo2_backend runtime_backend,
    uint32_t context_thread_count,
    const RuntimeEnvironment& runtime_environment,
    uint64_t multiplier_checks,
    uint64_t digest,
    const std::vector<CaseResult>& cases,
    const std::vector<BenchmarkResult>& benchmarks)
{
    output << std::fixed << std::setprecision(3);
    output << "{\n"
           << "  \"schema_version\": \"leopard2-c2-cpp-v2\",\n"
           << "  \"status\": \"pass\",\n"
           << "  \"wire_identity\": \"existing padded dyadic parent\",\n"
           << "  \"requested_backend\": \"" << json_escape(requested_backend) << "\",\n"
           << "  \"runtime_backend\": \"" << backend_name(runtime_backend) << "\",\n"
           << "  \"mode\": \"" << json_escape(mode) << "\",\n"
           << "  \"context_thread_count\": " << context_thread_count << ",\n"
           << "  \"shard_bytes\": " << sizeof(Shard) << ",\n"
           << "  \"pointer_bytes\": " << sizeof(void*) << ",\n"
           << "  \"pair_staging_shards\": 4,\n"
           << "  \"runtime_environment\": {\n"
           << "    \"affinity_api\": \""
           << json_escape(runtime_environment.affinity_api) << "\",\n"
           << "    \"process_affinity_cpus\": [";
    for (size_t i = 0; i < runtime_environment.affinity_cpus.size(); ++i)
    {
        if (i != 0)
            output << ", ";
        output << runtime_environment.affinity_cpus[i];
    }
    output << "],\n"
           << "    \"omp_num_threads_env\": \""
           << json_escape(runtime_environment.omp_num_threads) << "\",\n"
           << "    \"openmp_max_threads\": "
           << runtime_environment.openmp_max_threads << ",\n"
           << "    \"hostname\": \""
           << json_escape(runtime_environment.hostname) << "\",\n"
           << "    \"machine\": \""
           << json_escape(runtime_environment.machine) << "\",\n"
           << "    \"cpu_model\": \""
           << json_escape(runtime_environment.cpu_model) << "\"\n"
           << "  },\n"
           << "  \"build_binding\": {\n"
           << "    \"source_sha256\": \"" << LEO2_C2_SOURCE_SHA256 << "\",\n"
           << "    \"core_git_sha\": \"" << LEO2_C2_CORE_GIT_SHA << "\",\n"
           << "    \"core_matrix_sha256\": \""
           << LEO2_C2_CORE_MATRIX_SHA256 << "\",\n"
           << "    \"linked_library_sha256\": \""
           << LEO2_C2_LIBRARY_SHA256 << "\",\n"
           << "    \"sanitizer_mode\": \"" << LEO2_C2_SANITIZER_MODE << "\",\n"
           << "    \"compiler\": \"" << json_escape(__VERSION__) << "\"\n"
           << "  },\n"
           << "  \"production_multiplier_checks\": " << multiplier_checks << ",\n"
           << "  \"correctness_digest_fnv1a64\": \"0x"
           << std::hex << std::setw(16) << std::setfill('0') << digest
           << std::dec << std::setfill(' ') << "\",\n"
           << "  \"cases\": [\n";
    for (size_t i = 0; i < cases.size(); ++i)
    {
        const CaseResult& item = cases[i];
        output << "    {\"name\": \"" << json_escape(item.name)
               << "\", \"field_bits\": " << item.field_bits
               << ", \"parent_size\": " << item.parent_size
               << ", \"shift\": " << item.shift
               << ", \"direction\": \"" << item.direction
               << "\", \"active\": " << item.active
               << ", \"requested\": " << item.requested
               << ", \"compared_bytes\": " << item.compared_bytes
               << ", \"pair_operations\": " << item.pair_operations
               << ", \"complete_operations\": " << item.complete_operations
               << ", \"maximum_complete_block\": " << item.maximum_complete_block
               << ", \"scratch_slots\": " << item.scratch_slots
               << ", \"peak_live_slots\": " << item.peak_live_slots
               << ", \"execution_scratch_bytes\": " << item.execution_scratch_bytes
               << ", \"padded_scratch_bytes\": " << item.padded_scratch_bytes
               << ", \"serialized_schedule_bytes\": " << item.schedule_bytes
               << ", \"resident_plan_bytes_excluding_allocator_headers\": "
               << item.resident_plan_bytes << "}"
               << (i + 1 == cases.size() ? "\n" : ",\n");
    }
    output << "  ],\n  \"benchmarks\": [\n";
    for (size_t i = 0; i < benchmarks.size(); ++i)
    {
        const BenchmarkResult& item = benchmarks[i];
        output << "    {\"name\": \"" << json_escape(item.name)
               << "\", \"field_bits\": " << item.field_bits
               << ", \"parent_size\": " << item.parent_size
               << ", \"direction\": \"" << item.direction
               << "\", \"repetitions\": " << item.repetitions
               << ", \"setup_median_us\": " << item.setup_median_us
               << ", \"setup_mad_us\": " << item.setup_mad_us
               << ", \"candidate_median_us\": " << item.candidate_median_us
               << ", \"candidate_mad_us\": " << item.candidate_mad_us
               << ", \"padded_median_us\": " << item.padded_median_us
               << ", \"padded_mad_us\": " << item.padded_mad_us
               << ", \"padded_over_candidate\": " << item.speedup
               << ", \"candidate_scratch_bytes\": " << item.candidate_scratch_bytes
               << ", \"padded_scratch_bytes\": " << item.padded_scratch_bytes
               << "}" << (i + 1 == benchmarks.size() ? "\n" : ",\n");
    }
    output << "  ]\n}\n";
}

struct Arguments
{
    std::string output_path;
    std::string backend_label;
    std::string mode;
    unsigned benchmark_scale;

    Arguments() : backend_label("unspecified"), mode("all"), benchmark_scale(1) {}
};

Arguments parse_arguments(int argc, char** argv)
{
    Arguments result;
    for (int i = 1; i < argc; ++i)
    {
        const std::string argument(argv[i]);
        if ((argument == "--output" || argument == "--backend-label" ||
             argument == "--mode" || argument == "--benchmark-scale") &&
            i + 1 >= argc)
            throw std::invalid_argument(argument + " requires a value");
        if (argument == "--output")
            result.output_path = argv[++i];
        else if (argument == "--backend-label")
            result.backend_label = argv[++i];
        else if (argument == "--mode")
            result.mode = argv[++i];
        else if (argument == "--benchmark-scale")
            result.benchmark_scale = static_cast<unsigned>(std::strtoul(argv[++i], NULL, 10));
        else
            throw std::invalid_argument("unknown argument: " + argument);
    }
    require(result.mode == "all" || result.mode == "correctness",
        "mode must be all or correctness");
    require(result.benchmark_scale >= 1, "benchmark scale must be positive");
    return result;
}

void add_prefix_cases(
    std::vector<Case>* cases,
    const std::string& prefix,
    const BinaryField& field,
    unsigned parent,
    unsigned active,
    unsigned requested,
    unsigned shift,
    bool benchmark)
{
    for (unsigned direction = 0; direction < 2; ++direction)
    {
        Case test_case;
        std::ostringstream name;
        name << prefix << "_n" << parent << "_s" << shift << "_"
             << (direction == 0 ? "forward" : "inverse");
        test_case.name = name.str();
        test_case.field = &field;
        test_case.parent_size = parent;
        test_case.shift = shift;
        test_case.direction = direction == 0 ? kForward : kInverse;
        test_case.input_mask = prefix_mask(parent, active);
        test_case.output_mask = prefix_mask(parent, requested);
        test_case.benchmark = benchmark;
        cases->push_back(test_case);
    }
}

} // namespace

int main(int argc, char** argv)
{
    leo2_context* context = NULL;
    try
    {
        const Arguments arguments = parse_arguments(argc, argv);
        leo2_context_options options;
        std::memset(&options, 0, sizeof(options));
        options.struct_size = sizeof(options);
        options.backend = LEO2_BACKEND_AUTO;
        options.thread_count = 1;
        const leo2_result create_result = leo2_context_create(&options, &context);
        require(create_result == LEO2_SUCCESS,
            std::string("context create failed: ") + leo2_result_string(create_result));
        const leo2_backend runtime_backend = leo2_context_backend(context);
        const RuntimeEnvironment runtime_environment = read_runtime_environment();

        const BinaryField gf8 = make_gf8();
        const BinaryField gf16 = make_gf16();
        const std::vector<Element> skew8 = make_independent_skew(gf8);
        const std::vector<Element> skew16 = make_independent_skew(gf16);
        uint64_t multiplier_checks = 0;
        for (unsigned i = 0; i < skew8.size(); ++i)
        {
            require(skew8[i] == production_multiplier(gf8, i),
                "independent GF8 skew differs from production");
            ++multiplier_checks;
        }
        for (unsigned i = 0; i < skew16.size(); ++i)
        {
            require(skew16[i] == production_multiplier(gf16, i),
                "independent GF16 skew differs from production");
            ++multiplier_checks;
        }

        std::vector<Case> cases;
        add_prefix_cases(&cases, "gf8_prefix", gf8, 128, 65, 33, 0, false);
        add_prefix_cases(&cases, "gf8_last_coset", gf8, 128, 65, 33, 128, false);
        add_prefix_cases(&cases, "gf8_benchmark", gf8, 256, 129, 65, 0, true);

        add_prefix_cases(&cases, "gf16_small", gf16, 256, 129, 65, 0, false);
        add_prefix_cases(&cases, "gf16_small_last", gf16, 256, 129, 65, 65280, false);
        add_prefix_cases(&cases, "gf16_medium", gf16, 1024, 513, 257, 0, true);
        add_prefix_cases(&cases, "gf16_medium_last", gf16, 1024, 513, 257, 64512, false);
        add_prefix_cases(&cases, "gf16_large", gf16, 4096, 2049, 1025, 0, true);
        add_prefix_cases(&cases, "gf16_large_last", gf16, 4096, 2049, 1025, 61440, false);
        add_prefix_cases(&cases, "gf16_deep", gf16, 8192, 4097, 2049, 57344, false);

        for (unsigned direction = 0; direction < 2; ++direction)
        {
            Case irregular;
            irregular.name = std::string("gf16_irregular_") +
                (direction == 0 ? "forward" : "inverse");
            irregular.field = &gf16;
            irregular.parent_size = 1024;
            irregular.shift = 32768;
            irregular.direction = direction == 0 ? kForward : kInverse;
            irregular.input_mask = irregular_mask(1024, 5, 3, 11);
            irregular.output_mask = irregular_mask(1024, 7, 1, 13);
            irregular.benchmark = false;
            cases.push_back(irregular);
        }

        uint64_t digest = UINT64_C(14695981039346656037);
        std::vector<CaseResult> case_results;
        std::vector<BenchmarkResult> benchmark_results;
        for (size_t i = 0; i < cases.size(); ++i)
        {
            case_results.push_back(verify_case(cases[i], &digest));
            if (arguments.mode == "all" && cases[i].benchmark)
                benchmark_results.push_back(benchmark_case(
                    cases[i], arguments.benchmark_scale));
        }

        std::unique_ptr<std::ofstream> file;
        std::ostream* output = &std::cout;
        if (!arguments.output_path.empty())
        {
            file.reset(new std::ofstream(arguments.output_path.c_str(),
                std::ios::out | std::ios::trunc));
            require(file->good(), "could not open result path");
            output = file.get();
        }
        write_results(*output, arguments.backend_label, arguments.mode,
            runtime_backend, leo2_context_thread_count(context),
            runtime_environment,
            multiplier_checks, digest,
            case_results, benchmark_results);
        leo2_context_destroy(context);
        return 0;
    }
    catch (const std::exception& error)
    {
        if (context)
            leo2_context_destroy(context);
        std::cerr << "FAIL c2_truncated_cpp: " << error.what() << std::endl;
        return 1;
    }
}

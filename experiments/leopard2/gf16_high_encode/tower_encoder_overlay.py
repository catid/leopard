#!/usr/bin/env python3
"""Strict isolated field overlay; never changes the production worktree file."""
import hashlib

BASE = 'fd27e72cb09068d83bc32b23acf49c6d133382579f109d6cfa227361c43d8cbf'


def overlay(source):
    if hashlib.sha256(source.encode()).hexdigest() != BASE:
        raise ValueError('unqualified field base')
    begin = source.index('static void IFFT_DIT_Encoder_Impl(')
    end = source.index('\nstatic void IFFT_DIT_Encoder(', begin)
    encoder = source[begin:end]
    old = 'memcpy(work[i], data[i], bytes);'
    if encoder.count(old) != 2:
        raise ValueError('source boundary anchors')
    encoder = encoder.replace(old, 'tower_encoder::CopySource(ops, work[i], data[i], bytes);')
    source = source[:begin] + encoder + source[end:]
    begin = source.index('void ReedSolomonEncodeWithSourcePolicy(')
    end = source.index('\n\nvoid ReedSolomonEncode(', begin)
    entry = source[begin:end]
    entry = entry.replace('const backend::Ops& ops,', 'const backend::Ops& original_ops,', 1)
    anchor = '{\n#if !defined(LEO2_ENABLE_TEST_HOOKS)'
    if entry.count(anchor) != 1 or not entry.endswith('\n}\n'):
        raise ValueError('field entry anchors')
    entry = entry.replace(anchor, '''{
    const bool use_tower = tower_encoder::Select(original_ops, buffer_bytes,
        source_policy_bytes, m, sparse_plans);
    const backend::Ops& ops = use_tower
        ? tower_encoder::GetOps(original_ops) : original_ops;
#if !defined(LEO2_ENABLE_TEST_HOOKS)''', 1)
    entry = entry[:-2] + '''    if (use_tower)
        tower_encoder::Finish(work, recovery_count, buffer_bytes);
}
'''
    return '#include "tower_encoder.h"\n' + source[:begin] + entry + source[end:]

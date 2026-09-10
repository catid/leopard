#!/usr/bin/env python3
"""Independent retained-record oracle for the three-epoch clock-free frontend.

Never launches a codec, imports a collector, or reads benchmark clocks. Reuse
only the qualified one-epoch identity/ELF/geometry primitives; cumulative
public calls, epoch marks, arithmetic-fault boundaries and inventories are
reconstructed here.
"""
import ast
import copy
import json
from pathlib import Path
import sys

import paired_epoch_overlay as overlay
import verify_paired_metadata as prior
from verify_paired_r19932 import ARCHIVES, CELLS, expected, valid
from verify_paired_metadata import (sha, compare, parse, keys, integer, equal,
    require, elf, split_scope, scope_command, check_inventory, build_resources)

BEAD = 'leopard-79h.38.5.4.19.1.4.3'
REFERENCE_PIN = 'f8f2529034e462ccbcfcc94ec3951182a3497a34e1258e01b32b49e3354ae543'
REFERENCE_BUILD_ROOT = '/tmp/leopard-paired-metadata.MnmrEa/build'
DRIVER = 'leopard-paired-epoch-r19932/v1'
META = 'paired-epoch-runtime-metadata/v1'
MARKS = 'paired-epoch-public-marks/v1'
CLOCK = 'paired-epoch-synthetic-clock/v1'
PROGRESS = 'paired-epoch-progress/v1'
ACCOUNT = 'paired-epoch-accounting/v1'
WITNESS = 'leopard-paired-witness/v1'
BASE_NAMES = ('paired_timer_r19932.cpp', 'PairedRuntimeMetadata.h', 'paired_timer_clock.cpp',
              'paired_timer_witness.cpp', 'paired_public_witness.cpp', 'PairedGroupTiming.h',
              'clock_guard.cpp', 'tower_public_scope.sh')
BUILD_TOOLS = frozenset(('build_paired_epoch.py','paired_epoch_overlay.py','paired_metadata_overlay.py',
    'verify_auto_gfni_boundary_checks.py','verify_auto_r19932_checks.py','verify_paired_metadata.py',
    'verify_paired_r19932.py','verify_paired_timer_r19932.py'))
RUN_TOOLS = BUILD_TOOLS | {'run_paired_epoch_checks.py','verify_paired_epoch.py','verify_paired_epoch_units.py'}


def inventory():
    rows = []
    for old in prior.inventory():
        for epoch in (range(3) if old['fault'] else (None,)):
            row = dict(old, fault_epoch=epoch)
            if epoch is not None: row['label'] += '-epoch' + str(epoch)
            rows.append(row)
    equal(len(rows), 366)
    equal(len({r['label'] for r in rows}), len(rows))
    return rows


def public_oracle(profile, cell, schedule, group, exercise, completed=3, tail=None):
    """Ordered calls and boundary marks, including partial clock-failure epochs."""
    valid(profile, cell, schedule, group)
    integer(completed, 0, 3)
    require(tail in (None, 'abort', 'fault'), 'public oracle tail')
    require(tail is None or (completed < 3 and exercise), 'partial epoch domain')
    api = 2 if profile == 'native' else 1 if cell == 1 else 0
    states = [2 if s == 'N' else int(s) for s in schedule]
    counts, apis, calls, digest, marks = [0]*3, [0]*3, 0, 14695981039346656037, []
    def observe(state):
        nonlocal calls, digest
        calls += 1; counts[state] += 1; apis[api] += 1
        digest = ((digest ^ (state + 4*api)) * 1099511628211) & ((1 << 64)-1)
    def mark():
        marks.append(dict(endpoint=len(marks), calls=calls, states=list(counts),
                          apis=list(apis), order_hash=f'{digest:016x}'))
    for epoch in range(completed + (tail is not None)):
        mark()
        for state in states: observe(state)
        partial = epoch == completed
        passes = 4 if partial else 25 if exercise else 0
        for _ in range(passes):
            for state in states:
                for _ in range(group): observe(state)
        if partial:
            if tail == 'fault':
                for _ in range(group): observe(states[0])
        else: mark()
    return (dict(schema=WITNESS, calls=calls, states=counts, apis=apis, order_hash=f'{digest:016x}'),
            dict(schema=MARKS, marks=marks, timed=False))


def clocks(group, count=504):
    require(type(group) is int and group in (1,256), 'clock group')
    integer(count, 0, 504)
    calls = []
    for endpoint in range(count):
        epoch, local = divmod(endpoint, 168)
        span, end = divmod(local, 2)
        calls.append(epoch*(4+100*group) + 4+16*group + (span+end)*group)
    return dict(schema=CLOCK, clock_calls=count, public_calls_at_clock=calls, timed=False)


def validate_metadata(meta, profile, cell, schedule, exercise, image):
    keys(meta, 'schema bead timed observation pointer_bytes preflight_gfni_counts selections native_route_label snapshots')
    equal([meta['schema'], meta['bead'], meta['timed'], meta['observation']],
          [META, BEAD, False, 'new_three_epoch_endpoints'])
    stride = 104 if exercise else 4
    require(type(meta['selections']) is list and len(meta['selections']) == 3*stride, 'epoch selections')
    require(type(meta['snapshots']) is list and len(meta['snapshots']) == 6, 'six snapshots')
    require(type(meta['preflight_gfni_counts']) is list and len(meta['preflight_gfni_counts']) == 3, 'epoch probes')
    origin = None
    for epoch in range(3):
        normal = copy.deepcopy(meta)
        normal.update(schema=prior.META, bead=prior.BEAD, observation='new_frontend_endpoints',
                      preflight_gfni_counts=meta['preflight_gfni_counts'][epoch], selections=[], snapshots=[])
        for local, record in enumerate(meta['selections'][epoch*stride:(epoch+1)*stride]):
            equal(record['epoch'], epoch); equal(record['index'], epoch*stride+local)
            selection = dict(record); selection.pop('epoch'); selection['index'] = local
            normal['selections'].append(selection)
        for phase in range(2):
            record = meta['snapshots'][epoch*2+phase]
            equal([record['endpoint'], record['epoch'], record['phase']],
                  [epoch*2+phase, epoch, 'after' if phase else 'before'])
            snapshot = {k:v for k,v in record.items() if k not in ('endpoint','epoch','phase')}
            if origin is None: origin = snapshot
            equal(snapshot, origin)  # All six endpoints, NOT only three local pairs.
            normal['snapshots'].append(dict(snapshot, endpoint=phase))
        prior.validate_metadata(normal, profile, cell, schedule, exercise, image)


def records(output):
    result = {}
    for line in output.splitlines():
        row = parse(line)
        require(type(row) is dict and type(row.get('schema')) is str, 'schema record')
        require(row['schema'] not in result, 'duplicate schema')
        result[row['schema']] = row
    return result


def verify_record(row, output, error, image):
    actual = records(output)
    p, kind, args = row['profile'], row['binary'], row['args']
    progress = dict(schema=PROGRESS, selection_count=0, snapshot_count=0, timed=False)
    if '-bad-' in row['label']:
        w, m = public_oracle(p, 0, 'NNNN' if p=='native' else '0110', 1, False, completed=0)
        wanted = {WITNESS:w, MARKS:m, PROGRESS:progress}
        if kind == 'synthetic': wanted[CLOCK] = clocks(1, 0)
        equal(actual, wanted)
        require(bool(error), 'CLI refusal diagnostic')
        if args and args[0]=='--measure': equal(error, 'metadata refuses real timing\n')
        return
    mode, cell, schedule, group = args
    cell, group = int(cell), int(group)
    exercise = mode != '--check'
    wanted = {PROGRESS:progress}
    if row['code'] == 0:
        equal(error, '')
        base = expected(p, cell, schedule, group, exercise)
        per_epoch = dict(base)
        base.update(schema=DRIVER, clock_source=kind, encode_calls=3*base['encode_calls'],
                    selections=3*base['selections'], per_slot_calls=[3*n for n in base['per_slot_calls']],
                    samples=[[257+17*i, (257+17*i)/group] for i in range(252)] if kind=='synthetic' else [])
        for key in ('warmup_passes', 'exercise_passes'): base[key+'_per_epoch'] = base.pop(key)
        wanted[DRIVER] = base
        wanted[ACCOUNT] = dict(schema=ACCOUNT, timed=False, epochs=[dict(epoch=e,
            calls=per_epoch['encode_calls'], selections=per_epoch['selections'],
            per_slot_calls=per_epoch['per_slot_calls'], probes=per_epoch['probes'],
            sample_begin=e*84 if kind=='synthetic' else 0, sample_count=84 if kind=='synthetic' else 0)
            for e in range(3)])
        w, m = public_oracle(p, cell, schedule, group, exercise)
        progress.update(selection_count=base['selections'], snapshot_count=6)
        require(META in actual, 'complete metadata missing')
        validate_metadata(actual.pop(META), p, cell, schedule, exercise, image)
        if kind == 'synthetic': wanted[CLOCK] = clocks(group)
    elif row['code'] == 86:
        equal(error, 'unexpected driver benchmark clock\n')
        w, m = public_oracle(p, cell, schedule, group, True, completed=0, tail='abort')
        progress.update(selection_count=21, snapshot_count=1)
    else:
        equal(row['code'], 1); integer(row['fault_epoch'], 0, 2)
        equal(error, ('group duration exceeds exact binary64 integer range' if row['fault']=='huge'
                      else 'nonpositive or reversed grouped clock interval')+'\n')
        epoch = row['fault_epoch']
        w, m = public_oracle(p, cell, schedule, group, True, completed=epoch, tail='fault')
        progress.update(selection_count=104*epoch+21, snapshot_count=2*epoch+1)
        wanted[CLOCK] = clocks(group, 168*epoch+2)
    wanted.update({WITNESS:w, MARKS:m})
    equal(actual, wanted)  # No success schema is accepted after any failure.


def tool_closure(folder, roots):
    seen = set()
    def visit(name):
        if name in seen: return
        seen.add(name)
        for node in ast.walk(ast.parse((folder/name).read_text())):
            modules = ([node.module] if isinstance(node, ast.ImportFrom) else
                       [alias.name for alias in node.names] if isinstance(node, ast.Import) else [])
            for module in modules:
                if module and (folder/(module+'.py')).is_file(): visit(module+'.py')
    for root in roots: visit(root)
    return seen


def verify_tools(folder, pins, required, roots):
    # Discovery alone treats a removed local import as an external dependency.
    # Require the authoritative complete names even if both file and digest
    # entry were removed; then independently check transitive import coverage.
    equal(sorted(pins), sorted(required))
    check_inventory(folder, pins)
    equal(sorted(tool_closure(folder, roots)), sorted(pins))


def verify_build(root):
    folder = root/'build'; build = parse((folder/'build.json').read_text())
    equal([build['bead'], build['completed'], build['timed'], build['reference_build_sha256']],
          [BEAD, True, False, REFERENCE_PIN])
    equal(sha(folder/'reference_build.json'), REFERENCE_PIN)
    reference = parse((folder/'reference_build.json').read_text())
    # This pinned older manifest records full commands but has no root field.
    original = REFERENCE_BUILD_ROOT
    current = str(Path(build['root'])/'build')
    wanted = [[a.replace(original+'/', current+'/') for a in command] for command in reference['commands']]
    equal(build['commands'], wanted)
    baseline = {name:reference['artifacts'][name] for name in BASE_NAMES}
    for name, digest in reference['artifacts'].items():
        if any(name == p+'/codec.a' or name.startswith(p+'/include/') for p in ARCHIVES): baseline[name] = digest
    equal(build['baseline'], baseline)
    for name, digest in baseline.items():
        equal(sha(folder/('baseline/'+name if name in BASE_NAMES else name)), digest)
    equal((folder/'paired_metadata.cpp').read_text(), overlay.driver((folder/'baseline/paired_timer_r19932.cpp').read_text()))
    adapters = {'PairedRuntimeMetadata.h':overlay.header, 'paired_timer_clock.cpp':overlay.clock,
                'paired_timer_witness.cpp':overlay.witness}
    for name in BASE_NAMES:
        if name in adapters: equal((folder/name).read_text(), adapters[name]((folder/'baseline'/name).read_text()))
        else: equal(sha(folder/name), baseline[name])
    verify_tools(folder/'source_tools', build['tools'], BUILD_TOOLS, ['build_paired_epoch.py'])
    names = set(BASE_NAMES) | {'baseline/'+n for n in BASE_NAMES} | {'reference_build.json', 'paired_metadata.cpp'}
    names |= {n for n in baseline if n not in BASE_NAMES} | {'source_tools/'+n for n in build['tools']}
    images = {}
    for p, digest in ARCHIVES.items():
        equal(sha(folder/p/'codec.a'), digest)
        names |= {p+'/'+n for n in ('driver.o','witness.o','abort-clock.o','synthetic-clock.o','abort','synthetic')}
        for kind in ('abort', 'synthetic'): images[p,kind] = elf(folder/p/kind, p=='native')
    equal(sorted(build['artifacts']), sorted(names))
    check_inventory(folder, build['artifacts'], exclude=('build.json',))
    return build, images, build_resources((root/'build.log').read_text())


def replay(root):
    build, images, build_peak = verify_build(root)
    checks = parse((root/'checks/checks.json').read_text())
    equal([checks['bead'], checks['completed'], checks['timed'], checks['root']], [BEAD, True, False, build['root']])
    verify_tools(root/'tools', checks['tools'], RUN_TOOLS, ['run_paired_epoch_checks.py'])
    rows = inventory()
    equal([r['label'] for r in checks['records']], [r['label'] for r in rows])
    totals = dict(positive=0, clock_aborts=0, clock_faults=0, cli_refusals=0, parity_comparisons=0,
                  parity_bytes=0, selections=0, snapshots=0, public_encode_calls=0, synthetic_spans=0)
    peaks, files = [], {'checks.json'}
    for record, row in zip(checks['records'], rows):
        label, p, kind = row['label'], row['profile'], row['binary']
        equal([record['returncode'], record['fault'], record['fault_epoch']], [row['code'], row['fault'], row['fault_epoch']])
        args = list(row['args']); parity = root/'checks'/(label+'.parity') if row['parity'] else None
        if parity: args.append(str(Path(checks['root'])/'checks'/parity.name)); files.add(parity.name)
        equal(record['args'], [p, kind, *args])
        equal(record['scope_command'], scope_command(checks['root'], p, kind, args))
        stdout, stderr = root/'checks'/(label+'.stdout'), root/'checks'/(label+'.stderr')
        files |= {stdout.name, stderr.name}
        equal(sha(stdout), record['stdout_sha256']); equal(sha(stderr), record['stderr_sha256'])
        output, error, peak = split_scope(stdout.read_text(), stderr.read_text(), row['code'])
        equal(record['memory_peak'], peak); peaks.append(peak)
        verify_record(row, output, error, images[p,kind])
        if parity:
            totals['positive'] += 1; totals['snapshots'] += 6
            c, g = int(row['args'][1]), int(row['args'][3]); exercise = row['args'][0]!='--check'
            totals['selections'] += 312 if exercise else 12
            totals['public_encode_calls'] += 3*expected(p,c,row['args'][2],g,exercise)['encode_calls']
            totals['synthetic_spans'] += 252 if kind=='synthetic' else 0
            equal(sha(parity), record['parity_sha256']); equal(parity.stat().st_size, CELLS[c][1]*CELLS[c][2])
            baseline = root/'checks'/f'native-{c}-NNNN-1-check.parity'
            if parity != baseline:
                totals['parity_bytes'] += compare(parity, baseline); totals['parity_comparisons'] += 1
        else: totals['cli_refusals' if '-bad-' in label else 'clock_aborts' if row['code']==86 else 'clock_faults'] += 1
    equal(sorted(p.name for p in (root/'checks').iterdir()), sorted(files))
    from verify_paired_epoch_units import replay as unit_replay
    unit_result = unit_replay(root)
    return dict(bead=BEAD, timed=False, default_enabled=False, totals=totals, build_peak=build_peak,
                maximum_native_peak=max(peaks), all_native_memory_events_zero=True, swap_bytes=0,
                boundary_units=unit_result, qualification_complete=False,
                remaining='real-record adversarial checks, final review and sealed delivery',
                historical_shift_cause_established=False)


if __name__ == '__main__': print(json.dumps(replay(Path(sys.argv[1]).resolve(strict=True)), sort_keys=True))

#!/usr/bin/env python3
"""Single preregistered paired successor; leopard-79h.38.5.4.18.4.3.

Never use this collector for dry runs: creating attempt1 consumes the attempt.
Pure analysis and validation functions do not execute binaries or read clocks.
"""
import fcntl
import json
import math
import os
from pathlib import Path
import statistics
import subprocess
import sys

from run_split_cache_screen import check_passive, digest, host_identity, sibling_ticks
from verify_paired_r19932 import equal, parse, require
from verify_tower_public import expected

PLAN = 'tower_screen_plan.json'
HOST = dict(hostname='work', kernel='6.8.0-137-generic', vendor_id='AuthenticAMD',
            **{'cpu family':'26', 'model':'8', 'model name':'AMD Ryzen Threadripper 9980X 64-Cores'})
PAIRS = ('0110', '1001')
TARGETS = (0,1,2,4,8)
AFFECTED = ()
UNCHANGED = (3,5,6,7)
CROSS = {'same_off': ['0000']*4, 'same_on': ['1111']*4,
         'original_on': ['PPPP','1111','1111','PPPP'],
         'original_off': ['PPPP','0000','0000','PPPP'],
         'same_original': ['PPPP']*4}
NATIVE_CROSS = {'native_on': ['NNNN','1111','1111','NNNN'], 'same_native':['NNNN']*4}
PINS = {
    "native": "9f777340c054041870f15c916374ae5728ea3c4cb6a5291c12ffe1cd0e942115",
    "original": "dd0c82b408e00744716f809c1c87a7f972eb5d7be010079b5154c44cbd290f68",
    "current": "45ba62393647dcc64873eb592402e551b33de37719245540310ca0bf988920b9",
    "native.a": "3f83c55599cbfbeea5d15adf3627405bcde1cbed2ddf4a672fbc75b186f5c4c1",
    "original.a": "89f33d3d43f792de078c38469f9406dded2b69126af2c6fc8db6e2acb77cd334",
    "current.a": "24a78109405a25a5164e656bcaabbcc9aff586849a7899492cbaf290b3018ba9",
    "build.json": "d304b72038eb460efd69ccd2526bd26c187fec083be7ba2f2724c7085a479818",
    "checks.json": "71c9d4d37e96322f5381558a452c3fd7e51da24a7b4f9ac9c3066dc7f0e7320b",
    "qualification.log": "e796f82d638e888d142fc8ff7335d97711299471973b79de00eb90c64cf200ec",
    "qualification.manifest": "c52a32f869cdd5c8aa4886fcc9f4bc2ee77fa7f1f9a2bc27cd295600b29a1c20",
    "l2-driver.o": "99223acc12f2f616c0e4e5c0cb39d7a1a5207aa27570eedf3f08a1ec9290aa75",
    "native-driver.o": "e718e89c0d9b454862251034adef1dd598e21d242959af5bdc1941aaa0d93ba4",
    "leopard2.cpp": "93412028272454b2457f9928a6bc4c34684591f5bee2c882a3ad0000fd50cc70",
    "Leopard2Direct.h": "da4b39084cf7be12e824f46e799386d742120fe50a852c8b6a12dd9fe0c6b423"
}
PROTOCOL = dict(
    schema='leopard-tower-plan/v1', bead='leopard-79h.38.5.4.18.4.3',
    host=HOST, cpu=26, sibling=90, controller_cpu=0, passive_seconds=10,
    condition='slipgate-disabled-20260909', attempt_budget=1,
    attempt_root='/tmp/leopard-tower-screen.0mjTlC/attempt1',
    rounds=3, sample_passes=21, spans_per_process=84, groups=[1]*7+[256,1],
    paired_orders=list(PAIRS), cross_orders=CROSS, native_cross_orders=NATIVE_CROSS,
    target_cells=list(TARGETS), affected_neighbors=list(AFFECTED), unchanged_neighbors=list(UNCHANGED),
    native_cells=list(TARGETS), cells=[[1000,200,32768,2,3,"encode","target"],[1000,199,65536,2,3,"encode","target"],[1000,200,65536,2,3,"encode","target"],[4096,512,4096,2,3,"encode","unchanged_neighbor"],[1000,199,32768,2,0,"encode","target"],[1000,200,32768,2,6,"encode","unchanged_neighbor"],[1000,200,32768,2,0,"encode","unchanged_neighbor"],[17,7,64,1,3,"encode","unchanged_neighbor"],[1000,200,32768,2,3,"one_item_batch","target"]],
    timed_processes=714, preflights=36, cross_control_aggregates=32, within_control_aggregates=128,
    paired_estimator='median21_sqrt_off_product_over_on_product',
    within_estimator='median21_sqrt_outer_product_over_inner_product_each_samepath_process_slot',
    process_estimator='median84_all_group_averages', cross_estimator='sqrt_cost0_cost3_over_cost1_cost2',
    round_aggregation='geometric_mean_three_rounds_each_comparison_separately',
    minimum_gain=1.10, original_minimum_gain=1.10, native_minimum_gain=1.10, equivalence_bound=1.02,
    affected_neighbor_minimum=1/1.02, original_off_diagnostic_only=True,
    every_target_round_positive=True, production_promotion=False, confidence_intervals=False,
    authoritative_v19=False, pristine_off_claim=False,
    temperature='warm_after_four_preflights', extra_table_bytes=6291456,
    boundary_conversions_inside_encode=True, cold_initialization_measured=False,
    codec_commit='45e2effd869859c9b3aa48190eff6f4738817c61',
    native_commit='6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198',
    runtime_commit='6917d49425841d656416820e13ff2b77e6bdc507',
    timer_commit='7d17cda578b0a5e991d10f0c956fa1803bd09f62', artifact_sha256=PINS)
SOURCE_FILES = {
    'verify_tower_public.py',
    'build_tower_public.py',
    'tower_public_overlay.py',
    'tower_public_link.h',
    'tower_public_link.cpp',
    'tower_encoder.h',
    'tower_encoder.cpp',
    'tower_encoder_overlay.py',
    'tower_butterfly_probe.h',
    'tower_butterfly_probe.cpp',
    'tower_public_scope.sh',
    'run_tower_screen.py',
    'replay_tower_screen.py',
    'test_tower_screen.py',
    'freeze_tower_screen.py',
    'tower_screen_method.md',
    'run_split_cache_screen.py',
    'verify_paired_r19932.py',
    'verify_auto_r19932_checks.py',
    'verify_auto_gfni_boundary_checks.py',
    'verify_avx2_adjacent_public.py',
    'build_avx2_adjacent_public.py',
    'avx2_adjacent_counts.py',
    'verify_gf16_callback_probe.py',
    'verify_gfni_source_stage.py',
    'avx2_adjacent_public.cpp',
    'avx2_adjacent_public_link.cpp',
    'avx2_adjacent_public_link.h',
    'avx2_adjacent_control.h',
    'avx2_adjacent_control.cpp',
    'avx2_adjacent_runtime.py',
    'avx2_adjacent_schedule.h',
    'PairedGroupTiming.h',
    'paired_timer_clock.cpp'}
FROZEN_FILES = SOURCE_FILES | set(PINS) | {PLAN, 'check.sh', 'check-condition.sh'}


def read(path):
    require(path.stat().st_size < 4*1024**2, 'oversized JSON')
    return parse(path.read_text())


def validate_plan(plan):
    equal(plan, PROTOCOL)


def validate_pins(pins):
    equal(sorted(pins), ['files','schema'])
    equal(pins['schema'], 'leopard-tower-pins/v1')
    equal(sorted(pins['files']), sorted(FROZEN_FILES))
    for name, value in pins['files'].items():
        require(type(value) is str and len(value)==64 and
                all(c in '0123456789abcdef' for c in value), 'invalid SHA-256')
        if name in PINS: equal(value, PINS[name])


def expected_record(cell, order, measured):
    require(type(measured) is bool, 'measured flag')
    profile = 'native' if order=='NNNN' else 'original' if order=='PPPP' else 'release'
    row = expected(profile, cell, order, PROTOCOL['groups'][cell],
                   '--exercise' if measured else '--check', 'steady')
    del row['samples']
    row['timed'] = measured
    return row


def validate(row, cell, order, measured):
    equal({k:v for k,v in row.items() if k!='samples'}, expected_record(cell,order,measured))
    samples = row['samples']
    require(type(samples) is list and len(samples)==(84 if measured else 0), 'sample count')
    for sample in samples:
        require(type(sample) is list and len(sample)==2, 'sample pair')
        elapsed, average = sample
        require(type(elapsed) is int and 0 < elapsed <= 2**53-1, 'elapsed bounds')
        require(type(average) is float and math.isfinite(average) and
                average == elapsed / PROTOCOL['groups'][cell], 'exact group normalization')


def comparisons(cell):
    return dict(CROSS, **NATIVE_CROSS) if cell in TARGETS else CROSS


def schedule():
    for cell in range(9):
        for round_id in range(3):
            for order in PAIRS:
                yield dict(cell=cell, round=round_id, comparison='paired_'+order, slot=0, order=order)
            for comparison, orders in comparisons(cell).items():
                for slot, order in enumerate(orders):
                    yield dict(cell=cell, round=round_id, comparison=comparison, slot=slot, order=order)


def label(item):
    return f"cell-{item['cell']}-round-{item['round']}-{item['comparison']}-slot-{item['slot']}"


def within(row, reverse=False):
    values = [v[1] for v in row['samples']]
    ratios = []
    for i in range(0,84,4):
        a,b,c,d = values[i:i+4]
        ratio = math.sqrt((a/b)*(d/c))
        ratios.append(1/ratio if reverse else ratio)
    return statistics.median(ratios)


def analyze(rows):
    require(type(rows) is list and len(rows)==714, 'partial attempt has no analysis')
    ratios = [dict() for _ in range(9)]
    costs = []
    def add(cell, key, value):
        ratios[cell].setdefault(key, []).append(value)
    for row, item in zip(rows, schedule()):
        equal(sorted(row), sorted([*item, 'sibling_delta', 'record']))
        equal({k:row[k] for k in item}, item)
        equal(row['sibling_delta'], 0)
        cell, comp, slot = item['cell'], item['comparison'], item['slot']
        validate(row['record'], cell, item['order'], True)
        if comp.startswith('paired_'):
            add(cell, comp, within(row['record'], item['order']=='1001'))
        else:
            costs.append(statistics.median(v[1] for v in row['record']['samples']))
            if comp.startswith('same_'):
                add(cell, 'within_'+comp+'_'+str(slot), within(row['record']))
            if slot==3:
                add(cell, comp, math.sqrt((costs[0]/costs[1])*(costs[3]/costs[2])))
                costs = []
    cells = [dict(cell=i, role=PROTOCOL['cells'][i][6], round_ratios=values,
        ratios={k:math.exp(statistics.mean(math.log(v) for v in data)) for k,data in values.items()})
        for i,values in enumerate(ratios)]
    controls = all(1/1.02<=v<=1.02 for c in cells for k,v in c['ratios'].items()
                   if k.startswith(('same_', 'within_')))
    neighbor_keys = ('paired_0110','paired_1001','original_on')
    unchanged = all(1/1.02<=cells[i]['ratios'][key]<=1.02 for i in UNCHANGED for key in neighbor_keys)
    affected = all(cells[i]['ratios'][key]>=1/1.02 for i in AFFECTED for key in neighbor_keys)
    def gains(keys):
        return all(cells[i]['ratios'][key]>=1.10 and min(cells[i]['round_ratios'][key])>1
                   for i in TARGETS for key in keys)
    targets = gains(('paired_0110','paired_1001'))
    original = gains(('original_on',))
    native = gains(('native_on',))
    decision = ('inconclusive_controls' if not controls else 'reject_unchanged_neighbor' if not unchanged else
                'reject_affected_neighbor' if not affected else 'reject_target_gate' if not targets else
                'reject_original_gate' if not original else 'reject_native_gate' if not native else
                'qualify_production_candidate')
    return dict(cells=cells, controls_pass=controls, unchanged_neighbors_pass=unchanged,
                affected_neighbors_pass=affected, targets_pass=targets, original_pass=original,
                native_pass=native, decision=decision, production_promotion=False,
                confidence_intervals=False, authoritative_v19=False)


def verify(bundle, pins):
    validate_pins(pins)
    for name, value in pins['files'].items():
        path = bundle/name
        require(path.is_file() and not path.is_symlink() and not path.stat().st_mode & 0o222,
                'mutable input: '+name)
        equal(digest(path), value)


def check_preregistration(bundle, commit):
    require(len(commit)==40 and all(c in '0123456789abcdef' for c in commit), 'commit')
    for name in SOURCE_FILES | {PLAN}:
        require(subprocess.check_output(['git','show',commit+':experiments/leopard2/gf16_high_encode/'+name])
                == (bundle/name).read_bytes(), 'preregistration differs: '+name)
    subprocess.run(['git','merge-base','--is-ancestor',commit,
                    'origin/codex/claude-fable-5-1-audit'], check=True)


def run(bundle, output, commit):
    require(output==bundle.parent/'attempt1'==Path(PROTOCOL['attempt_root']), 'fixed attempt path')
    output.mkdir(mode=0o700)  # Exclusive even when preflight fails. No retries.
    state = dict(schema='leopard-tower-attempt/v1', preregistration=commit,
                 preflight=[], invocations=[], complete=False)
    locks = []
    try:
        plan, pins = read(bundle/PLAN), read(bundle/'pins.json')
        validate_plan(plan); verify(bundle, pins); check_preregistration(bundle, commit)
        equal(host_identity(), HOST)
        for cpu in (26,90):
            equal(Path(f'/sys/devices/system/cpu/cpu{cpu}/topology/thread_siblings_list')
                  .read_text().strip(), '26,90')
        os.sched_setaffinity(0,{0})
        state.update(pins=pins, plan_sha256=digest(bundle/PLAN), host=host_identity())
        lease = Path(f'/run/user/{os.getuid()}/leopard2-cpu-leases')
        require(lease.is_dir() and not lease.is_symlink() and lease.stat().st_uid==os.getuid()
                and lease.stat().st_mode & 0o777==0o700, 'lease directory')
        for path in (Path('/tmp/leopard-gf8-authoritative.lock'),
                     lease/f'leopard2-cpu-pair-{os.getuid()}-26-90.lock'):
            fd = os.open(path, os.O_RDONLY|os.O_CREAT|os.O_NOFOLLOW, 0o600)
            locks.append(fd); fcntl.flock(fd, fcntl.LOCK_EX|fcntl.LOCK_NB)
        verify(bundle,pins)
        env = dict(PATH='/usr/bin:/bin', LANG='C', LC_ALL='C', OMP_NUM_THREADS='1',
                   OMP_DYNAMIC='FALSE', OMP_THREAD_LIMIT='1')

        def condition(name):
            with (output/(name+'.stdout')).open('xb') as out, (output/(name+'.stderr')).open('xb') as err:
                result = subprocess.run(['/bin/bash',str(bundle/'check-condition.sh')],
                                        stdout=out,stderr=err,timeout=30)
            require(result.returncode==0 and (output/(name+'.stderr')).stat().st_size==0, 'condition')

        def invoke(cell, order, name, measured):
            command = ['/usr/bin/taskset','-c','26','/usr/bin/prlimit','--cpu=30:30','--core=0:0',
                       '--fsize=1048576:1048576','--',str(bundle/('native' if order=='NNNN' else 'original' if order=='PPPP' else 'current')),
                       '--measure' if measured else '--check',str(cell),order,str(plan['groups'][cell])]
            verify(bundle,pins)
            before = sibling_ticks(90)
            outpath, errpath = output/(name+'.stdout'), output/(name+'.stderr')
            with outpath.open('xb') as out, errpath.open('xb') as err:
                result = subprocess.run(command,stdout=out,stderr=err,env=env,timeout=60)
            delta = sibling_ticks(90)-before
            require(result.returncode==0 and errpath.stat().st_size==0, 'child: '+name)
            row = read(outpath)
            verify(bundle,pins)
            return row,delta

        condition('condition-before')
        for cell in range(9):
            for order in ('NNNN','PPPP','0000','1111'):
                row,_ = invoke(cell,order,f'check-{cell}-{order}',False)
                validate(row,cell,order,False); state['preflight'].append(row)
        check_passive(state,plan)
        for item in schedule():
            row,delta = invoke(item['cell'],item['order'],label(item),True)
            state['invocations'].append(dict(item,record=row,sibling_delta=delta))
            validate(row,item['cell'],item['order'],True); equal(delta,0)
            if item['slot']==3 and item['comparison']==('same_native' if item['cell'] in TARGETS else 'same_original'):
                print(f"cell {item['cell']} round {item['round']}: paired and cross-process controls retained",flush=True)
        condition('condition-after')
        equal((output/'condition-before.stdout').read_text(),(output/'condition-after.stdout').read_text())
        verify(bundle,pins)
        state.update(analysis=analyze(state['invocations']),complete=True)
    except BaseException as error:
        state['failure'] = f'{type(error).__name__}: {error}'
        raise
    finally:
        with (output/'attempt.json').open('x') as stream:
            # Keep the complete 714-process record within the unchanged 4 MiB
            # parser bound. Only whitespace changes; no samples are omitted.
            json.dump(state,stream,separators=(',',':'),allow_nan=False); stream.write('\n')
        for fd in reversed(locks): os.close(fd)


if __name__=='__main__':
    require(len(sys.argv)==4,'usage: run_tower_screen.py frozen attempt1 pushed_commit')
    run(Path(sys.argv[1]).resolve(),Path(sys.argv[2]).resolve(),sys.argv[3])

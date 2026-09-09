#!/usr/bin/env python3
"""Independent raw replay. Never imports collector, runs a codec, or reads a clock."""
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys

# Deliberately no collector or qualification-validator imports.
SPEC = json.loads(r'''
{
  "schema": "leopard-avx2-adjacent-plan/v1",
  "bead": "leopard-79h.38.5.4.18.3.2",
  "host": {
    "hostname": "work",
    "kernel": "6.8.0-137-generic",
    "vendor_id": "AuthenticAMD",
    "cpu family": "26",
    "model": "8",
    "model name": "AMD Ryzen Threadripper 9980X 64-Cores"
  },
  "cpu": 26,
  "sibling": 90,
  "controller_cpu": 0,
  "passive_seconds": 10,
  "condition": "slipgate-disabled-20260909",
  "attempt_budget": 1,
  "attempt_root": "/tmp/leopard-adjacent-screen.04vrZl/attempt1",
  "rounds": 3,
  "sample_passes": 21,
  "spans_per_process": 84,
  "groups": [
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    256,
    1
  ],
  "paired_orders": [
    "0110",
    "1001"
  ],
  "cross_orders": {
    "same_off": [
      "0000",
      "0000",
      "0000",
      "0000"
    ],
    "same_on": [
      "1111",
      "1111",
      "1111",
      "1111"
    ],
    "original_on": [
      "PPPP",
      "1111",
      "1111",
      "PPPP"
    ],
    "original_off": [
      "PPPP",
      "0000",
      "0000",
      "PPPP"
    ],
    "same_original": [
      "PPPP",
      "PPPP",
      "PPPP",
      "PPPP"
    ]
  },
  "native_cross_orders": {
    "native_on": [
      "NNNN",
      "1111",
      "1111",
      "NNNN"
    ],
    "same_native": [
      "NNNN",
      "NNNN",
      "NNNN",
      "NNNN"
    ]
  },
  "target_cells": [
    0,
    1,
    2,
    8
  ],
  "affected_neighbors": [
    3,
    4
  ],
  "unchanged_neighbors": [
    5,
    6,
    7
  ],
  "native_cells": [
    0,
    1,
    2,
    8
  ],
  "cells": [
    [
      1000,
      200,
      32768,
      2,
      3,
      "encode",
      "target"
    ],
    [
      1000,
      199,
      65536,
      2,
      3,
      "encode",
      "target"
    ],
    [
      1000,
      200,
      65536,
      2,
      3,
      "encode",
      "target"
    ],
    [
      4096,
      512,
      4096,
      2,
      3,
      "encode",
      "affected_neighbor"
    ],
    [
      1000,
      199,
      32768,
      2,
      0,
      "encode",
      "affected_neighbor"
    ],
    [
      1000,
      200,
      32768,
      2,
      6,
      "encode",
      "unchanged_neighbor"
    ],
    [
      1000,
      200,
      32768,
      2,
      0,
      "encode",
      "unchanged_neighbor"
    ],
    [
      17,
      7,
      64,
      1,
      3,
      "encode",
      "unchanged_neighbor"
    ],
    [
      1000,
      200,
      32768,
      2,
      3,
      "one_item_batch",
      "target"
    ]
  ],
  "timed_processes": 690,
  "preflights": 36,
  "cross_control_aggregates": 31,
  "within_control_aggregates": 124,
  "paired_estimator": "median21_sqrt_off_product_over_on_product",
  "within_estimator": "median21_sqrt_outer_product_over_inner_product_each_samepath_process_slot",
  "process_estimator": "median84_all_group_averages",
  "cross_estimator": "sqrt_cost0_cost3_over_cost1_cost2",
  "round_aggregation": "geometric_mean_three_rounds_each_comparison_separately",
  "minimum_gain": 1.05,
  "original_minimum_gain": 1.05,
  "native_minimum_gain": 1.05,
  "equivalence_bound": 1.02,
  "affected_neighbor_minimum": 0.9803921568627451,
  "original_off_diagnostic_only": true,
  "every_target_round_positive": true,
  "production_promotion": false,
  "confidence_intervals": false,
  "authoritative_v19": false,
  "pristine_off_claim": false,
  "codec_commit": "45e2effd869859c9b3aa48190eff6f4738817c61",
  "native_commit": "6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198",
  "runtime_commit": "bd4c17594fa289b4b4f517f0d113b0f74bd68ebb",
  "timer_commit": "a6634c781be0d4b31248de8d8d2a30a1e00511ef",
  "artifact_sha256": {
    "native": "fafc027e811613e5a5e82bdc0ad6d94e0081491dc6adc7b669c92a83a95b755e",
    "original": "e45a2d4b305f5d16a75c548e978b11af4a50e6d28482d05ae9c95b5e0ba95a71",
    "current": "b9089340fab39f6261dc2963b8fd0d7b14acf5dbeb8f0e556a1679302e74c3af",
    "native.a": "3f83c55599cbfbeea5d15adf3627405bcde1cbed2ddf4a672fbc75b186f5c4c1",
    "original.a": "89f33d3d43f792de078c38469f9406dded2b69126af2c6fc8db6e2acb77cd334",
    "current.a": "542c340e53c1c607feb44e4d8b048bdcef41a7a04c2829c9fb62d27acd624644",
    "build.json": "f02f61e9a9ced23552e2d720478d7699eba82ce38011efb73acdb3c475809592",
    "checks.json": "12f817fe7e0b1c9879e53f8b5ff683b9333b60303c9f820d8bee6e091aa7300f",
    "qualification.json": "4e328eb98f873e2b7abf8255676dc44b34dd27af3e7d16436933032ac816c307",
    "qualification.manifest": "f3e58727918ded2e06b62c250ff1b19bde5d7860c107dbf9a8166a67fa0f4dfa",
    "l2-driver.o": "94128723910cea77359a781767219dcc847fd15b14f40ed91972e172a4fec5da",
    "native-driver.o": "d6c7ca19e9454a71328f51b3abca781a7cd26240c6cc591ccca57eaee28f094f",
    "leopard2.cpp": "93412028272454b2457f9928a6bc4c34684591f5bee2c882a3ad0000fd50cc70",
    "Leopard2Direct.h": "da4b39084cf7be12e824f46e799386d742120fe50a852c8b6a12dd9fe0c6b423"
  }
}
''')
CODES = {'native': '3f83c55599cbfbeea5d15adf3627405bcde1cbed2ddf4a672fbc75b186f5c4c1', 'original': '89f33d3d43f792de078c38469f9406dded2b69126af2c6fc8db6e2acb77cd334',
         'release': '542c340e53c1c607feb44e4d8b048bdcef41a7a04c2829c9fb62d27acd624644'}
INPUTS = ['8b78decd04e67d27','ba910426bb91dec4','ba910426bb91dec4','8ab693ec532fd3fa',
          '8b78decd04e67d27','8b78decd04e67d27','8b78decd04e67d27','bef540029ca6f65b','8b78decd04e67d27']
OUTPUTS = ['50d5c45deb48367a','0b985bca4df709d5','b1d946fa500bbf12','12c692cff473f0ee',
           '00648d9dbf0f2b20','50d5c45deb48367a','50d5c45deb48367a','83958df25bd2c455','50d5c45deb48367a']

PLAN = 'avx2_adjacent_screen_plan.json'
SOURCES = {
    'run_avx2_adjacent_screen.py',
    'replay_avx2_adjacent_screen.py',
    'test_avx2_adjacent_screen.py',
    'freeze_avx2_adjacent_screen.py',
    'avx2_adjacent_screen_method.md',
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
ARTIFACTS = {"native","original","current","native.a","original.a","current.a","build.json","checks.json","qualification.json","qualification.manifest","l2-driver.o","native-driver.o","leopard2.cpp","Leopard2Direct.h"}
FILES = SOURCES | ARTIFACTS | {PLAN,'check.sh','check-condition.sh'}
GROUPS = [1]*7+[256,1]


def require(ok, message):
    if not ok: raise ValueError(message)


def same(a,b):
    require(json.dumps(a,sort_keys=True,allow_nan=False)==json.dumps(b,sort_keys=True,allow_nan=False),
            'typed identity mismatch')


def read(path):
    require(path.stat().st_size<4*1024**2,'oversized JSON')
    def unique(pairs):
        result = {}
        for k,v in pairs:
            require(k not in result,'duplicate JSON key'); result[k] = v
        return result
    return json.loads(path.read_text(),object_pairs_hook=unique)


def sha(path):
    with path.open('rb') as stream: return hashlib.file_digest(stream,'sha256').hexdigest()


def record(row, cell, order, measured):
    require(type(cell) is int and cell in range(9) and type(measured) is bool, 'cell/mode')
    require(order in ('NNNN','PPPP','0110','1001','0000','1111'), 'schedule')
    p = 'native' if order=='NNNN' else 'original' if order=='PPPP' else 'release'
    k,r,b,field,backend,_,_ = SPEC['cells'][cell]
    count = 1+25*GROUPS[cell] if measured else 1
    scratch = ([16777216,33554432,33554432,4194304,16777216,16777216,16777216,1024,16777216]
               if p=='native' else [16808512]*3+[4308992]+[16808512]*3+[1728,16808512])[cell]
    zero = dict(calls=[[0,0],[0,0]],blocks=[[0,0],[0,0]])
    fixed = dict(schema='leopard-adjacent-public/v1',
        codec=p+':'+CODES[p],cell=cell,k=k,r=r,bytes=b,field=field,backend=backend,
        api='leo_encode' if p=='native' else 'leo2_encode_batch_one_item' if cell==8 else 'leo2_encode',
        schedule=order,group=GROUPS[cell],warmup_passes=4 if measured else 0,
        exercise_passes=21 if measured else 0,encode_calls=count*4,selections=104 if measured else 4,
        per_slot_calls=[count]*4,probes=[int(p!='native' and cell==6)]*4,
        scratch_bytes=scratch,input_hash=INPUTS[cell],output_hash=OUTPUTS[cell],
        clock_source='steady',timed=measured,default_enabled=False,traced=False,
        pair_probes=[zero]*3,pair_totals=zero)
    same({k:v for k,v in row.items() if k!='samples'},fixed)
    samples = row['samples']
    require(type(samples) is list and len(samples)==(84 if measured else 0),'sample inventory')
    for sample in samples:
        require(type(sample) is list and len(sample)==2,'sample shape')
        ns,average = sample
        require(type(ns) is int and ns in range(1,2**53),'integer interval')
        require(type(average) is float and math.isfinite(average) and average*GROUPS[cell]==ns,
                'normalized sample')


def median(values):
    values = sorted(values)
    n = len(values)
    return values[n//2] if n%2 else (values[n//2-1]+values[n//2])/2


def internal(row, reverse=False):
    spans = [x[1] for x in row['samples']]
    values = []
    for i in range(21):
        a,b,c,d = spans[4*i:4*i+4]
        values.append(math.sqrt((b*c)/(a*d)) if reverse else math.sqrt((a*d)/(b*c)))
    return median(values)


def derive(rows):
    require(type(rows) is list and len(rows)==690,'complete ordered attempt required')
    cursor,cells = 0,[]
    for cell in range(9):
        rounds = {}
        for rnd in range(3):
            comparisons = [('paired_0110',['0110']),('paired_1001',['1001']),
                           ('same_off',['0000']*4),('same_on',['1111']*4),
                           ('original_on',['PPPP','1111','1111','PPPP']),
                           ('original_off',['PPPP','0000','0000','PPPP']),
                           ('same_original',['PPPP']*4)]
            if cell in (0,1,2,8): comparisons += [('native_on',['NNNN','1111','1111','NNNN']),('same_native',['NNNN']*4)]
            for name,orders in comparisons:
                costs = []
                for slot,order in enumerate(orders):
                    row = rows[cursor]; cursor += 1
                    same(sorted(row),['cell','comparison','order','record','round','sibling_delta','slot'])
                    same({k:v for k,v in row.items() if k!='record'},
                         dict(cell=cell,round=rnd,comparison=name,slot=slot,order=order,sibling_delta=0))
                    record(row['record'],cell,order,True)
                    if name.startswith('paired_'):
                        rounds.setdefault(name,[]).append(internal(row['record'],order=='1001'))
                    else:
                        costs.append(median([s[1] for s in row['record']['samples']]))
                        if name.startswith('same_'):
                            rounds.setdefault('within_'+name+'_'+str(slot),[]).append(internal(row['record']))
                if len(orders)==4:
                    a,b,c,d = costs
                    rounds.setdefault(name,[]).append(math.sqrt((a*d)/(b*c)))
        cells.append(dict(cell=cell,role='target' if cell in (0,1,2,8) else 'affected_neighbor' if cell in (3,4) else 'unchanged_neighbor',round_ratios=rounds,
                          ratios={k:math.prod(v)**(1/3) for k,v in rounds.items()}))
    controls = all(1/1.02<=v<=1.02 for c in cells for k,v in c['ratios'].items()
                   if k.startswith(('same_','within_')))
    keys = ('paired_0110','paired_1001','original_on')
    unchanged = all(1/1.02<=cells[i]['ratios'][key]<=1.02 for i in (5,6,7) for key in keys)
    affected = all(cells[i]['ratios'][key]>=1/1.02 for i in (3,4) for key in keys)
    passes = []
    for group in (('paired_0110','paired_1001'),('original_on',),('native_on',)):
        passes.append(all(cells[i]['ratios'][key]>=1.05 and all(x>1 for x in cells[i]['round_ratios'][key])
                          for i in (0,1,2,8) for key in group))
    targets,original,native = passes
    decision = ('inconclusive_controls' if not controls else 'reject_unchanged_neighbor' if not unchanged else
                'reject_affected_neighbor' if not affected else 'reject_target_gate' if not targets else
                'reject_original_gate' if not original else 'reject_native_gate' if not native else
                'qualify_production_candidate')
    return dict(cells=cells,controls_pass=controls,unchanged_neighbors_pass=unchanged,
                affected_neighbors_pass=affected,targets_pass=targets,original_pass=original,
                native_pass=native,decision=decision,production_promotion=False,
                confidence_intervals=False,authoritative_v19=False)


def compare_analysis(actual, derived):
    same({k:v for k,v in actual.items() if k!='cells'},{k:v for k,v in derived.items() if k!='cells'})
    same(len(actual['cells']),9)
    for a,b in zip(actual['cells'],derived['cells']):
        same(sorted(a),sorted(b)); same(a['cell'],b['cell']); same(a['role'],b['role'])
        for category in ('round_ratios','ratios'):
            same(sorted(a[category]),sorted(b[category]))
            for key in a[category]:
                left = a[category][key] if category=='round_ratios' else [a[category][key]]
                right = b[category][key] if category=='round_ratios' else [b[category][key]]
                same(len(left),len(right))
                require(all(type(v) is float and math.isfinite(v) and math.isclose(v,w,rel_tol=2e-14,abs_tol=0)
                            for v,w in zip(left,right)),'analysis differs from raw')


def resource(path):
    lines = path.read_text().splitlines()
    require(lines.count('memory.peak')==1 and '\tExit status: 0' in lines,'resource exit')
    i = lines.index('memory.peak'); peak = int(lines[i+1])
    require(0<peak<=268435456,'memory bound')
    same(lines[i+2:],['memory.max','268435456','memory.events','low 0','high 0','max 0','oom 0',
                     'oom_kill 0','oom_group_kill 0','memory.swap.current','0','memory.swap.max','0'])
    return peak


def inputs(root, commit=None):
    frozen = root/'frozen'
    plan,pins = read(frozen/PLAN),read(frozen/'pins.json')
    same(sorted(pins),['files','schema']); same(pins['schema'],'leopard-avx2-adjacent-pins/v1')
    same(sorted(pins['files']),sorted(FILES))
    same(sorted(p.name for p in frozen.iterdir()),sorted(FILES|{'pins.json'}))
    for name,value in pins['files'].items():
        path = frozen/name
        require(path.is_file() and not path.is_symlink() and not path.stat().st_mode & 0o222,'mutable input')
        same(sha(path),value)
    same(plan,SPEC)
    same(sorted(plan['artifact_sha256']),sorted(ARTIFACTS))
    for name,value in plan['artifact_sha256'].items(): same(pins['files'][name],value)
    # This exact qualified manifest independently binds all executable/record copies.
    same(sha(frozen/'qualification.manifest'),'f3e58727918ded2e06b62c250ff1b19bde5d7860c107dbf9a8166a67fa0f4dfa')
    manifest = dict(line.split('  ',1)[::-1] for line in (frozen/'qualification.manifest').read_text().splitlines())
    for name,original in [('native','build/native/plain'),('original','build/original/plain'),('current','build/release/plain'),
            ('native.a','build/native/codec.a'),('original.a','build/original/codec.a'),('current.a','build/release/codec.a'),
            ('l2-driver.o','build/l2-release-objects/driver.o'),('native-driver.o','build/native-objects/driver.o'),
            ('build.json','build/build.json'),('checks.json','checks/checks.json'),('qualification.json','replay.json')]:
        same(pins['files'][name],manifest[original])
    build = read(frozen/'build.json')
    for name in ('avx2_adjacent_public.cpp','avx2_adjacent_public_link.cpp','avx2_adjacent_public_link.h',
                 'avx2_adjacent_control.h','PairedGroupTiming.h','paired_timer_clock.cpp','Leopard2Direct.h'):
        matches = [v for k,v in build['inputs'].items() if Path(k).name==name]
        same(matches,[pins['files'][name]])
    same(pins['files']['leopard2.cpp'],'93412028272454b2457f9928a6bc4c34684591f5bee2c882a3ad0000fd50cc70')
    # The manifest binds the driver objects; link recipes prove all three
    # Release comparisons actually consume the one shared compiled driver.
    for profile in ('native','original','release'):
        shared = 'native-objects' if profile=='native' else 'l2-release-objects'
        matches = [c for c in build['commands'] if c[-1].endswith('/'+profile+'/plain')]
        require(len(matches)==1,'unique qualified plain link')
        wanted = [x for x in matches[0] if x.endswith('/'+shared+'/driver.o')]
        require(len(wanted)==1,'shared qualified driver')
        same(build['artifacts'][shared+'/driver.o'],pins['files']['native-driver.o' if profile=='native' else 'l2-driver.o'])
    if commit is not None:
        require(len(commit)==40 and all(c in '0123456789abcdef' for c in commit),'commit')
        for name in SOURCES|{PLAN}:
            same(subprocess.check_output(['git','show',commit+':experiments/leopard2/gf16_high_encode/'+name]).decode(),
                 (frozen/name).read_text())
    return plan,pins


def replay(root, commit):
    plan,pins = inputs(root,commit)
    attempt = root/'attempt1'
    state = read(attempt/'attempt.json')
    same(sorted(state),['analysis','complete','host','invocations','passive','pins','plan_sha256','preflight',
                        'preregistration','schema'])
    same(state['schema'],'leopard-avx2-adjacent-attempt/v1'); same(state['preregistration'],commit)
    same(state['pins'],pins); same(state['plan_sha256'],pins['files'][PLAN]); same(state['host'],plan['host'])
    same(state['complete'],True)
    passive = state['passive']
    same(sorted(passive),['after','before','elapsed_ns'])
    require(type(passive['before']) is int and passive['before']>=0 and type(passive['after']) is int and
            passive['before']==passive['after'] and type(passive['elapsed_ns']) is int and
            passive['elapsed_ns']>=10000000000,'passive isolation')
    same(len(state['preflight']),36)
    names = []
    for cell in range(9):
        for i,order in enumerate(('NNNN','PPPP','0000','1111')):
            name = f'check-{cell}-{order}'; names.append(name)
            row = read(attempt/(name+'.stdout'))
            same(row,state['preflight'][4*cell+i]); record(row,cell,order,False)
    derived = derive(state['invocations']); compare_analysis(state['analysis'],derived)
    for row in state['invocations']:
        name = f"cell-{row['cell']}-round-{row['round']}-{row['comparison']}-slot-{row['slot']}"
        names.append(name); same(read(attempt/(name+'.stdout')),row['record'])
    condition = []
    for unit in ('slipgate-catstream.service','slipgate-obs.service','slipgate-headless-x.service'):
        condition += ['MainPID=0','Id='+unit,'ActiveState=inactive','UnitFileState=disabled']
    for container in ('3fad3dfb247b3d4f3991e0a6801fba0805aae0901bdc89018841d4707a2d888c',
                      '002ba0a05bf860b12f2d799b49da1976aed663cb1494b2f8236751e18dc50182'):
        condition.append(container+' exited|false|no')
    for name in ('condition-before','condition-after'):
        names.append(name); same((attempt/(name+'.stdout')).read_text().splitlines(),condition)
    for name in names: same((attempt/(name+'.stderr')).read_text(),'')
    same(sorted(p.name for p in attempt.iterdir()),sorted(['attempt.json']+
         [name+'.'+suffix for name in names for suffix in ('stdout','stderr')]))
    return dict(preregistration=commit,frozen_inputs=len(FILES),preflights=36,timed_invocations=690,
        spans=690*84,cross_control_aggregates=31,within_control_aggregates=124,sibling_nonidle_jiffies=0,
        analysis=derived,memory_peak=resource(root/'attempt1.log'))


if __name__=='__main__':
    require(len(sys.argv)==3,'usage: replay_avx2_adjacent_screen.py ROOT PREREGISTRATION')
    print(json.dumps(replay(Path(sys.argv[1]).resolve(),sys.argv[2]),indent=2,sort_keys=True))

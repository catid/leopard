#!/usr/bin/env python3
"""Independent raw replay. Never imports collector, runs a codec, or reads a clock."""
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys

from verify_auto_r19932_checks import public_expected

PLAN = 'paired_r19932_screen_plan.json'
SOURCES = {'run_paired_r19932_screen.py','replay_paired_r19932_screen.py','test_paired_r19932_screen.py',
    'freeze_paired_r19932_screen.py','paired_r19932_screen_method.md','run_split_cache_screen.py',
    'verify_paired_r19932.py','verify_auto_r19932_checks.py','verify_auto_gfni_boundary_checks.py',
    'paired_timer_r19932.cpp','PairedGroupTiming.h','paired_timer_clock.cpp'}
ARTIFACTS = {'native','current','native.a','current.a','build.json','checks.json','qualification.json',
             'qualification.manifest','leopard2.cpp','Leopard2Direct.h'}
FILES = SOURCES | ARTIFACTS | {PLAN,'check.sh','check-condition.sh'}
GROUPS = [1]*8+[256]


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
    p = 'native' if order=='NNNN' else 'release'
    old = public_expected(p,0,cell)
    count = 1+25*GROUPS[cell] if measured else 1
    fixed = dict(schema='leopard-paired-timer-r19932/v1',
        codec=p+':'+old['codec_commit'].split(':')[-1],cell=cell,k=old['k'],r=old['r'],bytes=old['bytes'],
        api=old['api'],schedule=order,group=GROUPS[cell],warmup_passes=4 if measured else 0,
        exercise_passes=21 if measured else 0,encode_calls=count*4,selections=104 if measured else 4,
        per_slot_calls=[count]*4,probes=[0 if p=='native' else int(2<=cell<=4 or (cell<2 and s=='1')) for s in order],
        scratch_bytes=old['scratch_bytes'],input_hash=old['input_hash'],output_hash=old['output_hash'],
        clock_source='steady',timed=measured,default_enabled=False)
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
    require(type(rows) is list and len(rows)==318,'complete ordered attempt required')
    cursor,cells = 0,[]
    for cell in range(9):
        rounds = {}
        for rnd in range(3):
            comparisons = [('paired_0110',['0110']),('paired_1001',['1001']),
                           ('same_off',['0000']*4),('same_on',['1111']*4)]
            if cell<2: comparisons += [('native_on',['NNNN','1111','1111','NNNN']),('same_native',['NNNN']*4)]
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
        cells.append(dict(cell=cell,role='target' if cell<2 else 'unchanged_neighbor',round_ratios=rounds,
                          ratios={k:math.prod(v)**(1/3) for k,v in rounds.items()}))
    controls = all(1/1.02<=v<=1.02 for c in cells for k,v in c['ratios'].items()
                   if k.startswith(('same_','within_')))
    neighbors = all(1/1.02<=c['ratios'][k]<=1.02 for c in cells[2:] for k in ('paired_0110','paired_1001'))
    targets = all(c['ratios'][k]>=1.05 and all(v>1 for v in c['round_ratios'][k])
                  for c in cells[:2] for k in ('paired_0110','paired_1001'))
    native = all(c['ratios']['native_on']>=1.05 and all(v>1 for v in c['round_ratios']['native_on'])
                 for c in cells[:2])
    decision = ('inconclusive_controls' if not controls else 'reject_neighbor_gate' if not neighbors else
                'reject_target_gate' if not targets else 'reject_native_gate' if not native else
                'qualify_default_on_artifact')
    return dict(cells=cells,controls_pass=controls,neighbors_pass=neighbors,targets_pass=targets,
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
    same(sorted(pins),['files','schema']); same(pins['schema'],'leopard-paired-r19932-pins/v1')
    same(sorted(pins['files']),sorted(FILES))
    same(sorted(p.name for p in frozen.iterdir()),sorted(FILES|{'pins.json'}))
    for name,value in pins['files'].items():
        path = frozen/name
        require(path.is_file() and not path.is_symlink() and not path.stat().st_mode & 0o222,'mutable input')
        same(sha(path),value)
    same(sorted(plan['artifact_sha256']),sorted(ARTIFACTS))
    for name,value in plan['artifact_sha256'].items(): same(pins['files'][name],value)
    # This exact qualified manifest independently binds all executable/record copies.
    same(sha(frozen/'qualification.manifest'),'6c86e08b8fcef3a7d3a9c345b5820bf9e8af1310133b6c8835cd4adcc2fab37b')
    manifest = dict(line.split('  ',1)[::-1] for line in (frozen/'qualification.manifest').read_text().splitlines())
    for name,original in [('native','build/native/plain'),('current','build/release/plain'),
            ('native.a','build/native/codec.a'),('current.a','build/release/codec.a'),
            ('build.json','build/build.json'),('checks.json','checks/checks.json'),('qualification.json','replay.json')]:
        same(pins['files'][name],manifest[original])
    build = read(frozen/'build.json')
    for name in ('paired_timer_r19932.cpp','PairedGroupTiming.h','paired_timer_clock.cpp','Leopard2Direct.h'):
        matches = [v for k,v in build['inputs'].items() if Path(k).name==name]
        same(matches,[pins['files'][name]])
    same(pins['files']['leopard2.cpp'],'93412028272454b2457f9928a6bc4c34684591f5bee2c882a3ad0000fd50cc70')
    same([plan[k] for k in ('schema','bead','cpu','sibling','controller_cpu','rounds','sample_passes',
        'spans_per_process','groups','attempt_budget','attempt_root','passive_seconds','minimum_gain',
        'native_minimum_gain','equivalence_bound','timed_processes','preflights','cross_control_aggregates',
        'within_control_aggregates','every_target_round_positive','production_promotion','confidence_intervals',
        'authoritative_v19','pristine_off_claim')],
        ['leopard-paired-r19932-plan/v1','leopard-79h.38.5.4.19.1.2',26,90,0,3,21,84,GROUPS,1,
         '/tmp/leopard-paired-integration.la42Cz/attempt1',10,1.05,1.05,1.02,318,27,20,80,True,False,False,False,False])
    same(plan['paired_orders'],['0110','1001'])
    same(plan['cross_orders'],dict(same_off=['0000']*4,same_on=['1111']*4))
    same(plan['native_cross_orders'],dict(native_on=['NNNN','1111','1111','NNNN'],same_native=['NNNN']*4))
    same(plan['native_cells'],[0,1])
    wanted_cells = []
    for cell in range(9):
        old = public_expected('release',0,cell)
        api = 'one_item_batch' if cell==1 else 'explicit_avx2_encode' if cell==6 else 'gf8_encode' if cell==8 else 'encode'
        wanted_cells.append([old['k'],old['r'],old['bytes'],api,'target' if cell<2 else 'unchanged_neighbor'])
    same(plan['cells'],wanted_cells)
    same([plan[k] for k in ('paired_estimator','within_estimator','process_estimator','cross_estimator','round_aggregation')],
         ['median21_sqrt_off_product_over_on_product',
          'median21_sqrt_outer_product_over_inner_product_each_samepath_process_slot',
          'median84_all_group_averages','sqrt_cost0_cost3_over_cost1_cost2',
          'geometric_mean_three_rounds_each_comparison_separately'])
    same(plan['host'],dict(hostname='work',kernel='6.8.0-137-generic',vendor_id='AuthenticAMD',
        **{'cpu family':'26','model':'8','model name':'AMD Ryzen Threadripper 9980X 64-Cores'}))
    same(plan['condition'],'slipgate-disabled-20260909')
    same(plan['codec_commit'],'45e2effd869859c9b3aa48190eff6f4738817c61')
    same(plan['native_commit'],'6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198')
    same(plan['timer_commit'],'553d0df949abec4ca02cd067eabfed46fe92d1d7')
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
    same(state['schema'],'leopard-paired-r19932-attempt/v1'); same(state['preregistration'],commit)
    same(state['pins'],pins); same(state['plan_sha256'],pins['files'][PLAN]); same(state['host'],plan['host'])
    same(state['complete'],True)
    passive = state['passive']
    same(sorted(passive),['after','before','elapsed_ns'])
    require(type(passive['before']) is int and passive['before']>=0 and type(passive['after']) is int and
            passive['before']==passive['after'] and type(passive['elapsed_ns']) is int and
            passive['elapsed_ns']>=10000000000,'passive isolation')
    same(len(state['preflight']),27)
    names = []
    for cell in range(9):
        for i,order in enumerate(('NNNN','0000','1111')):
            name = f'check-{cell}-{order}'; names.append(name)
            row = read(attempt/(name+'.stdout'))
            same(row,state['preflight'][3*cell+i]); record(row,cell,order,False)
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
    return dict(preregistration=commit,frozen_inputs=len(FILES),preflights=27,timed_invocations=318,
        spans=318*84,cross_control_aggregates=20,within_control_aggregates=80,sibling_nonidle_jiffies=0,
        analysis=derived,memory_peak=resource(root/'attempt1.log'))


if __name__=='__main__':
    require(len(sys.argv)==3,'usage: replay_paired_r19932_screen.py ROOT PREREGISTRATION')
    print(json.dumps(replay(Path(sys.argv[1]).resolve(),sys.argv[2]),indent=2,sort_keys=True))

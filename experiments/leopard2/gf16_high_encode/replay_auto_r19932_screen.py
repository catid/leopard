#!/usr/bin/env python3
"""Collector-free retained raw replay; no native execution or timing calls."""
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys

PROFILES = ('native','off','on')
SHAPES = [(1000,199,32768),(1000,199,32768),(1000,200,32768),(1000,199,65536),
          (1000,200,65536),(1000,198,32768),(1000,199,32768),(4096,512,4096),(17,7,64)]
ORDERS = [('off_on',('off','on','on','off')),('same_off',('off',)*4),('same_on',('on',)*4)]
NATIVE_ORDERS = [('native_on',('native','on','on','native')),('same_native',('native',)*4)]
QUALIFICATION = Path('/home/catid/leopard/.research/leopard-79h/auto-r19932-qualified.515m59j2')
SOURCE_FILES = {'auto_r19932_screen.cpp','run_auto_r19932_screen.py','test_auto_r19932_screen.py',
    'replay_auto_r19932_screen.py','test_replay_auto_r19932_screen.py','run_split_cache_screen.py',
    'verify_auto_r19932_checks.py','verify_auto_gfni_boundary_checks.py','test_verify_auto_r19932_checks.py'}
ARTIFACTS = {'native','current','native.a','current.a','expected.json'}
PLAN = 'auto_r19932_screen_plan.json'
FROZEN_FILES = SOURCE_FILES | ARTIFACTS | {PLAN,'check-condition.sh','check.sh','build.json',
    'preparation.json','qualification.json','leopard2.cpp','Leopard2Direct.h'}


def require(value,message):
    if not value: raise ValueError(message)


def same(a,b):
    require(json.dumps(a,sort_keys=True,allow_nan=False)==json.dumps(b,sort_keys=True,allow_nan=False),
            'typed identity mismatch')


def sha(path):
    with path.open('rb') as f: return hashlib.file_digest(f,'sha256').hexdigest()


def read(path):
    require(path.stat().st_size < 4*1048576,'bounded JSON')
    return json.loads(path.read_text())


def record(value,expected,measured):
    same({k:v for k,v in value.items() if k not in ('encode_calls','samples_ns')},expected)
    same(value['encode_calls'],26 if measured else 1)
    samples = value['samples_ns']
    require(type(samples) is list and len(samples)==(21 if measured else 0) and
            all(type(v) is int and v>0 for v in samples),'invalid sample list')


def derive(rows,expected):
    require(type(rows) is list and len(rows)==372,'complete ordered attempt required')
    cursor,cells = 0,[]
    for cell in range(9):
        comparisons = ORDERS+NATIVE_ORDERS if cell<2 else ORDERS
        rounds = {name:[] for name,_ in comparisons}
        for round_id in range(3):
            for comparison,order in comparisons:
                medians = []
                for slot,profile in enumerate(order):
                    row = rows[cursor]; cursor += 1
                    same({k:row[k] for k in ('cell','round','comparison','slot','profile','sibling_delta')},
                         dict(cell=cell,round=round_id,comparison=comparison,slot=slot,
                              profile=profile,sibling_delta=0))
                    record(row['record'],expected[profile][cell],True)
                    medians.append(sorted(row['record']['samples_ns'])[10])
                a,b,c,d = medians
                rounds[comparison].append(math.sqrt((a*d)/(b*c)))
        cells.append(dict(cell=cell,role='target' if cell<2 else 'unchanged_neighbor',
            round_ratios=rounds,ratios={k:math.prod(values)**(1/3) for k,values in rounds.items()}))
    controls = all(1/1.02 <= v <= 1.02 for c in cells for k,v in c['ratios'].items() if k.startswith('same_'))
    neighbors = all(1/1.02 <= c['ratios']['off_on'] <= 1.02 for c in cells[2:])
    targets = all(c['ratios']['off_on']>=1.05 and all(v>1 for v in c['round_ratios']['off_on']) for c in cells[:2])
    native = all(c['ratios']['native_on']>=1.05 and all(v>1 for v in c['round_ratios']['native_on']) for c in cells[:2])
    decision = ('inconclusive_controls' if not controls else 'reject_neighbor_gate' if not neighbors else
                'reject_target_gate' if not targets else 'reject_native_gate' if not native else 'qualify_default_on_artifact')
    return dict(cells=cells,controls_pass=controls,neighbors_pass=neighbors,targets_pass=targets,
                native_pass=native,decision=decision,production_promotion=False,confidence_intervals=False,authoritative_v19=False)


def compare_analysis(a,b):
    same({k:v for k,v in a.items() if k!='cells'},{k:v for k,v in b.items() if k!='cells'})
    require(len(a['cells'])==len(b['cells'])==9,'cell count')
    for x,y in zip(a['cells'],b['cells']):
        same(sorted(x),sorted(y)); same(x['cell'],y['cell']); same(x['role'],y['role'])
        for category in ('ratios','round_ratios'):
            same(sorted(x[category]),sorted(y[category]))
            for key in x[category]:
                left = x[category][key] if category=='round_ratios' else [x[category][key]]
                right = y[category][key] if category=='round_ratios' else [y[category][key]]
                require(len(left)==len(right) and all(type(v) is float and math.isfinite(v) and
                        math.isclose(v,w,rel_tol=2e-14,abs_tol=0) for v,w in zip(left,right)),'derived ratios')


def resource(path,limit):
    lines = path.read_text().splitlines()
    require(lines.count('memory.peak')==1 and '\tExit status: 0' in lines,'scope exit')
    i = lines.index('memory.peak'); peak = int(lines[i+1])
    require(0<peak<=limit,'peak bound')
    same(lines[i+2:],['memory.max',str(limit),'memory.events','low 0','high 0','max 0','oom 0',
        'oom_kill 0','oom_group_kill 0','memory.swap.current','0','memory.swap.max','0'])
    return peak


def inputs(root,commit=None):
    from verify_auto_r19932_checks import public_expected, replay as qualification_replay
    frozen = root/'frozen'
    pins = read(frozen/'pins.json')
    same(sorted(pins),['files','schema']); same(pins['schema'],'leopard-auto-r19932-pins/v1')
    same(sorted(pins['files']),sorted(FROZEN_FILES))
    for name,digest in pins['files'].items():
        path = frozen/name
        require(path.is_file() and not path.is_symlink() and not path.stat().st_mode & 0o222,'mutable input')
        same(sha(path),digest)
    plan = read(frozen/PLAN)
    same(plan['schema'],'leopard-auto-r19932-plan/v1'); same(plan['bead'],'leopard-79h.38.5.4.19.1')
    same(plan['codec_commit'],'45e2effd869859c9b3aa48190eff6f4738817c61')
    same(plan['native_commit'],'6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198')
    same([(c['k'],c['r'],c['bytes']) for c in plan['cells']],SHAPES)
    same([c['role'] for c in plan['cells']],['target']*2+['unchanged_neighbor']*7)
    same(plan['orders'],{k:list(v) for k,v in ORDERS})
    same(plan['native_orders'],{k:list(v) for k,v in NATIVE_ORDERS})
    same(plan['native_cells'],[0,1])
    same([plan[k] for k in ('cpu','sibling','controller_cpu','rounds','samples_per_process','attempt_budget',
        'passive_seconds','minimum_gain','native_minimum_gain','equivalence_bound','production_promotion',
        'confidence_intervals','authoritative_v19')],[26,90,0,3,21,1,10,1.05,1.05,1.02,False,False,False])
    same(plan['attempt_root'],'/tmp/leopard-auto-r19932-screen.oSPt31/attempt1')
    same(plan['condition'],'slipgate-disabled-20260909')
    same(plan['host'],dict(hostname='work',kernel='6.8.0-137-generic',vendor_id='AuthenticAMD',
        **{'cpu family':'26','model':'8','model name':'AMD Ryzen Threadripper 9980X 64-Cores'}))
    same(sorted(plan['artifact_sha256']),sorted(ARTIFACTS))
    for name,digest in plan['artifact_sha256'].items(): same(digest,pins['files'][name])
    for name,key in (('leopard2.cpp','core_sha256'),('Leopard2Direct.h','header_sha256')):
        same(pins['files'][name],plan[key]); same(sha(QUALIFICATION/'build/source'/name),plan[key])
    same(plan['qualification_manifest_sha256'],'4799b38994df2804ea77df521200fb28e66ba4fb4f1a9b1bf600d4399b0a9f18')
    same(sha(QUALIFICATION/'SHA256SUMS'),plan['qualification_manifest_sha256'])
    qualification = qualification_replay(QUALIFICATION)
    same(read(frozen/'qualification.json'),qualification)
    build,prep = read(frozen/'build.json'),read(frozen/'preparation.json')
    same(build,read(QUALIFICATION/'build/build.json'))
    same(prep,read(QUALIFICATION/'preparation/preparation.json'))
    same(sha(frozen/'auto_r19932_screen.cpp'),build['driver_pins']['auto_r19932_screen.cpp'])
    expected = read(frozen/'expected.json'); same(sorted(expected),sorted(PROFILES))
    for variant in PROFILES:
        p = 'native' if variant=='native' else 'release'; mode = int(variant=='on')
        name = 'native' if variant=='native' else 'current'
        same(sha(frozen/name),build['profiles'][p]['files']['screen'])
        same(sha(frozen/(name+'.a')),build['profiles'][p]['archive_sha256'])
        same(len(expected[variant]),9)
        for cell in range(9):
            value = dict(expected[variant][cell],encode_calls=1,samples_ns=[])
            same(value,public_expected(p,mode,cell))
            same(value,prep['public_records'][f'{p}-{mode}-{cell}-check'])
    if commit is not None:
        require(len(commit)==40 and all(c in '0123456789abcdef' for c in commit),'commit')
        for name in SOURCE_FILES | {PLAN}:
            same(subprocess.check_output(['git','show',commit+':experiments/leopard2/gf16_high_encode/'+name]).decode(),
                 (frozen/name).read_text())
    return plan,pins,expected,qualification


def replay(root,commit):
    plan,pins,expected,qualification = inputs(root,commit)
    state = read(root/'attempt1/attempt.json')
    same(state['schema'],'leopard-auto-r19932-attempt/v1'); same(state['preregistration'],commit)
    same(state['pins'],pins); same(state['host'],plan['host'])
    same(state['plan_sha256'],pins['files'][PLAN])
    require(state['complete'] is True and 'failure' not in state,'complete attempt')
    passive = state['passive']
    require(type(passive['before']) is int and passive['before']>=0 and type(passive['after']) is int and
        passive['before']==passive['after'] and type(passive['elapsed_ns']) is int and
        passive['elapsed_ns']>=10000000000,'passive isolation')
    same(len(state['preflight']),27)
    for cell in range(9):
        for i,variant in enumerate(PROFILES):
            name = f'check-{cell}-{variant}'
            actual = read(root/'attempt1'/(name+'.stdout'))
            same(actual,state['preflight'][cell*3+i]); record(actual,expected[variant][cell],False)
            same((root/'attempt1'/(name+'.stderr')).read_text(),'')
    derived = derive(state['invocations'],expected)
    compare_analysis(state['analysis'],derived)
    for row in state['invocations']:
        name = f"cell-{row['cell']}-round-{row['round']}-{row['comparison']}-slot-{row['slot']}"
        same(read(root/'attempt1'/(name+'.stdout')),row['record'])
        same((root/'attempt1'/(name+'.stderr')).read_text(),'')
    condition = []
    for unit in ('slipgate-catstream.service','slipgate-obs.service','slipgate-headless-x.service'):
        condition += ['MainPID=0','Id='+unit,'ActiveState=inactive','UnitFileState=disabled']
    for container in ('3fad3dfb247b3d4f3991e0a6801fba0805aae0901bdc89018841d4707a2d888c',
                      '002ba0a05bf860b12f2d799b49da1976aed663cb1494b2f8236751e18dc50182'):
        condition.append(container+' exited|false|no')
    for name in ('condition-before','condition-after'):
        same((root/'attempt1'/(name+'.stdout')).read_text().splitlines(),condition)
        same((root/'attempt1'/(name+'.stderr')).read_text(),'')
    for suffix in ('stdout','stderr'):
        same(len(list((root/'attempt1').glob('*.'+suffix))),401)
    return dict(preregistration=commit,frozen_inputs=len(FROZEN_FILES),preflights=27,timed_invocations=372,
        sibling_nonidle_jiffies=0,qualification=qualification,analysis=derived,
        memory_peak=resource(root/'attempt1.log',268435456))


if __name__=='__main__':
    require(len(sys.argv)==3,'usage: replay_auto_r19932_screen.py ROOT PREREGISTRATION')
    print(json.dumps(replay(Path(sys.argv[1]).resolve(),sys.argv[2]),indent=2))

#!/usr/bin/env python3
"""Independent raw replay, no collector imports and no codec execution."""
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys

PROFILES = ('native','auto','gfni')
ORDERS = [('auto_over_native',('native','auto','auto','native')),
          ('gfni_over_auto',('auto','gfni','gfni','auto')),
          ('gfni_over_native',('native','gfni','gfni','native')),
          *[('same_'+p,(p,)*4) for p in PROFILES]]
CELLS = [(1000,199,32768)]



def require(value,message):
    if not value: raise ValueError(message)


def same(a,b):
    require(json.dumps(a,sort_keys=True,allow_nan=False) ==
            json.dumps(b,sort_keys=True,allow_nan=False),'value/type mismatch')


def sha(path):
    with path.open('rb') as stream: return hashlib.file_digest(stream,'sha256').hexdigest()


def read(path): return json.loads(path.read_text())


def record(value,expected,measured,calls=None):
    same({k:v for k,v in value.items() if k not in ('encode_calls','samples_ns')},expected)
    same(value['encode_calls'],calls if calls is not None else 26 if measured else 1)
    values=value['samples_ns']
    require(type(values) is list and len(values)==(21 if measured else 0) and
            all(type(v) is int and v>0 for v in values),'sample list')


def derive(rows,expected):
    require(type(rows) is list and len(rows)==72,'incomplete')
    cells=[]
    index=0
    for cell in range(1):
        rounds={name:[] for name,_ in ORDERS}
        for round_id in range(3):
            for comparison,order in ORDERS:
                medians=[]
                for slot,profile in enumerate(order):
                    row=rows[index]
                    index+=1
                    same({k:row[k] for k in ('cell','round','comparison','slot','profile','sibling_delta')},
                         dict(cell=cell,round=round_id,comparison=comparison,slot=slot,
                              profile=profile,sibling_delta=0))
                    record(row['record'],expected[profile][cell],True)
                    medians.append(sorted(row['record']['samples_ns'])[10])
                rounds[comparison].append(math.sqrt(medians[0]/medians[1]*medians[3]/medians[2]))
        cells.append(dict(cell=cell,round_ratios=rounds,ratios={k:math.exp(sum(map(math.log,v))/3)
                                                             for k,v in rounds.items()}))
    controls=all(1/1.02 <= c['ratios']['same_'+p] <= 1.02 for c in cells for p in PROFILES)
    gain=cells[0]['ratios']['gfni_over_auto']>=1.05 and all(
        v>1 for v in cells[0]['round_ratios']['gfni_over_auto'])
    return dict(cells=cells,controls_pass=controls,
                decision='inconclusive_controls' if not controls else
                    'qualify_bounded_auto_candidate' if gain else 'reject_for_this_screen',
                production_promotion=False,confidence_intervals=False,
                neighbor_qualification=False,authoritative_v19=False)


def compare_analysis(a,b):
    same({k:v for k,v in a.items() if k!='cells'},{k:v for k,v in b.items() if k!='cells'})
    require(len(a['cells'])==len(b['cells'])==1,'cell count')
    for x,y in zip(a['cells'],b['cells']):
        same(x['cell'],y['cell'])
        for category in ('ratios','round_ratios'):
            same(sorted(x[category]),sorted(y[category]))
            for key in x[category]:
                left=x[category][key] if category=='round_ratios' else [x[category][key]]
                right=y[category][key] if category=='round_ratios' else [y[category][key]]
                require(len(left)==len(right),'ratio count')
                require(all(type(v) is float and math.isfinite(v) for v in left),'ratio type')
                require(all(math.isclose(v,w,rel_tol=1e-14,abs_tol=0) for v,w in zip(left,right)),
                        'ratio mismatch')


def resource(path,limit):
    lines=path.read_text().splitlines()
    require(lines.count('memory.peak')==1,'scope count')
    i=lines.index('memory.peak')
    peak=int(lines[i+1])
    require(0<peak<=limit,'memory peak')
    same(lines[i+2:],['memory.max',str(limit),'memory.events','low 0','high 0','max 0',
                     'oom 0','oom_kill 0','oom_group_kill 0','memory.swap.current','0','memory.swap.max','0'])
    require('\tExit status: 0' in lines,'scope status')
    return peak



QUALIFICATION = Path('/home/catid/leopard/.research/leopard-79h/r199-boundary-qualified.s3qjbngr')
NORMALIZATION = Path('/home/catid/leopard/.research/leopard-79h/r199-boundary-delivery.F0mkol')
SOURCE_FILES = {'r199_boundary_screen.cpp','run_r199_boundary_screen.py',
    'run_split_cache_screen.py','test_r199_boundary_screen.py','replay_r199_boundary_screen.py',
    'test_replay_r199_boundary_screen.py','verify_r199_boundary_checks.py',
    'verify_gf16_callback_probe.py','verify_gfni_source_stage.py','audit_avx2_pair_schedule.py',
    'test_verify_r199_boundary_checks.py'}
FROZEN_FILES = SOURCE_FILES | {'native','current','native.a','current.a','expected.json',
    'r199_boundary_screen_plan.json','build.json','preflight.json','qualification.json',
    'normalization.json','check.sh','check-condition.sh'}


def inputs(root,commit=None):
    from verify_r199_boundary_checks import replay as qualification_replay, public_record
    frozen=root/'frozen'
    pins=read(frozen/'pins.json')
    same(sorted(pins),['files','schema'])
    same(pins['schema'],'leopard-r199-boundary-pins/v1')
    same(sorted(pins['files']),sorted(FROZEN_FILES))
    for name,digest in pins['files'].items():
        path=frozen/name
        require(path.is_file() and not path.is_symlink() and not path.stat().st_mode & 0o222,'mutable input')
        same(sha(path),digest)
    plan=read(frozen/'r199_boundary_screen_plan.json')
    same(plan['schema'],'leopard-r199-boundary-plan/v1')
    same(plan['bead'],'leopard-79h.38.5.4.19')
    same([(c['k'],c['r'],c['bytes']) for c in plan['cells']],CELLS)
    same(plan['orders'],{k:list(v) for k,v in ORDERS})
    same([plan[k] for k in ('cpu','sibling','controller_cpu','rounds','samples_per_process',
        'attempt_budget','passive_seconds','control_bound','future_candidate_minimum_gain',
        'production_promotion','neighbor_qualification','confidence_intervals')],
        [26,90,0,3,21,1,10,1.02,1.05,False,False,False])
    same(sorted(plan['artifact_sha256']),['current','current.a','expected.json','native','native.a'])
    for name,digest in plan['artifact_sha256'].items(): same(digest,pins['files'][name])
    for bundle,key,digest in (
        (QUALIFICATION,'qualification_manifest_sha256','a509866ab0dbf2967a8a25258816134927f4b29ac8ce50c7bb3c2ce54de2fad3'),
        (NORMALIZATION,'normalization_manifest_sha256','1746a82fe079ca4327717a5917e770cc403940219acb0db22c98bd91759e6037')):
        same(plan[key],digest); same(sha(bundle/'SHA256SUMS'),digest)
    qualification=qualification_replay(QUALIFICATION)
    same(read(frozen/'qualification.json'),qualification)
    build=read(frozen/'build.json')
    prep=read(frozen/'preflight.json')
    same(build,read(QUALIFICATION/'build/build.json'))
    same(prep,read(QUALIFICATION/'preparation/preparation.json'))
    normal=read(frozen/'normalization.json')
    same(normal,read(NORMALIZATION/'normalized-objects-samepath/report.json'))
    same(normal['completed'],True)
    same(normal['normalized_source_sha256'],sha(frozen/'r199_boundary_screen.cpp'))
    same(normal['compiled_source_sha256'],sha(frozen/'r199_boundary_screen.cpp'))
    require((frozen/'r199_boundary_screen.cpp').read_bytes()+b'\n' ==
            (QUALIFICATION/'build/drivers/r199_boundary_screen.cpp').read_bytes(),'EOF-only source normalization')
    same([r['profile'] for r in normal['records']],['native','release','sanitize'])
    for row in normal['records']:
        name=row['profile']
        same(sha(NORMALIZATION/'normalized-objects-samepath'/(name+'.o')),row['object_sha256'])
        same(row['object_sha256'],build['profiles'][name]['files']['screen.o'])
    expected=read(frozen/'expected.json')
    same(sorted(expected),sorted(PROFILES))
    for variant in PROFILES:
        profile='native' if variant=='native' else 'release'
        name='native' if variant=='native' else 'current'
        same(sha(frozen/name),build['profiles'][profile]['files']['screen'])
        same(sha(frozen/(name+'.a')),build['profiles'][profile]['archive_sha256'])
        same(len(expected[variant]),1)
        value=dict(expected[variant][0],encode_calls=1,samples_ns=[])
        public_record(value,profile,variant)
        same(value,prep['public_records'][profile+'-'+variant+'-check'])
    if commit is not None:
        require(len(commit)==40 and all(c in '0123456789abcdef' for c in commit),'commit')
        for name in SOURCE_FILES | {'r199_boundary_screen_plan.json'}:
            require(subprocess.check_output(['git','show',commit+':experiments/leopard2/gf16_high_encode/'+name]) ==
                    (frozen/name).read_bytes(),'preregistered source: '+name)
    return plan,pins,expected,qualification


def replay(root,commit):
    plan,pins,expected,qualification=inputs(root,commit)
    state=read(root/'attempt1/attempt.json')
    same(state['schema'],'leopard-r199-boundary-attempt/v1')
    same(state['preregistration'],commit)
    same(state['plan_sha256'],pins['files']['r199_boundary_screen_plan.json'])
    same(state['pins'],pins)
    same(state['host'],plan['host'])
    same(state['complete'],True)
    require('failure' not in state,'attempt failure')
    passive=state['passive']
    same(passive['before'],passive['after'])
    require(type(passive['before']) is int and passive['before']>=0,'passive ticks')
    require(type(passive['elapsed_ns']) is int and passive['elapsed_ns']>=10000000000,'passive duration')
    same(len(state['preflight']),3)
    labels=['condition-before','condition-after']
    for variant,value in zip(PROFILES,state['preflight']):
        record(value,expected[variant][0],False)
        label='check-0-'+variant
        same(read(root/'attempt1'/(label+'.stdout')),value)
        labels.append(label)
    derived=derive(state['invocations'],expected)
    compare_analysis(state['analysis'],derived)
    for row in state['invocations']:
        label=f"cell-{row['cell']}-round-{row['round']}-{row['comparison']}-slot-{row['slot']}"
        same(read(root/'attempt1'/(label+'.stdout')),row['record'])
        labels.append(label)
    for suffix in ('stdout','stderr'):
        same(sorted(p.name for p in (root/'attempt1').glob('*.'+suffix)),sorted(l+'.'+suffix for l in labels))
    for path in (root/'attempt1').glob('*.stderr'): same(path.stat().st_size,0)
    condition=[]
    for unit in ('slipgate-catstream.service','slipgate-obs.service','slipgate-headless-x.service'):
        condition+=['MainPID=0','Id='+unit,'ActiveState=inactive','UnitFileState=disabled']
    for container in ('3fad3dfb247b3d4f3991e0a6801fba0805aae0901bdc89018841d4707a2d888c',
                      '002ba0a05bf860b12f2d799b49da1976aed663cb1494b2f8236751e18dc50182'):
        condition.append(container+' exited|false|no')
    for name in ('condition-before','condition-after'):
        same((root/'attempt1'/(name+'.stdout')).read_text().splitlines(),condition)
    peaks={name:resource(root/(name+'.log'),268435456)
           for name in ('pure-tests','pure-tests-opt','freeze','preclock','preclock-opt','attempt1')}
    return dict(preregistration=commit,frozen_inputs=len(pins['files']),preflights=3,
        timed_invocations=72,sibling_nonidle_jiffies=0,round_ratios=18,aggregate_ratios=6,
        preparation_checks=qualification['positive_native_records'],
        full_parity_comparison_bytes=qualification['full_native_leopard1_parity_bytes'],
        memory_peaks=peaks,analysis=derived)


if __name__=='__main__':
    require(len(sys.argv)==3,'usage: replay_r199_boundary_screen.py ROOT PREREGISTERED_COMMIT')
    print(json.dumps(replay(Path(sys.argv[1]),sys.argv[2]),indent=2))

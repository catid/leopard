#!/usr/bin/env python3
"""Independent raw replay, no collector imports and no codec execution."""
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys

PROFILES = ('native','pure','current')
ORDERS = [('pure_over_native',('native','pure','pure','native')),
          ('current_over_native',('native','current','current','native')),
          ('current_over_pure',('pure','current','current','pure')),
          *[('same_'+p,(p,)*4) for p in PROFILES]]
CELLS = [(1000,200,65536),(1000,200,32768),(1000,199,65536),(4096,512,4096)]


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
    require(type(rows) is list and len(rows)==288,'incomplete')
    cells=[]
    index=0
    for cell in range(4):
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
    return dict(cells=cells,controls_pass=controls,
                decision='attribution_available' if controls else 'inconclusive_controls',
                production_promotion=False,confidence_intervals=False,authoritative_v19=False)


def compare_analysis(a,b):
    same({k:v for k,v in a.items() if k!='cells'},{k:v for k,v in b.items() if k!='cells'})
    require(len(a['cells'])==len(b['cells'])==4,'cell count')
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


def replay(root,commit):
    frozen=root/'frozen'
    plan=read(frozen/'avx2_isa_screen_plan.json')
    pins=read(frozen/'pins.json')
    expected=read(frozen/'expected.json')
    state=read(root/'attempt1/attempt.json')
    same(state['preregistration'],commit)
    same(state['plan_sha256'],sha(frozen/'avx2_isa_screen_plan.json'))
    same(state['pins'],pins)
    same(state['host'],plan['host'])
    same(state['complete'],True)
    require('failure' not in state,'attempt failure')
    for name,digest in pins['files'].items():
        require(Path(name).name==name,'unsafe pin')
        path=frozen/name
        require(path.is_file() and not path.is_symlink() and not path.stat().st_mode & 0o222,'mutable pin')
        same(sha(path),digest)
    for name,digest in plan['artifact_sha256'].items(): same(digest,pins['files'][name])
    for name in ('avx2_isa_screen_plan.json','run_avx2_isa_screen.py','run_split_cache_screen.py',
                 'avx2_isa_screen.cpp','test_avx2_isa_screen.py','replay_avx2_isa_screen.py'):
        require(subprocess.check_output(['git','show',commit+':experiments/leopard2/gf16_high_encode/'+name])==
                (frozen/name).read_bytes(),'preregistered source')
    same([(c['k'],c['r'],c['bytes']) for c in plan['cells']],CELLS)
    same(plan['orders'],{k:list(v) for k,v in ORDERS})
    same([plan[k] for k in ('cpu','sibling','rounds','samples_per_process','attempt_budget','passive_seconds')],
         [26,90,3,21,1,10])
    passive=state['passive']
    same(passive['before'],passive['after'])
    require(type(passive['elapsed_ns']) is int and passive['elapsed_ns']>=10000000000,'passive duration')
    require(len(state['preflight'])==12,'preflight count')
    for index,value in enumerate(state['preflight']):
        cell,profile=index//3,PROFILES[index%3]
        record(value,expected[profile][cell],False)
        same(read(root/'attempt1'/f'check-{cell}-{profile}.stdout'),value)
    for row in state['invocations']:
        label=f"cell-{row['cell']}-round-{row['round']}-{row['comparison']}-slot-{row['slot']}"
        same(read(root/'attempt1'/(label+'.stdout')),row['record'])
    for path in (root/'attempt1').glob('*.stderr'): require(path.stat().st_size==0,'native stderr')
    for suffix in ('stdout','stderr'):
        require(len(list((root/'attempt1').glob('*.'+suffix)))==302,'raw record inventory')
    condition=[]
    for unit in ('slipgate-catstream.service','slipgate-obs.service','slipgate-headless-x.service'):
        condition+=['MainPID=0','Id='+unit,'ActiveState=inactive','UnitFileState=disabled']
    for container in ('3fad3dfb247b3d4f3991e0a6801fba0805aae0901bdc89018841d4707a2d888c',
                      '002ba0a05bf860b12f2d799b49da1976aed663cb1494b2f8236751e18dc50182'):
        condition.append(container+' exited|false|no')
    for name in ('condition-before','condition-after'):
        same((root/'attempt1'/(name+'.stdout')).read_text().splitlines(),condition)
    derived=derive(state['invocations'],expected)
    compare_analysis(state['analysis'],derived)
    preparation=read(root/'preflight.json')
    same(preparation['complete'],True)
    same(preparation['timings'],False)
    same([len(preparation[k]) for k in ('records','full_comparisons','clock_guard_rejections','cli_rejections')],
         [48,16,4,32])
    build=read(root/'build.json')
    for name,digest in build['files'].items(): same(sha(root/'build'/name),digest)
    same(build['source_sha256'],sha(frozen/'avx2_isa_screen.cpp'))
    for profile in PROFILES:
        same(sha(root/'build'/profile),sha(frozen/profile))
        same(sha(root/'build'/(profile+'.a')),sha(frozen/(profile+'.a')))
    repo=Path(subprocess.check_output(['git','rev-parse','--show-toplevel'],text=True).strip())
    reference=repo/'.research/leopard-79h/gf16-current-route-failed.STc10h'
    same(sha(reference/'SHA256SUMS'),'e83a217e680a23088e53208d9ba7747742e7918a1377e74efe9c35399b4ba342')
    oldpins={str(Path(name)):digest for digest,name in (line.split('  ',1)
             for line in (reference/'SHA256SUMS').read_text().splitlines())}
    prep_records={r['label']:r for r in preparation['records']}
    require(len(prep_records)==48,'duplicate preparation records')
    parity_bytes=0
    for profile in (*PROFILES,'sanitize'):
        for cell,(k,r,size) in enumerate(CELLS):
            identity=expected['current' if profile=='sanitize' else profile][cell]
            for label,calls in ((f'{profile}-{cell}',1),(f'{profile}--check-{cell}',1),
                                (f'{profile}--exercise-{cell}',26)):
                value=read(root/'preflight'/(label+'.stdout'))
                same(value,prep_records[label]['result'])
                record(value,identity,False,calls)
                require((root/'preflight'/(label+'.stderr')).stat().st_size==0,'preparation stderr')
            old=f'preflight/{(0,3,4,5)[cell]}-main.parity'
            same(sha(reference/old),oldpins[old])
            actual=root/'preflight'/f'{profile}-{cell}.parity'
            require(actual.stat().st_size==r*size,'parity size')
            with actual.open('rb') as left,(reference/old).open('rb') as right:
                while True:
                    chunk=left.read(65536)
                    require(chunk==right.read(65536),'full parity mismatch')
                    if not chunk: break
            parity_bytes+=r*size
        require((root/'preflight'/(profile+'-clock-selfcheck.stdout')).stat().st_size==0,'guard stdout')
        same((root/'preflight'/(profile+'-clock-selfcheck.stderr')).read_text(),
             'unexpected driver benchmark clock\n')
    same(parity_bytes,139198464)
    peaks={name:resource(root/(name+'.log'),536870912 if name=='build' else 268435456)
           for name in ('build','preflight','pure-tests','pure-tests-opt','attempt1')}
    return dict(preregistration=commit,frozen_inputs=len(pins['files']),preflights=12,timed_invocations=288,
                sibling_nonidle_jiffies=0,round_ratios=72,aggregate_ratios=24,
                preparation_checks=48,full_parity_comparison_bytes=parity_bytes,
                clock_guard_rejections=4,cli_rejections=32,
                memory_peaks=peaks,analysis=derived)


if __name__=='__main__':
    require(len(sys.argv)==3,'usage: replay_avx2_isa_screen.py ROOT PREREGISTERED_COMMIT')
    print(json.dumps(replay(Path(sys.argv[1]),sys.argv[2]),sort_keys=True))

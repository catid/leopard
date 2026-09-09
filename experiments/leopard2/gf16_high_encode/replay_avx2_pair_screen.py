#!/usr/bin/env python3
"""Independent stdlib replay: no collector imports or codec execution."""
import hashlib
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys

ORDERS = [('on_over_off',('off','on','on','off')),('same_off',('off',)*4),
          ('same_on',('on',)*4),('on_over_native',('native','on','on','native')),
          ('same_native',('native',)*4)]
CELLS = [(1000,200,32768,2,3),(1000,199,65536,2,3),(1000,200,65536,2,3),
         (4096,512,4096,2,3),(1000,199,32768,2,0),(1000,200,32768,2,6),
         (1000,200,32768,2,0),(17,7,64,1,3)]


def require(value, message):
    if not value: raise ValueError(message)


def same(a, b):
    require(json.dumps(a,sort_keys=True,allow_nan=False)==json.dumps(b,sort_keys=True,allow_nan=False),
            'value/type mismatch')


def read(path): return json.loads(path.read_text())


def sha(path):
    with path.open('rb') as stream: return hashlib.file_digest(stream,'sha256').hexdigest()


def identity(value):
    return {k:v for k,v in value.items() if k not in ('encode_calls','samples_ns','route_counts')}


def record(value, expected, measured):
    same(identity(value),expected)
    same(value['traced'],False)
    same(value['route_counts'],[0,0,0,0])
    same(value['encode_calls'],26 if measured else 1)
    values=value['samples_ns']
    require(type(values) is list and len(values)==(21 if measured else 0) and
            all(type(v) is int and v>0 for v in values),'sample list')


def derive(rows, expected):
    require(type(rows) is list and len(rows)==360,'incomplete')
    cells=[]
    index=0
    for cell in range(8):
        orders=ORDERS if cell<3 else ORDERS[:3]
        rounds={name:[] for name,_ in orders}
        for round_id in range(3):
            for comparison,order in orders:
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
    controls=all(1/1.02 <= value <= 1.02 for c in cells for name,value in c['ratios'].items()
                 if name.startswith('same_'))
    unchanged=all(1/1.02 <= cells[i]['ratios']['on_over_off'] <= 1.02 for i in (5,6,7))
    neighbors=all(cells[i]['ratios']['on_over_off'] >= 1/1.02 for i in (3,4))
    targets=all(cells[i]['ratios']['on_over_off'] >= 1.05 for i in (0,1,2))
    decision='inconclusive_controls' if not controls or not unchanged else \
        'reject_neighbor_regression' if not neighbors else 'candidate_pass' if targets else 'below_target_gate'
    return dict(cells=cells,controls_pass=controls,unchanged_neighbors_pass=unchanged,
                affected_neighbors_pass=neighbors,targets_pass=targets,decision=decision,
                production_promotion=False,confidence_intervals=False,authoritative_v19=False)


def compare_analysis(a,b):
    same({k:v for k,v in a.items() if k!='cells'},{k:v for k,v in b.items() if k!='cells'})
    require(len(a['cells'])==len(b['cells'])==8,'cell count')
    for x,y in zip(a['cells'],b['cells']):
        same(x['cell'],y['cell'])
        for category in ('ratios','round_ratios'):
            same(sorted(x[category]),sorted(y[category]))
            for key in x[category]:
                left=x[category][key] if category=='round_ratios' else [x[category][key]]
                right=y[category][key] if category=='round_ratios' else [y[category][key]]
                require(len(left)==len(right),'ratio count')
                require(all(type(v) is float and math.isfinite(v) for v in left),'ratio type')
                require(all(math.isclose(v,w,rel_tol=1e-14) for v,w in zip(left,right)),'ratio mismatch')


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


def preparation(root):
    prep=read(root/'preparation/preparation.json')
    same(prep['completed'],True)
    same(prep['timed'],False)
    original_root=Path(next(iter(prep['pins']))).parent.parent
    def locate(value):
        path=Path(value)
        return root/path.relative_to(original_root) if path.is_relative_to(original_root) else path
    for path,digest in prep['pins'].items(): same(sha(locate(path)),digest)
    require(len(prep['records'])==211,'preparation record count')
    labels=[r['label'] for r in prep['records']]
    require(len(set(labels))==211,'duplicate preparation record')
    positives=0
    for row in prep['records']:
        label=row['label']
        code=86 if label in ('native-native-clock-refusal','release-off-clock-refusal','release-on-clock-refusal') \
             else 1 if 'clock-refusal' in label or '-bad-' in label else 0
        same([row['returncode'],row['expected']],[code,code])
        for ext in ('stdout','stderr'): same(sha(root/'preparation'/(label+'.'+ext)),row[ext+'_sha256'])
        stderr=(root/'preparation'/(label+'.stderr')).read_text()
        if code==0:
            require(not stderr,'preparation stderr')
            positives+=1
        elif code==86: same(stderr,'unexpected driver benchmark clock\n')
    same(positives,184)
    for profile in ('release','sanitize'):
        for mode in ('off','on'):
            prefix=root/'preparation'/(profile+'-'+mode)
            same(Path(str(prefix)+'-pairs.stdout').read_text(),'pair kernel cases: 66147\n')
            same(Path(str(prefix)+'-split.stdout').read_text(),'split range cases: 256\n')
            same(Path(str(prefix)+'-roundtrip.stdout').read_text(),
                 ''.join('roundtrip cell '+str(i)+' passed\n' for i in range(3)))
            concurrent=Path(str(prefix)+'-concurrent.stdout').read_text().splitlines()
            same([concurrent.count('roundtrip cell 1 passed'),concurrent.count('roundtrip cell 2 passed'),
                  concurrent.count('four-thread both-field roundtrips passed'),len(concurrent)],[9,8,1,18])
            for cell in range(8):
                guarded=read(Path(str(prefix)+'-'+str(cell)+'.stdout'))
                same([guarded['cell'],guarded['subset_masks'],guarded['timed']],[cell,6,False])
    for profile in ('native','release','trace','sanitize'):
        modes=('native',) if profile=='native' else ('off','on')
        for mode in modes:
            for cell in range(8):
                prefix=root/'preparation'/(profile+'-'+mode+'-'+str(cell))
                check=read(Path(str(prefix)+'-check.stdout'))
                exercise=read(Path(str(prefix)+'-exercise.stdout'))
                same(identity(check),prep['expected'][profile][mode][cell])
                same(identity(exercise),identity(check))
                same([check['encode_calls'],exercise['encode_calls']],[1,26])
                same([check['samples_ns'],exercise['samples_ns']],[[],[]])
                same(check['route_counts'],prep['routes'][profile][mode][cell])
                same(exercise['route_counts'],[26*v for v in check['route_counts']])
                if profile in ('native','release'): same(read(Path(str(prefix)+'-plain.stdout')),check)
    for profile in ('trace','sanitize'):
        for cell in range(8):
            off=prep['routes'][profile]['off'][cell]
            on=prep['routes'][profile]['on'][cell]
            same(off,[on[1],on[0],on[3],on[2]])
            same(off,prep['routes']['trace']['off'][cell])
            if cell<5: require(off[0]>0 and off[2]>0 and off[1]==off[3]==0,'unobserved mode')
            else: same(off,[0,0,0,0])
    total=0
    require(len(prep['comparisons'])==52,'comparison inventory')
    for item in prep['comparisons']:
        a,b=locate(item['path']),locate(item['original'])
        same(sha(a),item['sha256'])
        same(sha(b),item['original_sha256'])
        same([a.stat().st_size,b.stat().st_size],[item['bytes']]*2)
        with a.open('rb') as left,b.open('rb') as right:
            while True:
                chunk=left.read(65536)
                require(chunk==right.read(65536),'full parity comparison')
                if not chunk: break
                total+=len(chunk)
            # Reclaim only this experiment's owned parity copies, never host
            # caches globally or another worker's mutable inputs.
            if a.is_relative_to(root): os.posix_fadvise(left.fileno(),0,0,os.POSIX_FADV_DONTNEED)
            if b.is_relative_to(root): os.posix_fadvise(right.fileno(),0,0,os.POSIX_FADV_DONTNEED)
    same(total,361368192)
    same(total,prep['parity_bytes'])
    # The fresh native driver must still agree with the independently retained
    # original oracle; four comparisons above bind old cells 3,4,0,5.
    build=read(root/'build/build.json')
    for name,digest in build['source_pins'].items(): same(sha(root/'source'/name),digest)
    for profile,entry in build['profiles'].items():
        archive=root/'build'/profile/'candidate.a'
        original=Path(entry['original'])
        same(sha(archive),entry['archive_sha256'])
        same(sha(original),entry['original_sha256'])
        before=subprocess.check_output(['ar','t',str(original)],text=True).splitlines()
        after=subprocess.check_output(['ar','t',str(archive)],text=True).splitlines()
        same(after,before+['avx2_pair_control.cpp.o'])
        same(len(entry['unchanged_members']),23)
        for member,digest in entry['unchanged_members'].items():
            old=subprocess.check_output(['ar','p',str(original),member])
            new=subprocess.check_output(['ar','p',str(archive),member])
            require(old==new,'unrelated archive member changed')
            same(hashlib.sha256(new).hexdigest(),digest)
    return dict(positive_records=positives,records=211,full_parity_bytes=total,
                route_counts=prep['routes']['trace'])


def replay(root,commit):
    frozen=root/'frozen'
    plan=read(frozen/'avx2_pair_screen_plan.json')
    pins=read(frozen/'pins.json')
    expected=read(frozen/'expected.json')
    state=read(root/'attempt1/attempt.json')
    same(state['preregistration'],commit)
    same(state['plan_sha256'],sha(frozen/'avx2_pair_screen_plan.json'))
    same(state['pins'],pins)
    same(state['host'],plan['host'])
    same(state['complete'],True)
    require('failure' not in state,'attempt failure')
    for name,digest in pins['files'].items():
        require(Path(name).name==name,'unsafe pin name')
        path=frozen/name
        require(path.is_file() and not path.is_symlink() and not path.stat().st_mode & 0o222,'mutable pin')
        same(sha(path),digest)
    for name,digest in plan['artifact_sha256'].items(): same(digest,pins['files'][name])
    for name in ('avx2_pair_screen_plan.json','run_avx2_pair_screen.py','run_split_cache_screen.py',
                 'avx2_pair_screen.cpp','test_avx2_pair_screen.py','replay_avx2_pair_screen.py',
                 'avx2_pair_control.cpp','avx2_pair_control.h','avx2_pair_schedule.h','avx2_pair_runtime.patch'):
        require(subprocess.check_output(['git','show',commit+':experiments/leopard2/gf16_high_encode/'+name])==
                (frozen/name).read_bytes(),'preregistered source')
    same([(c['k'],c['r'],c['bytes'],c['field'],c['backend']) for c in plan['cells']],CELLS)
    same(plan['orders'],{k:list(v) for k,v in ORDERS})
    same([plan[k] for k in ('cpu','sibling','rounds','samples_per_process','attempt_budget','passive_seconds')],
         [26,90,3,21,1,10])
    same([plan[k] for k in ('target_cells','affected_neighbors','unchanged_neighbors',
                           'minimum_target_gain','control_bound','affected_neighbor_minimum')],
         [[0,1,2],[3,4],[5,6,7],1.05,1.02,1/1.02])
    passive=state['passive']
    same(passive['before'],passive['after'])
    require(type(passive['elapsed_ns']) is int and passive['elapsed_ns']>=10000000000,'passive duration')
    wanted=[(c,p) for c in range(8) for p in (('off','on','native') if c<3 else ('off','on'))]
    same(len(state['preflight']),19)
    for value,(cell,profile) in zip(state['preflight'],wanted):
        record(value,expected[profile][cell],False)
        same(read(root/'attempt1'/f'check-{cell}-{profile}.stdout'),value)
    for row in state['invocations']:
        label=f"cell-{row['cell']}-round-{row['round']}-{row['comparison']}-slot-{row['slot']}"
        same(read(root/'attempt1'/(label+'.stdout')),row['record'])
    for path in (root/'attempt1').glob('*.stderr'): require(path.stat().st_size==0,'native stderr')
    for suffix in ('stdout','stderr'): same(len(list((root/'attempt1').glob('*.'+suffix))),381)
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
    prep=preparation(root)
    peaks={name:resource(root/(name+'.log'),536870912 if name in ('build','drivers') else 268435456)
           for name in ('build','drivers','preparation','unit','unit-opt','attempt1')}
    return dict(preregistration=commit,frozen_inputs=len(pins['files']),preflights=19,timed_invocations=360,
                sibling_nonidle_jiffies=0,round_ratios=90,aggregate_ratios=30,
                preparation=prep,memory_peaks=peaks,analysis=derived)


if __name__=='__main__':
    require(len(sys.argv)==3,'usage: replay_avx2_pair_screen.py ROOT PREREGISTERED_COMMIT')
    print(json.dumps(replay(Path(sys.argv[1]),sys.argv[2]),sort_keys=True))

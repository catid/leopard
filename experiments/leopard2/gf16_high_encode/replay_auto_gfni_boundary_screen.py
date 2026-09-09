#!/usr/bin/env python3
"""Independent retained-only AUTO qualification replay; never imports the collector."""
import copy
import hashlib
import json
import math
import os
from pathlib import Path
import sys

COMMIT = '6bff9ecb8658c1313e88d680b47413d3f2531065'
SHAPES = ((1000,200,32768),(1000,199,65536),(1000,200,32768),(1000,199,65536),
          (1000,200,65536),(1000,199,32768),(1000,200,32768),(4096,512,4096))
ORDERS = (('off_on',('off','on','on','off')), ('same_on',('on_a','on_b','on_b','on_a')))


def require(value, message):
    if not value: raise ValueError(message)


def same(actual, expected):
    require(json.dumps(actual,sort_keys=True,allow_nan=False) ==
            json.dumps(expected,sort_keys=True,allow_nan=False), 'typed identity')


def sha(path):
    with path.open('rb') as stream:
        value = hashlib.file_digest(stream,'sha256').hexdigest()
    return value


def read(path):
    require(path.stat().st_size < 1048576,'bounded JSON')
    return json.loads(path.read_text())


def schedule():
    for cell in range(8):
        comparisons = ORDERS + (('main_on',('main','on','on','main')),) if cell < 2 else ORDERS
        for round_id in range(3):
            for comparison, order in comparisons:
                for slot, variant in enumerate(order):
                    yield dict(cell=cell,round=round_id,comparison=comparison,
                               slot=slot,variant=variant,sibling_delta=0)


def derive(rows, expected):
    require(type(rows) is list and len(rows) == 216,'complete ordered attempt')
    grouped = {}
    for row, wanted in zip(rows,schedule()):
        same({key:row[key] for key in wanted},wanted)
        record = row['record']
        variant = 'on' if wanted['variant'].startswith('on') else wanted['variant']
        same({key:value for key,value in record.items() if key != 'samples_ns'},
             expected[variant][wanted['cell']])
        values = record['samples_ns']
        require(type(values) is list and len(values) == 21 and
                all(type(x) is int and x > 0 for x in values),'positive integer samples')
        key = wanted['cell'],wanted['round'],wanted['comparison']
        grouped.setdefault(key,[]).append(sorted(values)[10])
    cells = []
    for cell in range(8):
        rounds = {}
        for name in ('off_on','same_on','main_on') if cell < 2 else ('off_on','same_on'):
            values = []
            for round_id in range(3):
                a,b,c,d = grouped[cell,round_id,name]
                values.append(math.sqrt((a*d)/(b*c)))
            rounds[name] = values
        cells.append(dict(cell=cell,role='target' if cell < 4 else 'unchanged_neighbor',
                          round_ratios=rounds,
                          ratios={name:math.prod(values)**(1/3) for name,values in rounds.items()}))
    controls = all(1/1.02 <= x['ratios']['same_on'] <= 1.02 for x in cells)
    neighbors = all(1/1.02 <= x['ratios']['off_on'] <= 1.02 for x in cells[4:])
    targets = all(x['ratios']['off_on'] >= 1.05 and
                  all(y > 1 for y in x['round_ratios']['off_on']) for x in cells[:4])
    decision = ('inconclusive_controls' if not controls else 'reject_neighbor_gate' if not neighbors
                else 'reject_target_gate' if not targets else 'continue_to_production_integration')
    return dict(cells=cells,controls_pass=controls,neighbors_pass=neighbors,targets_pass=targets,
                decision=decision,confidence_intervals=False,production_promotion=False,
                authoritative_v19=False)


def analysis_equal(actual, expected):
    require(set(actual) == set(expected),'analysis fields')
    for key in expected.keys()-{'cells'}: same(actual[key],expected[key])
    require(len(actual['cells']) == 8,'cell count')
    for a,b in zip(actual['cells'],expected['cells']):
        require(set(a) == set(b),'cell fields')
        same(a['cell'],b['cell']); same(a['role'],b['role'])
        for key in ('round_ratios','ratios'):
            require(set(a[key]) == set(b[key]),'comparison fields')
            for name in b[key]:
                left = a[key][name] if key == 'round_ratios' else [a[key][name]]
                right = b[key][name] if key == 'round_ratios' else [b[key][name]]
                require(len(left) == len(right) and all(type(x) is float and
                        math.isclose(x,y,rel_tol=2e-14,abs_tol=0) for x,y in zip(left,right)),
                        'independent ratio replay')


def resource(path, maximum):
    lines = path.read_text().splitlines()
    require(lines.count('memory.peak') == 1 and '\tExit status: 0' in lines,'scope exit')
    index = lines.index('memory.peak')
    peak = int(lines[index+1])
    require(0 < peak <= maximum,'memory peak')
    same(lines[index:],['memory.peak',str(peak),'memory.max',str(maximum),
        'memory.events','low 0','high 0','max 0','oom 0','oom_kill 0','oom_group_kill 0',
        'memory.swap.current','0','memory.swap.max','0'])
    return peak


def mutations(rows, expected, derived):
    rejected = 0
    for key,value in (('cell',1),('round',1),('slot',1),('variant','on'),
                      ('comparison','same_on'),('sibling_delta',1),('sibling_delta',False),
                      ('partial',None),('samples',[True]*21),('samples',[0]*21),
                      ('samples',[1]*20),('identity',0)):
        changed = copy.deepcopy(rows)
        if key == 'partial': changed.pop()
        elif key == 'samples': changed[0]['record']['samples_ns'] = value
        elif key == 'identity': changed[0]['record']['untimed_route_calls'] = True
        else: changed[0][key] = value
        try: derive(changed,expected)
        except ValueError: rejected += 1
    require(rejected == 12,'row mutation escaped')
    for key in ('production_promotion','authoritative_v19','confidence_intervals'):
        changed = copy.deepcopy(derived); changed[key] = True
        try: analysis_equal(changed,derived)
        except ValueError: rejected += 1
    require(rejected == 15,'claim mutation escaped')
    # Exact synthetic samples test all decision gates independently of the actual outcome.
    synthetic = copy.deepcopy(rows)
    for row in synthetic:
        row['record']['samples_ns'] = [140 if (row['cell'] < 4 and row['variant'] == 'off')
                                      or row['variant'] == 'main' else 100]*21
    require(derive(synthetic,expected)['decision'] == 'continue_to_production_integration','positive gate')
    for condition,decision in (('control','inconclusive_controls'),('neighbor','reject_neighbor_gate'),
                               ('target','reject_target_gate')):
        changed = copy.deepcopy(synthetic)
        for row in changed:
            if (condition == 'control' and row['cell'] == 7 and row['variant'] == 'on_a') or \
               (condition == 'neighbor' and row['cell'] == 5 and row['variant'] == 'off'):
                row['record']['samples_ns'] = [105]*21
            if condition == 'target' and row['cell'] == 3 and row['variant'] == 'off':
                row['record']['samples_ns'] = [104]*21
        require(derive(changed,expected)['decision'] == decision,'decision precedence')
    return rejected


def replay(root):
    for name,value in {
        'attempt1/attempt.json':'383790fd94221d107dc3b2fccc9396ac6404828a8811ea27a2f2b15fec3f8537',
        'attempt1.log':'5d644d0cde435b9453de3396af20335c987aaa16e114dfb7e727272cb7167a4e',
        'frozen/pins.json':'d5a53d20b603a498e1aa6cede047106769f0616b94acb14a8dd9c7a5af3888bf',
    }.items(): same(sha(root/name),value)
    pins = read(root/'frozen/pins.json')
    same(pins['source_commit'],COMMIT)
    require(len(pins['files']) == 14,'frozen inventory')
    for name,value in pins['files'].items():
        path = root/'frozen'/name
        require(Path(name).name == name and path.is_file() and not path.is_symlink() and
                not path.stat().st_mode & 0o222 and sha(path) == value,'frozen input: '+name)
    plan = read(root/'frozen/auto_gfni_boundary_screen_plan.json')
    state = read(root/'attempt1/attempt.json')
    expected = read(root/'frozen/expected.json')
    same(state['pins'],pins); same(state['host'],plan['host'])
    same(state['plan_sha256'],pins['files']['auto_gfni_boundary_screen_plan.json'])
    for name,value in plan['artifact_sha256'].items(): same(pins['files'][name],value)
    same(pins['files']['leopard2.cpp'],plan['core_sha256'])
    same(pins['files']['Leopard2Direct.h'],plan['header_sha256'])
    require(state['complete'] is True and 'failure' not in state,'complete attempt required')
    passive = state['passive']
    require(type(passive['before']) is int and passive['before'] >= 0 and
            type(passive['after']) is int and passive['before'] == passive['after'] and
            type(passive['elapsed_ns']) is int and passive['elapsed_ns'] >= 10000000000,'passive isolation')
    build_artifacts = 0
    for line in (root/'build.log').read_text().splitlines():
        words = line.split()
        if len(words) == 2 and len(words[0]) == 64 and \
                words[1].startswith('/tmp/leopard-auto-boundary-screen.2TfJyu/stage/'):
            same(sha(root/'stage'/Path(words[1]).name),words[0]); build_artifacts += 1
    require(build_artifacts == 7,'build artifact inventory')
    for name in ('main','current','main.a','current.a','expected.json','auto_gfni_boundary_screen.cpp'):
        same(sha(root/'stage'/name),pins['files'][name])
    for name in ('leopard2.cpp','Leopard2Direct.h'):
        same(sha(root/'codec-current'/name),pins['files'][name])
    require(len(state['preflight']) == 24,'preflight count')
    parity_bytes = 0
    for cell,(k,r,size) in enumerate(SHAPES):
        for index,variant in enumerate(('main','off','on')):
            identity = expected[variant][cell]
            same([identity[x] for x in ('cell','k','r','bytes')],[cell,k,r,size])
            same(identity['boundary_mode'],-1 if variant == 'main' else int(variant == 'on'))
            gfni = variant != 'main' and (cell == 4 or (cell < 4 and variant == 'on'))
            same(identity['execution_route'],'native-avx2' if variant == 'main' else 'gfni' if gfni else 'avx2')
            same(identity['untimed_route_calls'],int(gfni))
            same(identity['api'],'leo_encode' if variant == 'main' else
                 'leo2_encode_batch_one_item' if cell in (2,3) else 'leo2_encode')
            same(identity['codec_commit'],plan['main_commit'] if variant == 'main' else plan['codec_commit'])
            record = dict(identity,samples_ns=[])
            same(state['preflight'][cell*3+index],record)
            for folder,stem in (('attempt1',f'check-{cell}-{variant}'),('preflight',f'{cell}-{variant}')):
                same(read(root/folder/(stem+'.stdout')),record)
                require((root/folder/(stem+'.stderr')).stat().st_size == 0,'preflight stderr')
            if variant != 'main':
                same(read(root/'preflight'/f'{cell}-{variant}-sanitize.stdout'),record)
                require((root/'preflight'/f'{cell}-{variant}-sanitize.stderr').stat().st_size == 0,'sanitizer stderr')
                original = root/'preflight'/f'{cell}-main.parity'
                current = root/'preflight'/f'{cell}-{variant}.parity'
                require(original.stat().st_size == current.stat().st_size == r*size,'parity size')
                with original.open('rb') as a, current.open('rb') as b:
                    while True:
                        block = a.read(65536)
                        require(block == b.read(65536),'full Leopard1 parity')
                        if not block: break
                    for stream in (a,b): os.posix_fadvise(stream.fileno(),0,0,os.POSIX_FADV_DONTNEED)
                same(identity['input_hash'],expected['main'][cell]['input_hash'])
                same(identity['output_hash'],expected['main'][cell]['output_hash'])
                parity_bytes += r*size
    for row in state['invocations']:
        prefix = root/'attempt1'/f"cell-{row['cell']}-round-{row['round']}-{row['comparison']}-slot-{row['slot']}"
        same(read(prefix.with_suffix('.stdout')),row['record'])
        require(prefix.with_suffix('.stderr').stat().st_size == 0,'timed stderr')
    condition = []
    for unit in ('slipgate-catstream.service','slipgate-obs.service','slipgate-headless-x.service'):
        condition += ['MainPID=0','Id='+unit,'ActiveState=inactive','UnitFileState=disabled']
    for container in ('3fad3dfb247b3d4f3991e0a6801fba0805aae0901bdc89018841d4707a2d888c',
                      '002ba0a05bf860b12f2d799b49da1976aed663cb1494b2f8236751e18dc50182'):
        condition.append(container+' exited|false|no')
    for name in ('condition-before','condition-after'):
        same((root/'attempt1'/(name+'.stdout')).read_text().splitlines(),condition)
        require((root/'attempt1'/(name+'.stderr')).stat().st_size == 0,'condition stderr')
    for suffix in ('stdout','stderr'):
        require(len(list((root/'attempt1').glob('*.'+suffix))) == 242,'raw invocation count')
        require(len(list((root/'preflight').glob('*.'+suffix))) == 40,'preparation invocation count')
    peaks = {name:resource(root/(name+'.log'),536870912 if name == 'build' else 268435456)
             for name in ('build','preflight','pure-tests','pure-tests-opt','freeze','frozen-tests','bind-pins','attempt1')}
    derived = derive(state['invocations'],expected)
    analysis_equal(state['analysis'],derived)
    return dict(preregistration=COMMIT,frozen_inputs=14,build_artifacts=build_artifacts,
                preparation_checks=40,preflights=24,timed_invocations=216,
                sibling_nonidle_jiffies=0,full_parity_comparison_bytes=parity_bytes,
                rejected_mutations=mutations(state['invocations'],expected,derived),
                decision_gate_selfchecks=4,memory_peaks=peaks,analysis=derived)


if __name__ == '__main__':
    require(len(sys.argv) == 2,'usage: replay_auto_gfni_boundary_screen.py ROOT')
    print(json.dumps(replay(Path(sys.argv[1])),sort_keys=True))

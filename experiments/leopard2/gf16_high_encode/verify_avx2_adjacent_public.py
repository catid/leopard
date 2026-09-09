#!/usr/bin/env python3
"""Clock-free, collector-free public frontend replay and structural count oracle."""
import json
import os
from pathlib import Path
import sys

from build_avx2_adjacent_public import ARCHIVES, BEAD
from avx2_adjacent_counts import SHAPES, validate
from verify_gf16_callback_probe import model_shape
from verify_paired_r19932 import equal, require, read, parse, sha
from verify_auto_r19932_checks import scope

CELLS = [(1000,200,32768,2,3),(1000,199,65536,2,3),(1000,200,65536,2,3),
         (4096,512,4096,2,3),(1000,199,32768,2,0),(1000,200,32768,2,6),
         (1000,200,32768,2,0),(17,7,64,1,3),(1000,200,32768,2,3)]
INPUTS = ['8b78decd04e67d27','ba910426bb91dec4','ba910426bb91dec4','8ab693ec532fd3fa',
          '8b78decd04e67d27','8b78decd04e67d27','8b78decd04e67d27','bef540029ca6f65b','8b78decd04e67d27']
OUTPUTS = ['50d5c45deb48367a','0b985bca4df709d5','b1d946fa500bbf12','12c692cff473f0ee',
           '00648d9dbf0f2b20','50d5c45deb48367a','50d5c45deb48367a','83958df25bd2c455','50d5c45deb48367a']


def schedules(profile):
    require(profile in ARCHIVES, 'profile')
    return ('NNNN',) if profile == 'native' else ('PPPP',) if profile == 'original' else ('0110','1001','0000','1111')


def valid(profile, cell, schedule, group):
    require(type(cell) is int and 0 <= cell < 9, 'cell')
    require(schedule in schedules(profile), 'schedule')
    require(type(group) is int and (group == 1 or (cell == 7 and group == 256)), 'group')


def structural(cell):
    """Generate expected callback traversal, then apply the earlier edge model.

    This is independent of runtime observations and source-template codegen.
    Cell8 is one complete public encode in a one-item batch, hence shape0.
    """
    require(type(cell) is int and 0 <= cell < 9, 'cell')
    cell = 0 if cell == 8 else cell
    if cell == 7:
        row = dict(schema='gf16-callback-counts/v1', timed=False, calls=0, passes=[], buckets=[])
    else:
        k,r,b,kind,tile,passes = SHAPES[cell]
        traversal = model_shape(SHAPES[cell])
        row = dict(schema='gf16-callback-counts/v1', timed=False, calls=sum(traversal.values()),
            passes=[dict(kind=kind,k=k,r=r,requested=r,side=1 << (r-1).bit_length(),
                         sparse_blocks=0,bytes=tile,source_policy=b)]*passes,
            buckets=[dict(op=op,distance=d,zero_mask=z,prefer_fused=f,bytes=b,calls=n)
                     for (op,d,z,f,b),n in traversal.items()])
    return validate(row, cell)


def pair_counts(profile, cell, states):
    row = dict(calls=[[0,0],[0,0]], blocks=[[0,0],[0,0]])
    if profile not in ('trace','sanitize'): return row
    model = structural(cell)
    for state in states:
        for family, name in enumerate(('forward','accumulating')):
            row['calls'][family][int(state)] += model[name+'_pairs']
            row['blocks'][family][int(state)] += model[name+'_blocks']
    return row


def expected(profile, cell, schedule, group, mode, clock):
    valid(profile, cell, schedule, group)
    require(mode in ('--check','--exercise','--clock-exercise'), 'successful mode')
    require(clock in ('steady','synthetic','abort') and
            (mode != '--clock-exercise' or clock == 'synthetic'), 'clock')
    exercise, synthetic = mode != '--check', mode == '--clock-exercise'
    k,r,b,field,backend = CELLS[cell]
    scratch = ([16777216,33554432,33554432,4194304,16777216,16777216,16777216,1024,16777216]
               if profile == 'native' else [16808512]*3+[4308992]+[16808512]*3+[1728,16808512])[cell]
    calls = 1 + (25*group if exercise else 0)
    return dict(schema='leopard-adjacent-public/v1',codec=profile+':'+ARCHIVES[profile],
        cell=cell,k=k,r=r,bytes=b,field=field,backend=backend,
        api='leo_encode' if profile=='native' else 'leo2_encode_batch_one_item' if cell==8 else 'leo2_encode',
        schedule=schedule,group=group,warmup_passes=4 if exercise else 0,exercise_passes=21 if exercise else 0,
        encode_calls=4*calls,selections=104 if exercise else 4,per_slot_calls=[calls]*4,
        probes=[int(profile!='native' and cell==6)]*4,scratch_bytes=scratch,
        input_hash=INPUTS[cell],output_hash=OUTPUTS[cell],clock_source=clock,
        samples=[[257+17*i,(257+17*i)/group] for i in range(84)] if synthetic else [],
        traced=profile in ('trace','sanitize'),
        pair_probes=[pair_counts(profile,cell,s) for s in schedule[1:]],
        pair_totals=pair_counts(profile,cell,schedule*(25*group if exercise else 0)),
        timed=False,default_enabled=False)


def witness(profile, cell, schedule, group, passes=25, extra=0):
    valid(profile,cell,schedule,group)
    api = 2 if profile=='native' else 1 if cell==8 else 0
    states = [2 if s in 'NP' else int(s) for s in schedule]
    sequence = states + [s for _ in range(passes) for s in states for _ in range(group)] + [states[0]]*extra
    counts, apis, h = [0]*3, [0]*3, 14695981039346656037
    for state in sequence:
        counts[state] += 1; apis[api] += 1
        h = ((h ^ (state+4*api))*1099511628211) % 2**64
    return dict(schema='leopard-paired-witness/v1',calls=len(sequence),states=counts,apis=apis,order_hash=f'{h:016x}')


def clocks(group, samples=84):
    return dict(schema='paired-synthetic-clock/v1',clock_calls=2*samples,
        public_calls_at_clock=[4+16*group+(i//2+i%2)*group for i in range(2*samples)],timed=False)


def bad_arguments(p):
    s = schedules(p)[0]
    return [[],['--check'],['--timing','0',s,'1'],['--exercise','0',s,'1','a','b'],
            ['--check','9',s,'1'],['--check','00',s,'1'],['--check','0',s,'256'],
            ['--check','8',s,'256'],['--check','7',s,'0256'],['--check','7',s,'0'],
            ['--check','0','0101','1'],['--check','0','111','1'],
            ['--check','0','0110' if p in ('native','original') else 'PPPP','1'],
            ['--clock-exercise','0',s,'1'],['--clock-guard','0',s,'1','forbidden'],
            ['--clock-guard','7',s,'256']]


def inventory():
    rows = []
    for p in ARCHIVES:
        for c in range(9):
            for s in schedules(p):
                for g in ((1,256) if c==7 else (1,)):
                    for binary,mode in (('plain','--exercise'),('synthetic','--clock-exercise')):
                        rows.append(dict(label=f'{p}-{c}-{s}-{g}-{binary}',profile=p,binary=binary,
                            args=[mode,str(c),s,str(g)],code=0,parity=True,fault=None))
        s = schedules(p)[0]
        for c,g in ((0,1),(7,256)):
            rows.append(dict(label=f'{p}-{c}-check',profile=p,binary='plain',
                args=['--check',str(c),s,str(g)],code=0,parity=True,fault=None))
        for c,g in ((0,1),(7,256),(8,1)):
            rows.append(dict(label=f'{p}-{c}-no-clock',profile=p,binary='abort',
                args=['--exercise',str(c),s,str(g)],code=0,parity=True,fault=None))
            for order in schedules(p):
                rows.append(dict(label=f'{p}-{c}-{order}-abort',profile=p,binary='abort',
                    args=['--clock-guard',str(c),order,str(g)],code=86,parity=False,fault=None))
        for fault in ('equal','reverse','negative','huge'):
            rows.append(dict(label=f'{p}-fault-{fault}',profile=p,binary='synthetic',
                args=['--clock-exercise','7',s,'256'],code=1,parity=False,fault=fault))
        for i,args in enumerate(bad_arguments(p)):
            rows.append(dict(label=f'{p}-bad-{i}',profile=p,binary='plain',args=args,code=1,parity=False,fault=None))
        if p in ('release','sanitize'):
            rows.append(dict(label=p+'-unit',profile=p,binary='group-unit',args=[],code=0,parity=False,fault=None))
    return rows


def compare(left, right):
    total = 0
    with left.open('rb') as a, right.open('rb') as b:
        while True:
            chunk = a.read(65536)
            require(chunk == b.read(65536), 'full native parity differs')
            if not chunk: break
            total += len(chunk)
        for stream in (a,b): os.posix_fadvise(stream.fileno(),0,0,os.POSIX_FADV_DONTNEED)
    return total


def verify_record(row, output, error):
    lines = [parse(line) for line in output.splitlines()]
    actual = {r['schema']:r for r in lines}
    require(len(lines)==len(actual), 'duplicate output schema')
    wanted = []
    p,binary,args = row['profile'],row['binary'],row['args']
    if '-bad-' in row['label']:
        require(bool(error), 'CLI refusal')
    elif binary == 'group-unit':
        equal(error,''); wanted = [dict(schema='paired-group-unit/v1',cases=37,timed=False)]
    else:
        mode,c,s,g = args; c,g = int(c),int(g)
        if row['code']==0:
            equal(error,'')
            clock = 'steady' if binary=='plain' else binary
            wanted = [expected(p,c,s,g,mode,clock)]
            if binary!='plain': wanted += [witness(p,c,s,g)]
            if mode=='--clock-exercise': wanted += [clocks(g)]
        elif row['code']==86:
            equal(error,'unexpected driver benchmark clock\n')
            wanted = [witness(p,c,s,g,4)]
        else:
            equal(error,('group duration exceeds exact binary64 integer range' if row['fault']=='huge'
                         else 'nonpositive or reversed grouped clock interval')+'\n')
            wanted = [witness(p,c,s,g,4,g),clocks(g,1)]
    equal(actual,{r['schema']:r for r in wanted})


def replay(root):
    build,checks = read(root/'build/build.json'),read(root/'checks/checks.json')
    for state in (build,checks): equal([state['bead'],state['completed'],state['timed']],[BEAD,True,False])
    for name,digest in build['inputs'].items(): equal(sha(Path(name)),digest)
    for name,digest in build['artifacts'].items(): equal(sha(root/'build'/name),digest)
    for p,digest in ARCHIVES.items():
        equal(sha(root/'build'/p/'codec.a'),digest)
        common = 'native' if p=='native' else 'sanitize' if p=='sanitize' else 'l2-release'
        for binary in ('plain','synthetic','abort'):
            links = [cmd for cmd in build['commands'] if cmd[-1].endswith('/'+p+'/'+binary)]
            require(len(links)==1,'unique link')
            equal([a for a in links[0] if a.endswith('/driver.o')],
                  [str(Path(checks['root'])/'build'/(common+'-objects')/'driver.o')])
            symbols = (root/'build'/p/(binary+'-undefined.txt')).read_text()
            equal('_ZNSt6chrono3_V212steady_clock3nowEv' in symbols,binary=='plain')
    rows = inventory()
    equal([r['label'] for r in checks['records']],[r['label'] for r in rows])
    total = comparisons = 0
    for record,row in zip(checks['records'],rows):
        label = row['label']; folder = root/'checks'
        out,err = folder/(label+'.stdout'),folder/(label+'.stderr')
        equal(record['returncode'],row['code']); equal(record['fault'],row['fault'])
        equal(sha(out),record['stdout_sha256']); equal(sha(err),record['stderr_sha256'])
        args = list(row['args'])
        if row['parity']: args.append(str(Path(checks['root'])/'checks'/(label+'.parity')))
        equal(record['args'],[row['profile'],row['binary'],*args])
        verify_record(row,out.read_text(),err.read_text())
        if row['parity']:
            cell = int(row['args'][1]); parity = folder/(label+'.parity')
            equal(parity.stat().st_size,CELLS[cell][1]*CELLS[cell][2])
            equal(sha(parity),record['parity_sha256'])
            baseline = folder/f'native-{cell}-NNNN-1-plain.parity'
            if parity != baseline:
                total += compare(parity,baseline); comparisons += 1
    equal(sorted(p.name for p in (root/'checks').glob('*.stdout')), sorted(r['label']+'.stdout' for r in rows))
    resources = {phase:scope((root/(phase+'.log')).read_text(),limit,max_events=0)
                 for phase,limit in (('build',512*1024**2),('checks',256*1024**2))}
    return dict(bead=BEAD,positive=sum(r['code']==0 for r in rows),
        clock_aborts=sum(r['code']==86 for r in rows),clock_faults=sum(bool(r['fault']) for r in rows),
        cli_refusals=sum('-bad-' in r['label'] for r in rows),
        synthetic_spans=sum(r['parity'] and r['binary']=='synthetic' for r in rows)*84,
        parity_files=comparisons,parity_bytes=total,resources=resources,
        isolated_model=[structural(c) for c in range(9)],timed=False,default_enabled=False)


if __name__ == '__main__': print(json.dumps(replay(Path(sys.argv[1]).resolve()),sort_keys=True))

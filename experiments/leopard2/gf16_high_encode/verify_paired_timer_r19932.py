#!/usr/bin/env python3
"""Read-only replay, no codec or clock execution; leopard-79h.38.5.4.19.1.1."""
import json
from pathlib import Path
import sys

from verify_paired_r19932 import ARCHIVES, BEAD, CELLS, equal, expected, parse, read, sha, witness, compare, require
from verify_auto_r19932_checks import scope


def driver_expected(profile,cell,schedule,group,synthetic):
    row = expected(profile,cell,schedule,group,True)
    row['schema'] = 'leopard-paired-timer-r19932/v1'
    row['clock_source'] = 'synthetic' if synthetic else 'steady'
    row['samples'] = [[257+17*i,(257+17*i)/group] for i in range(84)] if synthetic else []
    return row


def clock_expected(group, samples=84):
    begin = 4+16*group
    calls = [begin+(i//2 + i%2)*group for i in range(2*samples)]
    return dict(schema='paired-synthetic-clock/v1',clock_calls=2*samples,public_calls_at_clock=calls,timed=False)


def fault_witness(profile,group):
    # Four probes + sixteen warmup groups + the first failing measured group.
    state = 2 if profile=='native' else 0
    api = 2 if profile=='native' else 0
    schedule = 'NNNN' if profile=='native' else '0110'
    row = witness(profile,8,schedule,group,'clock-guard')
    digest = int(row['order_hash'],16)
    for _ in range(group):
        digest = ((digest ^ (state+4*api))*1099511628211) % (2**64)
    row['calls'] += group; row['states'][state] += group; row['apis'][api] += group
    row['order_hash'] = f'{digest:016x}'
    return row


def bad_args(profile):
    s = 'NNNN' if profile=='native' else '0110'
    # Never invoke --measure under this qualification task, even as bad CLI.
    return [[],['--check'],['--timing','0',s,'1'],['--exercise','0',s,'1','a','b'],
            ['--check','9',s,'1'],['--check','00',s,'1'],['--check','0',s,'256'],
            ['--check','8',s,'0256'],['--check','8',s,'0'],['--check','0','0101','1'],
            ['--check','0',('0110' if profile=='native' else 'NNNN'),'1'],
            ['--clock-exercise','0',s,'1'],['--clock-guard','0',s,'1','forbidden'],
            ['--clock-guard','8',s,'256']]


def inventory():
    rows = []
    for p in ARCHIVES:
        schedules = ('NNNN',) if p=='native' else ('0110','1001','0000','1111')
        for c in range(9):
            for s in schedules:
                for g in ((1,256) if c==8 else (1,)):
                    for kind in ('plain','synthetic'):
                        rows.append((f'{p}-{c}-{s}-{g}-{kind}',p,c,s,g,kind,0,None))
        for c,g in ((0,1),(8,256)):
            for s in schedules:
                rows.append((f'{p}-{c}-{s}-{g}-abort',p,c,s,g,'abort',86,None))
        s = schedules[0]
        for fault in ('equal','reverse','negative','huge'):
            rows.append((f'{p}-fault-{fault}',p,8,s,256,'fault',1,fault))
        for i,args in enumerate(bad_args(p)):
            rows.append((f'{p}-bad-{i}',p,None,None,None,'bad',1,args))
        if p!='native':
            rows.append((f'{p}-unit',p,None,None,None,'unit',0,None))
    return rows


def replay(root, max_events):
    build,checks = read(root/'build/build.json'),read(root/'checks/checks.json')
    for state in (build,checks):
        equal([state['bead'],state['completed'],state['timed']],[BEAD,True,False])
    for name,digest in build['inputs'].items(): equal(sha(Path(name)),digest)
    for name,digest in build['artifacts'].items(): equal(sha(root/'build'/name),digest)
    for profile,digest in ARCHIVES.items():
        equal(sha(root/'build'/profile/'codec.a'),digest)
        # Every linked clock variant consumes the exact same driver object.
        for variant in ('plain','abort','synthetic'):
            matches = [c for c in build['commands'] if c[-1].endswith('/'+profile+'/'+variant)]
            require(len(matches)==1,'unique variant link')
            driver = [a for a in matches[0] if a.endswith('/'+profile+'/driver.o')]
            require(len(driver)==1,'same compiled driver object')
            symbols = (root/'build'/profile/(variant+'-undefined.txt')).read_text()
            equal('_ZNSt6chrono3_V212steady_clock3nowEv' in symbols,variant=='plain')
    rows = inventory()
    equal([r['label'] for r in checks['records']],[r[0] for r in rows])
    total,comparisons = 0,0
    for record,(label,p,c,s,g,kind,code,extra) in zip(checks['records'],rows):
        equal(record['returncode'],code)
        out,err = root/'checks'/(label+'.stdout'),root/'checks'/(label+'.stderr')
        equal(sha(out),record['stdout_sha256']); equal(sha(err),record['stderr_sha256'])
        lines = [parse(line) for line in out.read_text().splitlines()]
        by_schema = {r['schema']:r for r in lines}
        require(len(by_schema)==len(lines),'duplicate output schema')
        wanted = []
        binary = 'synthetic' if kind=='fault' else 'plain' if kind=='bad' else 'group-unit' if kind=='unit' else kind
        args = []
        parity = None
        if kind in ('plain','synthetic'):
            args = ['--clock-exercise' if kind=='synthetic' else '--exercise',str(c),s,str(g)]
            parity = root/'checks'/(label+'.parity')
            args.append(str(Path(checks['root'])/'checks'/parity.name))
            wanted.append(driver_expected(p,c,s,g,kind=='synthetic'))
            if kind=='synthetic':
                wanted += [clock_expected(g),witness(p,c,s,g,'exercise')]
        elif kind=='abort':
            args = ['--clock-guard',str(c),s,str(g)]
            wanted = [witness(p,c,s,g,'clock-guard')]
        elif kind=='fault':
            args = ['--clock-exercise','8',s,'256']
            wanted = [clock_expected(256,1),fault_witness(p,256)]
        elif kind=='bad': args = extra
        else: wanted = [dict(schema='paired-group-unit/v1',cases=37,timed=False)]
        equal(record['args'],[p,binary,*args])
        equal(record['fault'],extra if kind=='fault' else None)
        equal(by_schema,{r['schema']:r for r in wanted})
        if code==0: equal(err.read_text(),'')
        elif kind=='abort': equal(err.read_text(),'unexpected driver benchmark clock\n')
        elif kind=='fault':
            equal(err.read_text(), ('group duration exceeds exact binary64 integer range' if extra=='huge'
                  else 'nonpositive or reversed grouped clock interval')+'\n')
        else: require(bool(err.read_text()),'CLI refusal')
        if parity:
            equal(parity.stat().st_size,CELLS[c][1]*CELLS[c][2])
            equal(sha(parity),record['parity_sha256'])
            baseline = root/'checks'/f'native-{c}-NNNN-{g}-plain.parity'
            if p!='native' or kind!='plain':
                total += compare(parity,baseline); comparisons += 1
    equal(sorted(x.name for x in (root/'checks').glob('*.stdout')),
          sorted(r[0]+'.stdout' for r in rows))
    resources = {phase:scope((root/(phase+'.log')).read_text(),limit,max_events=max_events[phase])
                 for phase,limit in (('build',512*1024**2),('checks',256*1024**2))}
    return dict(bead=BEAD,positive=182,clock_aborts=18,clock_faults=12,cli_refusals=42,
                synthetic_spans=90*84,per_profile_group_unit_cases=37,
                parity_files=comparisons,parity_bytes=total,resources=resources,timed=False,default_enabled=False)


if __name__=='__main__':
    print(json.dumps(replay(Path(sys.argv[1]).resolve(),dict(build=int(sys.argv[2]),checks=int(sys.argv[3]))),sort_keys=True))

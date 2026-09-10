"""Retained-only verifier for same-header epoch boundary and synthetic units."""
from pathlib import Path
from verify_paired_epoch import (BEAD, ARCHIVES, PROGRESS, MARKS, WITNESS, CLOCK, public_oracle,
    parse, equal, require, sha, check_inventory, build_resources, split_scope, scope_command, records)

SOURCE_NAMES = ('test_paired_epoch_native.cpp', 'test_paired_epoch_clock.cpp')


def inventory():
    rows = []
    for p in ARCHIVES:
        for mode,code in (('geometry',0),('marks',0),('mark-order',87),('mark-duplicate',87),('mark-capacity',87)):
            rows.append(dict(label=p+'-'+mode, profile=p, kind='metadata', args=[] if mode=='geometry' else ['--'+mode],
                             code=code, mode=mode))
        for mode,code in (('capacity',0),('overflow',89),('unknown-fault',89),('epoch2-no-fault',0),
                          ('-1',89),('3',89),('00',89),('',89),('x',89)):
            rows.append(dict(label=p+'-clock-'+str(len(rows)), profile=p, kind='clock', args=[mode], code=code, mode=mode))
    return rows


def verify_record(row, output, error):
    actual = records(output); mode = row['mode']
    if row['kind'] == 'clock':
        count = 504 if mode in ('capacity','overflow','epoch2-no-fault') else 1
        wanted = {CLOCK:dict(schema=CLOCK, clock_calls=count, public_calls_at_clock=[0]*count, timed=False)}
        if row['code']==0:
            wanted['paired-epoch-clock-unit/v1'] = dict(schema='paired-epoch-clock-unit/v1',endpoints=504,timed=False)
            equal(error,'')
        else: equal(error, ('too many synthetic clocks' if mode=='overflow' else
                           'unknown synthetic clock fault' if mode=='unknown-fault' else 'invalid synthetic fault epoch')+'\n')
    else:
        w,_ = public_oracle('native',0,'NNNN',1,False,completed=0)
        count = 6 if mode in ('marks','mark-capacity') else 1 if mode=='mark-duplicate' else 0
        mark = {k:v for k,v in w.items() if k!='schema'}
        wanted = {WITNESS:w, MARKS:dict(schema=MARKS, marks=[dict(mark,endpoint=i) for i in range(count)],timed=False),
                  PROGRESS:dict(schema=PROGRESS,selection_count=0,snapshot_count=0,timed=False)}
        if row['code']==0:
            wanted['paired-epoch-unit/v1'] = dict(schema='paired-epoch-unit/v1',
                cases=(33 if row['profile']=='native' else 45) if mode=='geometry' else 1,timed=False)
            equal(error,'')
        else: equal(error,'epoch witness mark order/capacity\n')
    equal(actual,wanted)


def recipes(root, frontend):
    result=[]; unit=Path(root)/'units'
    for p in ARCHIVES:
        directory=Path(root)/'build'/p
        compiles=[c for c in frontend['commands'] if c[-1]==str(directory/'driver.o')]
        links=[c for c in frontend['commands'] if c[-1]==str(directory/'abort')]
        require(len(compiles)==len(links)==1,'unit recipe reference')
        command=compiles[0]; flags=command[command.index('c++')+1:command.index('-c')]
        for kind,source in zip(('metadata','clock'),SOURCE_NAMES):
            obj=str(unit/(p+'-'+kind+'.o')); binary=str(unit/(p+'-'+kind))
            result.append(['prlimit','--cpu=120:120','--','c++',*flags,'-I'+str(Path(root)/'build'),
                           '-c',str(unit/source),'-o',obj])
            if kind=='metadata':
                link=[obj if a==str(directory/'driver.o') else a for a in links[0]]
                link[-1]=binary
            else: link=['prlimit','--cpu=120:120','--','c++',*flags,obj,'-o',binary]
            result.append(link)
    return result


def replay(root):
    folder=root/'units'; build=parse((folder/'build.json').read_text()); checks=parse((folder/'checks.json').read_text())
    frontend=parse((root/'build/build.json').read_text())
    for state in (build,checks): equal([state['bead'],state['completed'],state['timed']],[BEAD,True,False])
    equal(build['commands'],recipes(frontend['root'],frontend))
    names=set(SOURCE_NAMES)|{'build_paired_epoch_units.py','verify_paired_epoch_units.py'}
    names|={p+'-'+k+s for p in ARCHIVES for k in ('metadata','clock') for s in ('','.o')}
    equal(sorted(build['artifacts']),sorted(names))
    for name,digest in build['artifacts'].items(): equal(sha(folder/name),digest)
    equal(sorted(build['source_sha256']),sorted(SOURCE_NAMES))
    for name,digest in build['source_sha256'].items(): equal(sha(folder/name),digest)
    rows=inventory(); equal([r['label'] for r in checks['records']],[r['label'] for r in rows])
    files=names|{'build.json','checks.json'}; peaks=[]
    for row,record in zip(rows,checks['records']):
        p,kind,label=row['profile'],row['kind'],row['label']
        command=scope_command(frontend['root'],p,'abort',row['args'])
        command[len(command)-len(row['args'])-1]=str(Path(frontend['root'])/'units'/(p+'-'+kind))
        equal(record['command'],command); equal(record['returncode'],row['code'])
        stdout,stderr=folder/(label+'.stdout'),folder/(label+'.stderr')
        files|={stdout.name,stderr.name}
        equal(sha(stdout),record['stdout_sha256']); equal(sha(stderr),record['stderr_sha256'])
        output,error,peak=split_scope(stdout.read_text(),stderr.read_text(),row['code'])
        equal(record['memory_peak'],peak); peaks.append(peak)
        verify_record(row,output,error)
    equal(sorted(p.name for p in folder.iterdir()),sorted(files))
    return dict(profiles=3,geometry_and_probe_cases=dict(native=33,release=45,sanitize=45),mark_successes=3,mark_refusals=9,
                clock_successes=6,clock_refusals=21,public_encode_calls=0,maximum_native_peak=max(peaks),
                build_peak=build_resources((root/'units-build.log').read_text()))

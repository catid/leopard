#!/usr/bin/env python3
"""Retained-only metadata qualification. No native execution or clock reads."""
import json
import hashlib
import os
from pathlib import Path
import re
import struct
import sys

from paired_metadata_overlay import BEAD, adapt
from verify_paired_r19932 import ARCHIVES, CELLS, equal, expected, witness, require
from verify_paired_timer_r19932 import clock_expected, fault_witness

META = 'paired-runtime-metadata/v1'
DRIVER = 'leopard-paired-timer-r19932/v1'
UINT64 = (1 << 64) - 1


def sha(path):
    # These are this lane's retained files, not global cache eviction. Drop
    # streamed read pages so the 1.58GB parity replay fits the unchanged cap.
    with path.open('rb') as stream:
        try: return hashlib.file_digest(stream,'sha256').hexdigest()
        finally: os.posix_fadvise(stream.fileno(),0,0,os.POSIX_FADV_DONTNEED)


def compare(left,right):
    total=0
    with left.open('rb') as a,right.open('rb') as b:
        try:
            while True:
                chunk=a.read(65536); require(chunk==b.read(65536),'full native parity differs')
                if not chunk: return total
                total+=len(chunk)
        finally:
            for stream in (a,b): os.posix_fadvise(stream.fileno(),0,0,os.POSIX_FADV_DONTNEED)


def parse(text):
    def unique(pairs):
        obj = {}
        for key,value in pairs:
            require(key not in obj,'duplicate JSON field')
            obj[key] = value
        return obj
    def constant(value):
        raise ValueError('nonfinite JSON number: '+value)
    return json.loads(text,object_pairs_hook=unique,parse_constant=constant)


def keys(obj, names):
    require(type(obj) is dict and set(obj) == set(names.split()),'metadata object fields')


def integer(value, minimum=0, maximum=UINT64):
    require(type(value) is int and minimum <= value <= maximum,'metadata integer/range')
    return value


def end(address, size):
    integer(address); integer(size)
    require(size <= UINT64-address,'metadata address overflow')
    return address+size


def elf(path, native):
    """Read actual ELF64 program/symbol tables, independently of saved nm text."""
    with path.open('rb') as stream:
        length = path.stat().st_size
        def read(offset, size):
            require(0 <= offset <= length and 0 <= size <= min(8*1024**2,length-offset),'ELF range')
            stream.seek(offset); data = stream.read(size)
            require(len(data)==size,'short ELF read')
            return data
        header = struct.unpack('<16sHHIQQQIHHHHHH',read(0,64))
        require(header[0][:7] == b'\x7fELF\x02\x01\x01' and header[1] in (2,3) and header[2]==62,'ELF64 x86-64')
        phoff,shoff,phsize,phnum,shsize,shnum = header[5],header[6],header[9],header[10],header[11],header[12]
        require(phsize==56 and shsize==64 and 0 < phnum <= 64 and 0 < shnum <= 256,'ELF table geometry')
        segments = []
        for i in range(phnum):
            p = struct.unpack('<IIQQQQQQ',read(phoff+i*phsize,phsize))
            if p[0] == 1:
                segments.append(dict(address=p[3],bytes=p[6],flags=p[1]))
        require(0 < len(segments) <= 16,'ELF load count')
        sections = [struct.unpack('<IIQQQQIIQQ',read(shoff+i*shsize,shsize)) for i in range(shnum)]
        tables = [s for s in sections if s[1]==2]
        require(len(tables)==1,'ELF static symbol table')
        table = tables[0]; require(table[9]==24 and table[5]%24==0 and table[6]<shnum,'ELF symbols')
        strings = sections[table[6]]; require(strings[1]==3,'ELF string table')
        names = read(strings[4],strings[5]); data = read(table[4],table[5])
        functions,undefined = {},[]
        desired = {'main':'main', 'leo_encode' if native else 'leo2_encode':'encode'}
        if not native: desired['leo2_encode_batch']='batch'
        for i in range(0,len(data),24):
            name,info,_,index,value,size = struct.unpack_from('<IBBHQQ',data,i)
            require(name<len(names),'ELF symbol name')
            stop = names.find(b'\0',name); require(stop>=name,'ELF symbol terminator')
            name = names[name:stop].decode('ascii')
            if index==0: undefined.append(name)
            label = desired.get(name)
            if 'paired_metadata' in name and name.endswith('6AnchorEv'): label='metadata_anchor'
            if label:
                require(label not in functions and info&15==2 and index!=0 and size>0,'defined function symbol')
                functions[label]=dict(value=value,size=size)
        require(set(functions)==set(desired.values())|{'metadata_anchor'},'function inventory')
        require('_ZNSt6chrono3_V212steady_clock3nowEv' not in undefined,'real benchmark clock import')
        return dict(type=header[1],segments=segments,functions=functions)


def validate_metadata(meta, profile, cell, schedule, exercise, image):
    base = expected(profile,cell,schedule,1,exercise)
    keys(meta,'schema bead timed observation pointer_bytes preflight_gfni_counts selections native_route_label snapshots')
    native = profile=='native'; k,r,size = CELLS[cell]
    equal([meta[n] for n in ('schema','bead','timed','observation','pointer_bytes','preflight_gfni_counts','native_route_label')],
          [META,BEAD,False,'new_frontend_endpoints',8,base['probes'],
           'original_native_compiler_policy' if native else 'not_native'])
    want = []
    for i in range(104 if exercise else 4):
        slot = i%4
        want.append(dict(index=i,phase='preflight' if i<4 else 'exercise',pass_=-1 if i<4 else (i-4)//4,
            slot=slot,requested_backend=-1 if native else 3 if cell==6 else 0,
            context_backend=-1 if native else 3,candidate_state=-1 if native else int(schedule[slot]),
            operation_gfni=-1 if native else int(2<=cell<=4 or (cell<2 and schedule[slot]=='1'))))
        want[-1]['pass']=want[-1].pop('pass_')
    equal(meta['selections'],want)
    require(type(meta['snapshots']) is list and len(meta['snapshots'])==2,'two endpoint snapshots')
    for endpoint,snapshot in enumerate(meta['snapshots']):
        keys(snapshot,'endpoint allocations spans inputs outputs functions dladdr_image_base load_bias page_bytes segments')
        equal(snapshot['endpoint'],endpoint)
        keys(snapshot['allocations'],'source reference scratch output')
        sizes = dict(source=k*size,reference=r*size,scratch=base['scratch_bytes'],output=0 if native else r*size)
        ranges=[]
        for name,allocation in snapshot['allocations'].items():
            keys(allocation,'raw data bytes'); equal(allocation['bytes'],sizes[name])
            raw,data = integer(allocation['raw']),integer(allocation['data'])
            if native and name=='output': equal(allocation,dict(raw=0,data=0,bytes=0)); continue
            require(raw>0 and raw%64==0 and data==end(raw,64),'allocation alignment/guard')
            ranges.append((raw,end(end(data,sizes[name]),64)))
        keys(snapshot['spans'],'parity input_array output_array')
        for span in snapshot['spans'].values():
            keys(span,'address bytes'); integer(span['address'],1); end(span['address'],span['bytes'])
        parity = snapshot['spans']['parity']
        parity_data = snapshot['allocations']['scratch' if native else 'output']['data']
        equal(parity,dict(address=parity_data,bytes=r*size))
        output_count = base['scratch_bytes']//size if native else r
        require(output_count<=1024 and k<=4096,'pointer capacity')
        for name,count in (('input_array',k),('output_array',output_count)):
            span = snapshot['spans'][name]
            require(span['address']%8==0,'pointer array alignment'); equal(span['bytes'],count*8)
            ranges.append((span['address'],end(span['address'],span['bytes'])))
        for i,a in enumerate(ranges):
            for b in ranges[i+1:]: require(a[1]<=b[0] or b[1]<=a[0],'allocation/array overlap')
        equal(snapshot['inputs'],[snapshot['allocations']['source']['data']+i*size for i in range(k)])
        equal(snapshot['outputs'],[parity_data+i*size for i in range(output_count)])
        page=integer(snapshot['page_bytes'],1); equal(page,4096)
        bias=integer(snapshot['load_bias']); require(bias%page==0,'load bias alignment')
        if image['type']==2: equal(bias,0)
        else: require(bias>0,'PIE placement')
        loads=[dict(address=end(bias,s['address']),bytes=s['bytes'],flags=s['flags']) for s in image['segments']]
        equal(snapshot['segments'],loads)
        equal(snapshot['dladdr_image_base'],bias+min(s['address']//page*page for s in image['segments']))
        functions={name:end(bias,symbol['value']) for name,symbol in image['functions'].items()}
        if native: functions['batch']=0
        equal(snapshot['functions'],functions)
        for name,symbol in image['functions'].items():
            address=functions[name]
            require(any(s['flags']&1 and s['address']<=address and end(address,symbol['size'])<=end(s['address'],s['bytes'])
                        for s in loads),'function outside executable mapping')
    equal({k:v for k,v in meta['snapshots'][0].items() if k!='endpoint'},
          {k:v for k,v in meta['snapshots'][1].items() if k!='endpoint'})


def inventory():
    rows=[]
    def add(p,c,s,g,mode,binary,code=0,fault=None,extra=None):
        label=f'{p}-{c}-{s}-{g}-{mode.lstrip("-")}' + (('-'+fault) if fault else '')
        if extra is not None: label=f'{p}-bad-{len(rows)}'
        rows.append(dict(label=label,profile=p,binary=binary,args=extra if extra is not None else [mode,str(c),s,str(g)],
                         code=code,fault=fault,parity=code==0))
    for p in ARCHIVES:
        schedules=('NNNN',) if p=='native' else ('0110','1001','0000','1111')
        for c in range(9):
            for s in schedules:
                for g in ((1,256) if c==8 else (1,)):
                    for mode,binary in (('--check','abort'),('--exercise','abort'),('--clock-exercise','synthetic')):
                        add(p,c,s,g,mode,binary)
        for c,g in ((0,1),(8,256)):
            for s in schedules: add(p,c,s,g,'--clock-guard','abort',86)
        for fault in ('equal','reverse','negative','huge'):
            add(p,8,schedules[0],256,'--clock-exercise','synthetic',1,fault)
        s=schedules[0]
        bad=[[],['--check'],['--measure','0',s,'1'],['--timing','0',s,'1'],
             ['--check','9',s,'1'],['--check','00',s,'1'],['--check','0',s,'256'],
             ['--check','8',s,'0256'],['--check','8',s,'0'],['--check','0','0101','1'],
             ['--check','0',('0110' if p=='native' else 'NNNN'),'1'],
             ['--clock-exercise','0',s,'1'],['--clock-guard','0',s,'1','forbidden']]
        for args in bad: add(p,None,None,None,'bad','abort',1,extra=args)
        add(p,None,None,None,'bad','synthetic',1,extra=['--clock-guard','0',s,'1'])
    require(len({r['label'] for r in rows})==len(rows),'unique inventory')
    return rows


def split_scope(output,error,code):
    require(output.count('TOWER_CHILD_EXIT=')==1,'scope footer')
    native,footer=output.split('TOWER_CHILD_EXIT=')
    lines=footer.splitlines()
    equal(lines[:2],[str(code),'memory.peak']); peak=integer(int(lines[2]),1,256*1024**2)
    equal(lines[3:],['memory.max','268435456','memory.events','low 0','high 0','max 0','oom 0',
                    'oom_kill 0','oom_group_kill 0','memory.swap.current','0','memory.swap.max','0'])
    header,stderr=error.split('\n',1)
    require(re.fullmatch(r'Running as unit: run-[\w-]+\.scope; invocation ID: [0-9a-f]{32}',header) is not None,'scope header')
    return native,stderr,peak


def scope_command(root,profile,kind,args):
    root=Path(root)
    return ['systemd-run','--user','--scope','--expand-environment=no','-p','MemoryMax=256M','-p','MemorySwapMax=0',
        'bash',str(root/'build/tower_public_scope.sh'),'flock','-n','/tmp/leopard-gf8-authoritative.lock',
        'timeout','--signal=TERM','--kill-after=5','120','prlimit','--cpu=60:60','--core=0:0','--',
        str(root/'build'/profile/kind),*args]


def check_inventory(directory,pins,exclude=()):
    require(type(pins) is dict and bool(pins),'nonempty digest inventory')
    actual=[]
    for path in sorted(directory.rglob('*')):
        require(not path.is_symlink(),'no retained artifact links')
        if path.is_file() and str(path.relative_to(directory)) not in exclude:
            actual.append(str(path.relative_to(directory)))
    equal(sorted(pins),actual)
    for name,digest in pins.items():
        require(type(digest) is str and re.fullmatch('[0-9a-f]{64}',digest) is not None,'digest syntax')
        equal(sha(directory/name),digest)


def build_resources(text):
    marker='TOWER_CHILD_EXIT=0\nmemory.peak\n'
    require(text.count(marker)==1 and '\tExit status: 0\n' in text,'successful build scope')
    rest=text.split(marker)[1]; value,rest=rest.split('\n',1)
    peak=integer(int(value),1,512*1024**2)
    require(rest.startswith('memory.max\n536870912\nmemory.events\nlow 0\nhigh 0\nmax 0\noom 0\n'
        'oom_kill 0\noom_group_kill 0\nmemory.swap.current\n0\nmemory.swap.max\n0\n'),'build resources')
    return peak


def check_build_inputs(root,inputs):
    require(type(inputs) is dict and bool(inputs),'nonempty source inputs')
    # Reconstruct from the fixed builder's source contract, not map coverage alone.
    source_names={'build_paired_metadata.py','paired_metadata_overlay.py','PairedRuntimeMetadata.h',
                  'paired_timer_r19932.cpp','paired_timer_witness.cpp','paired_public_witness.cpp',
                  'paired_timer_clock.cpp','PairedGroupTiming.h','tower_public_scope.sh'}
    seen=set()
    for name,digest in inputs.items():
        path=Path(name)
        if path.parent==Path('/home/catid/leopard/experiments/leopard2/gf16_high_encode'):
            require(path.name in source_names,'known frontend input')
            targets=[root/'build'/path.name]; identity=('source',path.name)
        elif path.name=='clock_guard.cpp':
            targets=[root/'build/clock_guard.cpp']; identity=('guard',path.name)
        elif path.suffix=='.a':
            profiles=[p for p,d in ARCHIVES.items() if d==digest]
            require(len(profiles)==1,'known unchanged archive input')
            targets=[root/'build'/profiles[0]/'codec.a']; identity=('archive',profiles[0])
        elif path.suffix=='.h':
            native='/pure-checks-v2/source/' in name
            profiles=['native'] if native else ['release','sanitize']
            targets=[root/'build'/p/'include'/path.name for p in profiles]
            identity=('native-header' if native else 'l2-header',path.name)
        else: raise ValueError('unknown build input')
        require(identity not in seen,'duplicate logical build input'); seen.add(identity)
        for target in targets: equal(sha(target),digest)
    required={('source',n) for n in source_names}|{('guard','clock_guard.cpp')}|{('archive',p) for p in ARCHIVES}
    for profile,label in (('native','native-header'),('release','l2-header')):
        required|={(label,p.name) for p in (root/'build'/profile/'include').glob('*.h')}
    equal(sorted(seen),sorted(required))


def units(root,original_root):
    folder=root/'units'; build=parse((folder/'build.json').read_text()); checks=parse((folder/'checks.json').read_text())
    for state in (build,checks): equal([state['bead'],state['completed'],state['timed']],[BEAD,True,False])
    names={'test_paired_runtime_metadata.cpp','build_paired_metadata_units.py'}|{
        p+suffix for p in ARCHIVES for suffix in ('','.o')}
    equal(sorted(build['artifacts']),sorted(names))
    for name,digest in build['artifacts'].items(): equal(sha(folder/name),digest)
    equal(sha(folder/'test_paired_runtime_metadata.cpp'),build['source_sha256'])
    original=parse((root/'build/build.json').read_text())
    recipes=[]
    for profile in ARCHIVES:
        directory=Path(original_root)/'build'/profile
        compiles=[c for c in original['commands'] if c[-1]==str(directory/'driver.o')]
        links=[c for c in original['commands'] if c[-1]==str(directory/'abort')]
        require(len(compiles)==len(links)==1,'original unit recipe reference')
        command=compiles[0]; flags=command[command.index('c++')+1:command.index('-c')]
        unit=Path(original_root)/'units'
        recipes.append(['prlimit','--cpu=120:120','--','c++',*flags,'-I'+str(Path(original_root)/'build'),
                        '-c',str(unit/'test_paired_runtime_metadata.cpp'),'-o',str(unit/(profile+'.o'))])
        link=[str(unit/(profile+'.o')) if a==str(directory/'driver.o') else a for a in links[0]]
        link[-1]=str(unit/profile); recipes.append(link)
    equal(build['commands'],recipes)
    equal([r['profile'] for r in checks['records']],list(ARCHIVES))
    peaks=[]
    for record in checks['records']:
        p=record['profile']; equal(record['returncode'],0)
        command=scope_command(original_root,p,'abort',[]); command[-1]=str(Path(original_root)/'units'/p)
        equal(record['command'],command)
        out,err=folder/(p+'.stdout'),folder/(p+'.stderr')
        equal(sha(out),record['stdout_sha256']); equal(sha(err),record['stderr_sha256'])
        native,error,peak=split_scope(out.read_text(),err.read_text(),0); equal(error,''); peaks.append(peak)
        equal(record['memory_peak'],peak)
        equal([parse(line) for line in native.splitlines()],
              [dict(schema='paired-metadata-unit/v1',cases=22,timed=False),
               dict(schema='leopard-paired-witness/v1',calls=0,states=[0,0,0],apis=[0,0,0],order_hash='cbf29ce484222325')])
        links=[c for c in build['commands'] if c[-1]==str(Path(original_root)/'units'/p)]
        require(len(links)==1,'unique unit link')
        require(str(Path(original_root)/'build'/p/'codec.a') in links[0],'same unit archive')
    return dict(profiles=3,cases_per_profile=22,public_encode_calls=0,maximum_native_peak=max(peaks),
                build_peak=build_resources((root/'units-build.log').read_text()))


def verify_record(row,output,error,image):
    lines=[parse(line) for line in output.splitlines()]
    actual={r['schema']:r for r in lines}; require(len(actual)==len(lines),'duplicate schemas')
    p,kind,args=row['profile'],row['binary'],row['args']
    if '-bad-' in row['label']:
        equal(actual,{'leopard-paired-witness/v1':dict(schema='leopard-paired-witness/v1',calls=0,
              states=[0,0,0],apis=[0,0,0],order_hash='cbf29ce484222325'),
              **({'paired-synthetic-clock/v1':clock_expected(1,0)} if kind=='synthetic' else {})})
        require(bool(error),'CLI refusal diagnostic')
        if args and args[0]=='--measure': equal(error,'metadata refuses real timing\n')
        return
    mode,c,s,g=args; c,g=int(c),int(g)
    wanted=[]
    if row['code']==0:
        equal(error,''); exercise=mode!='--check'
        base=expected(p,c,s,g,exercise)
        base.update(schema=DRIVER,clock_source=kind,
            samples=[[257+17*i,(257+17*i)/g] for i in range(84)] if kind=='synthetic' else [])
        wanted=[base,witness(p,c,s,g,'exercise' if exercise else 'check')]
        validate_metadata(actual.pop(META),p,c,s,exercise,image)
        if kind=='synthetic': wanted.append(clock_expected(g))
    elif row['code']==86:
        equal(error,'unexpected driver benchmark clock\n'); wanted=[witness(p,c,s,g,'clock-guard')]
    else:
        equal(error,('group duration exceeds exact binary64 integer range' if row['fault']=='huge'
                     else 'nonpositive or reversed grouped clock interval')+'\n')
        wanted=[fault_witness(p,g),clock_expected(g,1)]
    equal(actual,{r['schema']:r for r in wanted})


def replay(root):
    build=parse((root/'build/build.json').read_text()); checks=parse((root/'checks/checks.json').read_text())
    for state in (build,checks): equal([state['bead'],state['completed'],state['timed']],[BEAD,True,False])
    check_inventory(root/'build',build['artifacts'],exclude=('build.json',))
    check_inventory(root/'tools',checks['tools'])
    require({'run_paired_metadata_checks.py','verify_paired_metadata.py','paired_metadata_overlay.py'} <= set(checks['tools']),
            'required qualification tools')
    check_build_inputs(root,build['inputs'])
    build_peak=build_resources((root/'build.log').read_text())
    equal((root/'build/paired_metadata.cpp').read_text(),adapt((root/'build/paired_timer_r19932.cpp').read_text()))
    images={}
    for p,digest in ARCHIVES.items():
        equal(sha(root/'build'/p/'codec.a'),digest)
        for kind in ('abort','synthetic'):
            matches=[c for c in build['commands'] if c[-1].endswith('/'+p+'/'+kind)]
            require(len(matches)==1,'unique link')
            command=matches[0]
            equal(sum(a.endswith('/'+p+'/driver.o') for a in command),1)
            require('-ldl' in command and '-Wl,--wrap=_ZNSt6chrono3_V212steady_clock3nowEv' in command,'metadata linkage')
            images[p,kind]=elf(root/'build'/p/kind,p=='native')
    unit_result=units(root,checks['root'])
    rows=inventory(); equal([r['label'] for r in checks['records']],[r['label'] for r in rows])
    totals=dict(positive=0,clock_aborts=0,clock_faults=0,cli_refusals=0,parity_comparisons=0,parity_bytes=0,
                selections=0,snapshots=0,public_encode_calls=0,synthetic_spans=0)
    peaks=[]
    for record,row in zip(checks['records'],rows):
        label,p,kind=row['label'],row['profile'],row['binary']
        equal(record['returncode'],row['code']); equal(record['fault'],row['fault'])
        args=list(row['args']); parity=root/'checks'/(label+'.parity') if row['parity'] else None
        if parity: args.append(str(Path(checks['root'])/'checks'/parity.name))
        equal(record['args'],[p,kind,*args])
        equal(record['scope_command'],scope_command(checks['root'],p,kind,args))
        out,err=root/'checks'/(label+'.stdout'),root/'checks'/(label+'.stderr')
        equal(sha(out),record['stdout_sha256']); equal(sha(err),record['stderr_sha256'])
        native,error,peak=split_scope(out.read_text(),err.read_text(),row['code']); peaks.append(peak)
        equal(record['memory_peak'],peak)
        verify_record(row,native,error,images[p,kind])
        if row['code']==0:
            totals['positive']+=1; totals['snapshots']+=2
            totals['selections']+=4 if row['args'][0]=='--check' else 104
            equal(sha(parity),record['parity_sha256']); c=int(row['args'][1])
            totals['public_encode_calls']+=expected(p,c,row['args'][2],int(row['args'][3]),row['args'][0]!='--check')['encode_calls']
            totals['synthetic_spans']+=84 if kind=='synthetic' else 0
            equal(parity.stat().st_size,CELLS[c][1]*CELLS[c][2])
            baseline=root/'checks'/f'native-{c}-NNNN-1-check.parity'
            if parity!=baseline:
                totals['parity_bytes']+=compare(parity,baseline); totals['parity_comparisons']+=1
        else: totals['cli_refusals' if '-bad-' in label else 'clock_aborts' if row['code']==86 else 'clock_faults']+=1
    equal(sorted(p.name for p in (root/'checks').glob('*.stdout')),sorted(r['label']+'.stdout' for r in rows))
    final=parse((root/'final-tools.json').read_text())
    equal([final['bead'],final['timed']],[BEAD,False]); check_inventory(root/'final_tools',final['files'])
    require({'verify_paired_metadata.py','test_paired_metadata.py','paired_metadata_overlay.py'} <= set(final['files']),
            'final replay tools')
    return dict(bead=BEAD,timed=False,default_enabled=False,totals=totals,maximum_native_peak=max(peaks),build_peak=build_peak,
                boundary_units=unit_result,
                all_native_memory_events_zero=True,swap_bytes=0,endpoint_stability_only=True,
                historical_shift_cause_established=False)


if __name__=='__main__': print(json.dumps(replay(Path(sys.argv[1]).resolve(strict=True)),sort_keys=True))

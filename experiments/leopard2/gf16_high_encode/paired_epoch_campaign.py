"""Pure campaign identities and validation; no acquisition or codec launch."""
import json
import os
from pathlib import Path
import re
import sys

import verify_paired_epoch_timing as qualified
from verify_paired_metadata import sha,parse,equal,require

BEAD=qualified.BEAD
ROOT=Path('/home/catid/leopard/.research/leopard-79h/paired-epoch-campaign.fSHScp')
QUALIFIED=Path('/home/catid/leopard/.research/leopard-79h/paired-epoch-timing-qualified.8oRYqG')
PLAN='paired_epoch_campaign_plan.json'
MANIFEST_PIN='b2970bbcc95e0bc9eaf63a96b53bd20039c809187a5930c7426d593871ef4abd'
BUILD_PIN='532104bfea50f997fd2bb258a0d1940c71903b874d10020dbc94850dc3ad360b'
REPLAY_PIN='03e8feb5755acf8ef2cc844d9ca003991ff25f33d5f8c3c79d9502cb29d4ba14'
FINAL_PIN='f310055a8bdad13604a370b01602133167b36896681cd9db1e16ca16854184a0'
CONDITION_PIN='88caa4c2df54c11dae472d5768eae2e0696474c3a7686619d1615ec35a65f944'
NEW_SOURCES={'paired_epoch_campaign.py','freeze_paired_epoch_campaign.py','run_paired_epoch_campaign.py',
             'replay_paired_epoch_campaign.py','test_paired_epoch_campaign.py'}
SOURCES=qualified.FINAL_PYTHON | NEW_SOURCES
SOURCE_ROOTS=qualified.FINAL_ROOTS+tuple(sorted(NEW_SOURCES))
ASSETS=set(qualified.FINAL_ASSETS)
PROOFS={'qualification.manifest':'SHA256SUMS','qualification.build.json':'build/build.json',
        'qualification.replay.json':'replay-normal.json','qualification.tools.json':'final-tools.json'}
EXTRA={PLAN,'paired_epoch_campaign_method.md','paired_epoch_condition.sh','paired_epoch_scope.sh','tower_public_scope.sh'}
FILES=SOURCES | ASSETS | set(PROOFS) | EXTRA | {'native','current'}
ENV=dict(PATH='/usr/bin:/bin',LANG='C',LC_ALL='C',OMP_NUM_THREADS='1',OMP_DYNAMIC='FALSE',OMP_THREAD_LIMIT='1')
HOST=dict(hostname='work',kernel='6.8.0-139-generic',vendor_id='AuthenticAMD',
          **{'cpu family':'26','model':'8','model name':'AMD Ryzen Threadripper 9980X 64-Cores'})


def plan(ready=False):
    return dict(schema='leopard-epoch-campaign-plan/v1',bead=BEAD,ready_to_time=ready,
        host=HOST,cpu=26,sibling=90,controller_cpu=0,passive_seconds=10,attempt_budget=1,
        attempt_root=str(ROOT/'attempt1'),qualification_path=str(ROOT/'qualification-v2'),timing_path=str(ROOT/'frozen'),
        qualification_output=str(ROOT/'qualification-checks-v2'),qualification_resource=str(ROOT/'qualification-resource-v2.log'),
        epochs=3,rounds=3,sample_passes=21,warmup_passes=4,spans_per_process=252,groups=[1]*8+[256],
        timed_processes=318,preflights=27,controls_per_epoch=100,homogeneous_processes=264,
        minimum_gain=1.05,native_minimum_gain=1.05,equivalence_bound=1.02,every_target_round_positive=True,
        diagnostic_only=True,production_promotion=False,default_enabled=False,confidence_intervals=False,
        epoch_aggregation='none_keep_all_three',condition='slipgate-disabled-20260909',
        codec_commit='45e2effd869859c9b3aa48190eff6f4738817c61',native_commit='6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198',
        frontend_commit='b64efb7957e1d99bea70435dd509f3bf4e80a48f',qualification_manifest_sha256=MANIFEST_PIN,
        cells=[[1000,199,32768,'encode','target'],[1000,199,32768,'one_item_batch','target'],
               [1000,200,32768,'encode','unchanged_neighbor'],[1000,199,65536,'encode','unchanged_neighbor'],
               [1000,200,65536,'encode','unchanged_neighbor'],[1000,198,32768,'encode','unchanged_neighbor'],
               [1000,199,32768,'explicit_avx2_encode','unchanged_neighbor'],[4096,512,4096,'encode','unchanged_neighbor'],
               [17,7,64,'gf8_AUTO_encode','unchanged_neighbor']])


def validate_plan(value,timing=False):
    require(type(value) is dict and type(value.get('ready_to_time')) is bool,'readiness boolean')
    equal(value,plan(value['ready_to_time']))
    if timing: require(value['ready_to_time'] is True,'collector not qualified/preregistered for timing')


def read(path):
    return parse(bounded_text(path,4*1024**2))


def bounded_text(path,limit=1024**2):
    require(path.is_file() and not path.is_symlink() and path.stat().st_size<limit,'bounded regular file')
    with path.open('rb') as stream:
        data=stream.read(limit)
    require(len(data)<limit,'bounded read')
    return data.decode('utf-8')


def condition_env():
    return dict(ENV,XDG_RUNTIME_DIR=f'/run/user/{os.getuid()}',
                DBUS_SESSION_BUS_ADDRESS=f'unix:path=/run/user/{os.getuid()}/bus')


def manifest(path):
    require(path.is_file() and not path.is_symlink() and path.stat().st_size<1024**2,'bounded manifest')
    equal(sha(path),MANIFEST_PIN); result={}
    for line in path.read_text().splitlines():
        digest,name=line.split('  ',1)
        require(re.fullmatch('[0-9a-f]{64}',digest) is not None and name not in result,'manifest identity')
        require(not Path(name).is_absolute() and '..' not in Path(name).parts and str(Path(name))==name,'manifest path')
        result[name]=digest
    return result


def executing_sources(pins):
    for module_name,module in tuple(sys.modules.items()):
        source=getattr(module,'__file__',None); name=Path(source).name if source else None
        expected=module_name+'.py'
        if expected not in SOURCES and name not in SOURCES: continue
        require(source is not None and name in SOURCES,'executing source: '+module_name)
        if expected in SOURCES: equal(name,expected)
        require(sha(Path(source))==pins['files'][name],'executing dependency hash: '+name)


def inputs(bundle, *, executing=True):
    require(bundle.is_dir() and not bundle.is_symlink() and not bundle.stat().st_mode&0o222,'immutable bundle')
    pins=read(bundle/'pins.json'); equal(sorted(pins),['files','schema'])
    equal(pins['schema'],'leopard-epoch-campaign-pins/v1'); equal(sorted(pins['files']),sorted(FILES))
    equal(sorted(p.name for p in bundle.iterdir()),sorted(FILES|{'pins.json'}))
    for name in FILES|{'pins.json'}:
        path=bundle/name
        require(path.is_file() and not path.is_symlink() and not path.stat().st_mode&0o222,'immutable file: '+name)
        if name in ('native','current'): equal(path.stat().st_mode&0o777,0o555)
        if name!='pins.json': equal(sha(path),pins['files'][name])
    value=read(bundle/PLAN); validate_plan(value)
    history=manifest(bundle/'qualification.manifest')
    for name,pin in (('qualification.build.json',BUILD_PIN),('qualification.replay.json',REPLAY_PIN),
                     ('qualification.tools.json',FINAL_PIN),('paired_epoch_condition.sh',CONDITION_PIN)):
        equal(pins['files'][name],pin)
    for name,original in PROOFS.items():
        if original!='SHA256SUMS': equal(pins['files'][name],history[original])
    old_tools=read(bundle/'qualification.tools.json')['files']
    for name in qualified.FINAL_PYTHON | ASSETS:
        equal(pins['files'][name],old_tools[name]); equal(pins['files'][name],history['final_tools/'+name])
    equal(pins['files']['tower_public_scope.sh'],history['build/tower_public_scope.sh'])
    build=read(bundle/'qualification.build.json'); equal([build['completed'],build['real_clocks_read']],[True,False])
    proof=read(bundle/'qualification.replay.json')
    equal([proof['bead'],proof['real_clocks_read'],proof['final_tools_sha256']],[BEAD,False,FINAL_PIN])
    images={}
    for name,profile in (('native','native'),('current','release')):
        original='measurement/'+profile+'/steady'
        equal(pins['files'][name],history['build/'+original]); equal(pins['files'][name],build['artifacts'][original])
        images[profile]=qualified.elf(bundle/name,profile=='native','steady')
        equal(images[profile],build['images'][original])
    equal(sorted(qualified.epoch.tool_closure(bundle,SOURCE_ROOTS)),sorted(SOURCES))
    if executing: executing_sources(pins)
    return value,pins,images


def preflights():
    return [dict(cell=c,order=s) for c in range(9) for s in ('NNNN','0000','1111')]


def name(item):
    if 'comparison' not in item: return f"check-{item['cell']}-{item['order']}"
    return f"cell-{item['cell']}-round-{item['round']}-{item['comparison']}-slot-{item['slot']}"


def command(bundle,cell,order,measured):
    require(type(cell) is int and 0<=cell<9 and order in ('NNNN','0110','1001','0000','1111'),'command workload')
    require(type(measured) is bool,'command mode')
    return ['/usr/bin/taskset','-c','26','/usr/bin/prlimit','--cpu=30:30','--core=0:0',
            '--fsize=1048576:1048576','--',str(bundle/('native' if order=='NNNN' else 'current')),
            '--measure' if measured else '--check',str(cell),order,'256' if cell==8 else '1']


def check_record(output,cell,order,image):
    profile='native' if order=='NNNN' else 'release'
    row=dict(profile=profile,variant='steady',args=['--check',str(cell),order,'256' if cell==8 else '1'],
             code=0,reason=None,fault=None,fault_epoch=None)
    qualified.verify_record(row,output,'',image)


def condition_lines():
    lines=[]
    for unit in ('slipgate-catstream.service','slipgate-obs.service','slipgate-headless-x.service'):
        lines+=['MainPID=0','Id='+unit,'ActiveState=inactive','UnitFileState=disabled']
    for container in ('3fad3dfb247b3d4f3991e0a6801fba0805aae0901bdc89018841d4707a2d888c',
                      '002ba0a05bf860b12f2d799b49da1976aed663cb1494b2f8236751e18dc50182'):
        lines.append(container+' exited|false|no')
    return lines


def controller_command(bundle,output,qualification,commit=None):
    script='freeze_paired_epoch_campaign.py' if qualification else 'run_paired_epoch_campaign.py'
    args=['check',str(bundle),str(output)] if qualification else [str(bundle),commit]
    return ['/usr/bin/timeout','--signal=TERM','--kill-after=5','1200','/usr/bin/prlimit','--cpu=600:600',
            '--core=0:0','--','/usr/bin/python3','-B',str(bundle/script),*args]


def disjoint_attempt(path):
    path=path.resolve(); attempt=(ROOT/'attempt1').resolve()
    require(path!=attempt and path not in attempt.parents and attempt not in path.parents,'reserved timing attempt overlap')


def resource(path,scope,command,wrapper):
    lines=bounded_text(path).splitlines()
    equal([line for line in lines if line.startswith('EPOCH_SCOPE=')],['EPOCH_SCOPE='+scope])
    equal([line for line in lines if line.startswith('EPOCH_WRAPPER=')],['EPOCH_WRAPPER='+str(wrapper)])
    equal(lines.count('EPOCH_COMMAND_BEGIN'),1); equal(lines.count('EPOCH_COMMAND_END'),1)
    equal(lines[lines.index('EPOCH_COMMAND_BEGIN')+1:lines.index('EPOCH_COMMAND_END')],command)
    equal([line for line in lines if line.startswith('TOWER_CHILD_EXIT=')],['TOWER_CHILD_EXIT=0'])
    equal([line for line in lines if line.startswith('\tExit status:')],['\tExit status: 0'])
    for field in ('memory.peak','memory.max','memory.events','memory.swap.current','memory.swap.max'):
        equal(lines.count(field),1)
    equal(lines.count('memory.peak'),1); i=lines.index('memory.peak'); peak=int(lines[i+1])
    require(0<peak<256*1024**2,'memory peak')
    equal(lines[i+2:i+15],['memory.max','268435456','memory.events','low 0','high 0','max 0',
         'oom 0','oom_kill 0','oom_group_kill 0','memory.swap.current','0','memory.swap.max','0'])
    return peak

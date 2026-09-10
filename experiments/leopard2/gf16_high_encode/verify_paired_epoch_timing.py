#!/usr/bin/env python3
"""Retained-only oracle for new frontend qualification; never reads clocks."""
import json
from pathlib import Path
import sys

import paired_epoch_timing_overlay as overlay
from paired_epoch_timing_elf import elf, VARIANTS, CLOCK
import verify_paired_epoch as epoch
from verify_paired_metadata import (sha, parse, equal, require, check_inventory, build_resources,
                                    split_scope, scope_command as old_scope)
from verify_paired_r19932 import ARCHIVES, CELLS, expected

BEAD = overlay.BEAD
REFERENCE_PIN = '7de095968e41289605e9164d29f585de934b822f570cca72d1bb4c3c44d5ca72'
NATIVE_CHECKS_PIN = '64a8218dcab9dff5e992c6b42c9f2e11b6e107c1d04a6b5474f160792cc5c1f6'
ASSETS = ('paired_metadata.cpp','PairedRuntimeMetadata.h','PairedGroupTiming.h',
          'paired_timer_clock.cpp','paired_timer_witness.cpp','paired_public_witness.cpp',
          'clock_guard.cpp','tower_public_scope.sh')
REUSED = ('codec.a','witness.o','abort-clock.o','synthetic-clock.o')
BUILD_TOOLS = epoch.BUILD_TOOLS | {'build_paired_epoch_timing.py','paired_epoch_timing_overlay.py',
                                  'paired_epoch_timing_elf.py'}
RUN_TOOLS = BUILD_TOOLS | {'run_paired_epoch_timing_checks.py','verify_paired_epoch_timing.py',
                         'verify_paired_epoch.py','verify_paired_epoch_units.py'}
FINAL_ROOTS = ('run_paired_epoch_timing_checks.py','test_paired_epoch_timing.py',
               'test_paired_epoch_analysis.py','freeze_paired_epoch_timing.py',
               'test_paired_epoch_timing_provenance.py','retain_paired_epoch_timing.py',
               'test_retain_paired_epoch_timing.py')
FINAL_PYTHON = RUN_TOOLS | set(FINAL_ROOTS) | {
    'paired_epoch_analysis.py','replay_paired_epoch_analysis.py','test_paired_r19932_screen.py',
    'run_paired_r19932_screen.py','replay_paired_r19932_screen.py','run_split_cache_screen.py'}
# Only these two original assets are read by the final adapter tests. Other
# generated/header/guard inputs remain bound in the unchanged build inventory.
FINAL_ASSETS = ('paired_timer_r19932.cpp','paired_timer_clock.cpp')


def executable(profile,variant):
    require(profile in ARCHIVES and variant in VARIANTS, 'executable selection')
    return ('measurement' if variant == 'steady' else 'qualification')+'/'+profile+'/'+variant


def qualification_args(variant,args):
    require(variant in VARIANTS, 'qualification variant')
    require(type(args) is list and args and all(type(a) is str for a in args), 'qualification args')
    # This is an execution firewall, not a claim based on stdout clock labels.
    require(variant != 'steady' or args[0] in ('--check','--exercise'),
            'real steady qualification is clock-free only')


def native_reference(root):
    folder = root/'native_reference'
    equal(sha(folder/'checks.json'),NATIVE_CHECKS_PIN)
    checks = parse((folder/'checks.json').read_text())
    files = {'checks.json'}; paths = {}
    for c in range(9):
        name = f'native-{c}-NNNN-1-check'
        rows = [r for r in checks['records'] if r['label']==name]
        require(len(rows)==1,'native reference row')
        path = folder/(name+'.parity'); files.add(path.name)
        equal(sha(path),rows[0]['parity_sha256'])
        equal(path.stat().st_size,CELLS[c][1]*CELLS[c][2])
        paths[c] = path
    equal(sorted(p.name for p in folder.iterdir()),sorted(files))
    return paths


def scope_command(root,profile,variant,args):
    qualification_args(variant,args)
    result = old_scope(root,profile,variant,args)
    result[-len(args)-1] = str(Path(root)/'build'/executable(profile,variant))
    return result


def inventory():
    rows = []
    def add(p,v,mode,c=8,s=None,g=1,code=0,fault=None,epoch_id=None,reason=None,args=None):
        if s is None: s = 'NNNN' if p == 'native' else '0110'
        rows.append(dict(label=f'{p}-{v}-{len(rows):03}',profile=p,variant=v,
                         args=args if args is not None else [mode,str(c),s,str(g)],code=code,
                         fault=fault,fault_epoch=epoch_id,reason=reason))
    for p in ARCHIVES:
        # All plain cells/APIs/profiles, plus grouped GF8. Reuse the old full
        # codec/guard matrix; these specifically qualify the new plain links.
        for c in range(9):
            for g in ((1,256) if c == 8 else (1,)):
                for mode in ('--check','--exercise'): add(p,'steady',mode,c,g=g)
                add(p,'fake-steady','--measure',c,g=g)
        # Clock-free plain-abort and witnessed-synthetic contrasts at target,
        # one-item batch, native1024-pointer geometry and grouped GF8 boundaries.
        for c,g in ((0,1),(1,1),(7,1),(8,256)):
            add(p,'abort','--exercise',c,g=g)
            add(p,'synthetic','--clock-exercise',c,g=g)
        if p != 'native':
            for s in ('1001','0000','1111'):
                for v,mode in (('steady','--exercise'),('fake-steady','--measure')):
                    add(p,v,mode,8,s=s,g=256)
        for e in range(3):
            for fault in ('equal','reverse','negative','huge'):
                add(p,'fake-steady','--measure',8,g=256,code=1,fault=fault,epoch_id=e)
        add(p,'abort','--clock-guard',code=86)
        for v,mode,reason in (
            ('abort','--measure','steady clock required'),
            ('synthetic','--measure','steady clock required'),
            ('fake-steady','--clock-exercise','synthetic clock required'),
            ('fake-steady','--clock-guard','abort clock required'),
            ('synthetic','--clock-guard','abort clock required')):
            add(p,v,mode,code=1,reason=reason)
        s = 'NNNN' if p == 'native' else '0110'
        for args,reason in ((['--check','9',s,'1'],'cell'),(['--check','0',s,'256'],'group'),
                            (['--measure','8',s,'256','forbidden'],'no measured or clock-guard parity dump')):
            add(p,'fake-steady',None,code=1,reason=reason,args=args)
    equal(len(rows),189)
    return rows


def verify_record(row,output,error,image):
    actual = epoch.records(output)
    p,v,args = row['profile'],row['variant'],row['args']
    witnessed = v in ('synthetic','fake-steady')
    wanted = {}
    calls,marks = epoch.public_oracle(p,0,'NNNN' if p=='native' else '0110',1,False,completed=0)
    progress = dict(schema=epoch.PROGRESS,selection_count=0,snapshot_count=0,timed=False)
    if row['reason']:
        equal(error,row['reason']+'\n')
        if witnessed: wanted[epoch.CLOCK] = epoch.clocks(1,0)
    else:
        mode,c,s,g = args[:4]; c,g = int(c),int(g)
        exercise = mode != '--check'
        if row['code'] == 0:
            equal(error,'')
            base = expected(p,c,s,g,exercise)
            per_epoch = dict(base)
            sampled = mode in ('--measure','--clock-exercise')
            base.update(schema=overlay.SCHEMA,clock_source='steady' if v in ('steady','fake-steady') else v,
                        encode_calls=base['encode_calls']*3,selections=base['selections']*3,
                        per_slot_calls=[n*3 for n in base['per_slot_calls']],timed=mode=='--measure',
                        samples=[[257+17*i,(257+17*i)/g] for i in range(252)] if sampled else [])
            for key in ('warmup_passes','exercise_passes'): base[key+'_per_epoch'] = base.pop(key)
            wanted[overlay.SCHEMA] = base
            wanted[epoch.ACCOUNT] = dict(schema=epoch.ACCOUNT,timed=False,epochs=[dict(epoch=e,
                calls=per_epoch['encode_calls'],selections=per_epoch['selections'],
                per_slot_calls=per_epoch['per_slot_calls'],probes=per_epoch['probes'],
                sample_begin=e*84 if sampled else 0,sample_count=84 if sampled else 0) for e in range(3)])
            require(epoch.META in actual,'metadata missing')
            epoch.validate_metadata(actual.pop(epoch.META),p,c,s,exercise,image)
            progress.update(selection_count=base['selections'],snapshot_count=6)
            calls,marks = epoch.public_oracle(p,c,s,g,exercise)
            if witnessed: wanted[epoch.CLOCK] = epoch.clocks(g)
        elif row['code'] == 86:
            equal(error,'unexpected driver benchmark clock\n')
            progress.update(selection_count=21,snapshot_count=1)
        else:
            equal(row['code'],1); require(row['fault_epoch'] in (0,1,2),'fault epoch')
            e = row['fault_epoch']
            equal(error,('group duration exceeds exact binary64 integer range' if row['fault']=='huge' else
                         'nonpositive or reversed grouped clock interval')+'\n')
            progress.update(selection_count=e*104+21,snapshot_count=2*e+1)
            calls,marks = epoch.public_oracle(p,c,s,g,True,completed=e,tail='fault')
            wanted[epoch.CLOCK] = epoch.clocks(g,168*e+2)
    if witnessed: wanted.update({epoch.WITNESS:calls,epoch.MARKS:marks})
    wanted[epoch.PROGRESS] = progress
    equal(actual,wanted)


def verify_build(root):
    out = root/'build'; state = parse((out/'build.json').read_text())
    equal([state['bead'],state['completed'],state['real_clocks_read'],state['reference_build_sha256']],
          [BEAD,True,False,REFERENCE_PIN])
    equal(sha(out/'reference_build.json'),REFERENCE_PIN)
    reference = parse((out/'reference_build.json').read_text())
    inputs = set(ASSETS)
    for p in ARCHIVES:
        inputs |= {p+'/'+n for n in REUSED}
        inputs |= {n for n in reference['artifacts'] if n.startswith(p+'/include/')}
    equal(state['inputs'],{n:reference['artifacts'][n] for n in inputs})
    names = {'baseline/'+n for n in ASSETS}|(set(ASSETS)-{'paired_metadata.cpp'})
    names |= inputs-set(ASSETS)
    names |= {'paired_epoch_diagnostic.cpp','fake_steady_clock.cpp','plain_marks.cpp','reference_build.json'}
    for name,digest in state['inputs'].items():
        equal(sha(out/('baseline/'+name if name in ASSETS else name)),digest)
    for name in set(ASSETS)-{'paired_metadata.cpp'}:
        equal(sha(out/name),state['inputs'][name])
    equal((out/'paired_epoch_diagnostic.cpp').read_text(),overlay.driver((out/'baseline/paired_metadata.cpp').read_text()))
    equal((out/'fake_steady_clock.cpp').read_text(),overlay.fake_steady_clock((out/'baseline/paired_timer_clock.cpp').read_text()))
    equal((out/'plain_marks.cpp').read_text(),overlay.PLAIN_MARK)
    commands,images = [],{}
    # Independently reconstruct each recipe from the pinned qualified driver
    # flags; never accept a recorded command merely because it produced a file.
    origin = str(Path(reference['root'])/'build')
    target = str(Path(state['root'])/'build')
    for p,digest in ARCHIVES.items():
        equal(sha(out/p/'codec.a'),digest)
        compile_ = [c for c in reference['commands'] if c[-1] == origin+'/'+p+'/driver.o']
        require(len(compile_) == 1,'reference compile identity')
        options = [a.replace(origin+'/',target+'/') for a in compile_[0][:-4]]
        obj = lambda n: target+'/'+p+'/'+n
        for filename,objectname,extra in (
            ('paired_epoch_diagnostic.cpp','driver.o',[]),('plain_marks.cpp','plain-marks.o',[]),
            ('paired_timer_clock.cpp','steady-clock.o',[]),
            ('fake_steady_clock.cpp','fake-steady-clock.o',['-DLEO_PAIRED_SYNTHETIC=1'])):
            commands.append(options+extra+['-c',target+'/'+filename,'-o',obj(objectname)])
            names.add(p+'/'+objectname)
        for v in VARIANTS:
            observed = v in ('synthetic','fake-steady')
            apis = ['leo_encode'] if p == 'native' else ['leo2_encode','leo2_encode_batch']
            flags = ['-Wl,--wrap='+api for api in apis] if observed else [
                '-Wl,--defsym=__real_'+api+'='+api for api in apis]
            command = options+[obj('driver.o'),obj('witness.o' if observed else 'plain-marks.o'),obj(v+'-clock.o')]
            if v == 'abort': command += [target+'/clock_guard.cpp']
            command += flags+([] if v == 'steady' else ['-Wl,--wrap='+CLOCK])
            name = executable(p,v); names.add(name)
            commands.append(command+[obj('codec.a'),'-ldl','-o',target+'/'+name])
            images[name] = elf(out/name,p=='native',v)
    equal(state['commands'],commands); equal(state['images'],images)
    epoch.verify_tools(out/'source_tools',state['tools'],BUILD_TOOLS,['build_paired_epoch_timing.py'])
    names |= {'source_tools/'+n for n in BUILD_TOOLS}
    equal(sorted(state['artifacts']),sorted(names))
    check_inventory(out,state['artifacts'],exclude=('build.json',))
    return state,images,build_resources((root/'build.log').read_text())


def verify_final_tools(root):
    """Bind final executing dependencies without rewriting collection tools."""
    final = parse((root/'final-tools.json').read_text())
    equal(sorted(final),['bead','files','real_clocks_read'])
    equal([final['bead'],final['real_clocks_read']],[BEAD,False])
    folder = root/'final_tools'
    require(folder.is_dir() and not folder.is_symlink(),'real final tools directory')
    equal(sorted(final['files']),sorted(FINAL_PYTHON | set(FINAL_ASSETS)))
    equal(sorted(p.name for p in folder.iterdir()),sorted(FINAL_PYTHON | set(FINAL_ASSETS)))
    check_inventory(folder,final['files'])
    equal(sorted(epoch.tool_closure(folder,FINAL_ROOTS)),sorted(FINAL_PYTHON))
    equal(sha(Path(__file__)),final['files']['verify_paired_epoch_timing.py'])
    for module_name,module in tuple(sys.modules.items()):
        source=getattr(module,'__file__',None)
        expected_name=module_name+'.py'; name=Path(source).name if source else None
        if expected_name not in FINAL_PYTHON and name not in FINAL_PYTHON: continue
        require(source is not None and name in FINAL_PYTHON,'executing dependency source: '+module_name)
        if expected_name in FINAL_PYTHON: equal(name,expected_name)
        require(sha(Path(source))==final['files'][name],'executing dependency hash: '+name)
    reference=root/'build/reference_build.json'; equal(sha(reference),REFERENCE_PIN)
    baseline=parse(reference.read_text())['baseline']
    for name in FINAL_ASSETS: equal(final['files'][name],baseline[name])
    return sha(root/'final-tools.json')


def replay(root, *, final=True):
    build,images,build_peak = verify_build(root)
    references = native_reference(root)
    state = parse((root/'checks/checks.json').read_text())
    equal([state['bead'],state['completed'],state['real_clocks_read'],state['root']],
          [BEAD,True,False,build['root']])
    epoch.verify_tools(root/'check_tools',state['tools'],RUN_TOOLS,['run_paired_epoch_timing_checks.py'])
    rows = inventory(); equal([r['label'] for r in state['records']],[r['label'] for r in rows])
    totals = dict(positive=0,clock_aborts=0,clock_faults=0,cli_refusals=0,parity_comparisons=0,parity_bytes=0)
    files,peaks = {'checks.json'},[]
    for row,record in zip(rows,state['records']):
        name = executable(row['profile'],row['variant']); label = row['label']
        args = list(row['args'])
        parity = root/'checks'/(label+'.parity') if row['code']==0 and row['args'][0]!='--measure' else None
        if parity: args.append(str(Path(build['root'])/'checks'/parity.name)); files.add(parity.name)
        equal(record['command'],scope_command(build['root'],row['profile'],row['variant'],args))
        equal([record['variant'],record['executable_sha256'],record['returncode'],record['fault'],record['fault_epoch']],
              [row['variant'],build['artifacts'][name],row['code'],row['fault'],row['fault_epoch']])
        stdout,stderr = root/'checks'/(label+'.stdout'),root/'checks'/(label+'.stderr')
        files |= {stdout.name,stderr.name}
        equal(sha(stdout),record['stdout_sha256']); equal(sha(stderr),record['stderr_sha256'])
        output,error,peak = split_scope(stdout.read_text(),stderr.read_text(),row['code'])
        equal(peak,record['memory_peak']); peaks.append(peak)
        verify_record(row,output,error,images[name])
        if row['code']==0:
            totals['positive'] += 1
            if parity:
                equal(sha(parity),record['parity_sha256'])
                reference = references[int(row['args'][1])]
                size = epoch.compare(parity,reference)
                equal(size,CELLS[int(row['args'][1])][1]*CELLS[int(row['args'][1])][2])
                totals['parity_comparisons'] += 1; totals['parity_bytes'] += size
        else: totals['cli_refusals' if row['reason'] else 'clock_aborts' if row['code']==86 else 'clock_faults'] += 1
    equal(sorted(p.name for p in (root/'checks').iterdir()),sorted(files))
    result = dict(bead=BEAD,real_clocks_read=False,default_enabled=False,diagnostic_only=True,
                  totals=totals,build_peak=build_peak,maximum_child_peak=max(peaks),collector_qualified=False)
    if final: result['final_tools_sha256']=verify_final_tools(root)
    return result


if __name__ == '__main__': print(json.dumps(replay(Path(sys.argv[1]).resolve(strict=True)),sort_keys=True))

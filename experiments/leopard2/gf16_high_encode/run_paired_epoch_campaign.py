#!/usr/bin/env python3
"""One separately preregistered diagnostic attempt. Never use for dry runs."""
import fcntl
import json
import os
from pathlib import Path
import subprocess
import sys
import stat

import paired_epoch_campaign as c
import paired_epoch_analysis as analysis
from run_split_cache_screen import host_identity,sibling_ticks,check_passive
from replay_paired_epoch_campaign import qualification_gate


def preregistration(bundle,commit):
    c.require(type(commit) is str and len(commit)==40 and all(v in '0123456789abcdef' for v in commit),'commit')
    for name in c.SOURCES | c.ASSETS | c.EXTRA:
        data=subprocess.check_output(['git','show',commit+':experiments/leopard2/gf16_high_encode/'+name],timeout=15)
        c.equal(data.decode(),(bundle/name).read_text())
    subprocess.run(['git','merge-base','--is-ancestor',commit,'origin/codex/claude-fable-5-1-audit'],check=True,timeout=15)


def scope_and_host():
    c.equal(host_identity(),c.HOST)
    for cpu in (26,90):
        c.equal(Path(f'/sys/devices/system/cpu/cpu{cpu}/topology/thread_siblings_list').read_text().strip(),'26,90')
    group=next(line.split(':',2)[2] for line in Path('/proc/self/cgroup').read_text().splitlines() if line.startswith('0::'))
    scope=Path('/sys/fs/cgroup')/group.lstrip('/')
    c.equal((scope/'memory.max').read_text().strip(),'268435456')
    c.equal((scope/'memory.swap.max').read_text().strip(),'0')
    c.equal((scope/'memory.events').read_text().splitlines(),['low 0','high 0','max 0','oom 0','oom_kill 0','oom_group_kill 0'])
    os.sched_setaffinity(0,{0})
    c.equal(sorted(os.sched_getaffinity(0)),[0])
    return dict(host=c.HOST,topology='26,90',controller_affinity=[0],memory_max=268435456,
                swap_max=0,initial_memory_events=[0]*6,cgroup=group)


def acquire_locks():
    lease=Path(f'/run/user/{os.getuid()}/leopard2-cpu-leases')
    c.require(lease.is_dir() and not lease.is_symlink() and lease.stat().st_uid==os.getuid()
              and lease.stat().st_mode&0o777==0o700,'lease directory')
    handles=[]
    try:
        for path in (Path('/tmp/leopard-gf8-authoritative.lock'),lease/f'leopard2-cpu-pair-{os.getuid()}-26-90.lock'):
            fd=os.open(path,os.O_RDONLY|os.O_CREAT|os.O_NOFOLLOW,0o600); handles.append(fd)
            info=os.fstat(fd)
            c.require(stat.S_ISREG(info.st_mode) and info.st_uid==os.getuid() and info.st_mode&0o777==0o600,
                      'lock owner/type/mode')
            fcntl.flock(fd,fcntl.LOCK_EX|fcntl.LOCK_NB)
        return handles
    except BaseException:
        for fd in reversed(handles): os.close(fd)
        raise


def launch(output,intents,label,command,env,observe=True):
    intent=dict(label=label,command=command,environment=env)
    intents.write(json.dumps(intent,sort_keys=True)+'\n'); intents.flush(); os.fsync(intents.fileno())
    record=dict(intent,returncode=None,timed_out=False,stdout=label+'.stdout',stderr=label+'.stderr',
                stdout_sha256=None,stderr_sha256=None,sibling_delta=None,failure=None)
    before=None
    try:
        before=sibling_ticks(90) if observe else None
        with (output/record['stdout']).open('xb') as out,(output/record['stderr']).open('xb') as err:
            try:
                record['returncode']=subprocess.run(command,stdout=out,stderr=err,env=env,timeout=60).returncode
            except subprocess.TimeoutExpired:
                record['timed_out']=True
    except BaseException as error:
        record['failure']=f'{type(error).__name__}: {error}'
    finally:
        try:
            if observe and before is not None: record['sibling_delta']=sibling_ticks(90)-before
            for field in ('stdout','stderr'):
                path=output/record[field]
                if path.is_file(): record[field+'_sha256']=c.sha(path)
        except BaseException as error:
            record['failure']=f'{type(error).__name__}: {error}'
        # A killed controller still leaves a flushed intent; an observed
        # failure also has its own outcome before any caller validation.
        with (output/(label+'.outcome.json')).open('x') as stream:
            json.dump(record,stream,sort_keys=True); stream.write('\n'); stream.flush(); os.fsync(stream.fileno())
    return record


def successful(output,row,quiet=True):
    c.equal([row['returncode'],row['timed_out'],row['failure']],[0,False,None])
    c.equal(c.bounded_text(output/row['stderr']),'')
    if quiet: c.equal(row['sibling_delta'],0)


def run(bundle,commit):
    # No attempt reservation, passive clock, host mutation or subprocess before
    # the immutable readiness and published source gates have passed.
    value,pins,images=c.inputs(bundle)
    c.validate_plan(value,timing=True); c.equal(str(bundle),value['timing_path'])
    qualification=qualification_gate(pins)
    preregistration(bundle,commit)
    output=Path(value['attempt_root']); output.mkdir(mode=0o700)
    state=dict(schema='leopard-epoch-campaign-attempt/v1',bead=c.BEAD,preregistration=commit,
               bundle_path=str(bundle),plan_sha256=pins['files'][c.PLAN],pins=pins,host=c.HOST,
               child_environment=c.ENV,preflight=[],invocations=[],conditions=[],complete=False,qualification=qualification)
    handles=[]
    try:
        handles=acquire_locks(); state['scope']=scope_and_host(); c.inputs(bundle)
        with (output/'launches.jsonl').open('x') as intents:
            def condition(label):
                record=launch(output,intents,label,['/bin/bash',str(bundle/'paired_epoch_condition.sh')],c.condition_env(),False)
                state['conditions'].append(record); successful(output,record,False)
                c.equal(c.bounded_text(output/record['stdout']).splitlines(),c.condition_lines())
            def invoke(item,measured):
                c.inputs(bundle)
                record=launch(output,intents,c.name(item),c.command(bundle,item['cell'],item['order'],measured),c.ENV)
                target=state['invocations' if measured else 'preflight']; target.append(dict(item,**record))
                successful(output,record,measured); c.inputs(bundle)
                profile='native' if item['order']=='NNNN' else 'release'
                text=c.bounded_text(output/record['stdout'])
                if measured: analysis.projections(c.qualified.epoch.records(text),item['cell'],item['order'],images[profile])
                else: c.check_record(text,item['cell'],item['order'],images[profile])
            condition('condition-before')
            for item in c.preflights(): invoke(item,False)
            check_passive(state,value)
            for item in analysis.inherited.schedule(): invoke(item,True)
            condition('condition-after')
        c.inputs(bundle)
        entries=[{k:r[k] for k in (*item,'sibling_delta','stdout','stdout_sha256')}
                 for item,r in zip(analysis.inherited.schedule(),state['invocations'])]
        state['analysis']=analysis.analyze(output,entries,images)
        c.executing_sources(pins)
        for fd in handles: os.fstat(fd)  # Descriptor leases remain held through analysis.
        c.require(len(handles)==2,'both campaign locks held')
        state['complete']=True
    except BaseException as error:
        state['complete']=False; state['failure']=f'{type(error).__name__}: {error}'; raise
    finally:
        try:
            with (output/'attempt.json').open('x') as stream:
                json.dump(state,stream,sort_keys=True,allow_nan=False); stream.write('\n'); stream.flush(); os.fsync(stream.fileno())
        finally:
            for fd in reversed(handles): os.close(fd)


if __name__=='__main__':
    c.require(len(sys.argv)==3,'usage: run_paired_epoch_campaign.py PREREGISTERED_FROZEN PUSHED_COMMIT')
    run(Path(sys.argv[1]).resolve(strict=True),sys.argv[2])

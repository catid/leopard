#!/usr/bin/env python3
"""Private campaign copies and clock-free checks. Never invokes --measure."""
import json
import os
from pathlib import Path
import shutil
import sys

import paired_epoch_campaign as c
import run_paired_epoch_campaign as collector


def freeze():
    source=Path(__file__).resolve().parent
    value=c.read(source/c.PLAN); c.validate_plan(value)
    output=Path(value['timing_path' if value['ready_to_time'] else 'qualification_path'])
    if value['ready_to_time']:
        prior=c.inputs(Path(value['qualification_path']))[1]
        collector.qualification_gate(prior)
    copy_bundle(output,source)
    if value['ready_to_time']:
        c.equal({k:v for k,v in c.inputs(output)[1]['files'].items() if k!=c.PLAN},
                {k:v for k,v in prior['files'].items() if k!=c.PLAN})
    return output


def copy_bundle(output,source):
    """Also used for disposable synthetic test copies; creates no attempt."""
    c.disjoint_attempt(output)
    history=c.manifest(c.QUALIFIED/'SHA256SUMS')
    # Pin the completed proof and old dependencies before creating a new copy.
    sources={name:source/name for name in c.NEW_SOURCES | c.EXTRA}
    sources.update({name:c.QUALIFIED/'final_tools'/name for name in c.qualified.FINAL_PYTHON | c.ASSETS})
    sources.update({name:c.QUALIFIED/path for name,path in c.PROOFS.items()})
    sources.update(native=c.QUALIFIED/'build/measurement/native/steady',
                   current=c.QUALIFIED/'build/measurement/release/steady')
    c.equal(sorted(sources),sorted(c.FILES))
    for name,path in sources.items():
        c.require(path.is_file() and not path.is_symlink(),'regular freeze source: '+name)
        if name in c.qualified.FINAL_PYTHON | c.ASSETS:
            c.equal(c.sha(path),history['final_tools/'+name])
    output.mkdir(mode=0o700)
    pins={}
    for name,path in sorted(sources.items()):
        before=c.sha(path); destination=output/name
        with path.open('rb') as src,destination.open('xb') as dst:
            shutil.copyfileobj(src,dst,65536)
        c.equal(c.sha(destination),before); c.equal(c.sha(path),before)
        c.require(destination.stat().st_ino!=path.stat().st_ino or destination.stat().st_dev!=path.stat().st_dev,
                  'private copy')
        destination.chmod(0o555 if name in ('native','current') else 0o444); pins[name]=before
    with (output/'pins.json').open('x') as stream:
        json.dump(dict(schema='leopard-epoch-campaign-pins/v1',files=pins),stream,sort_keys=True,indent=2)
        stream.write('\n')
    (output/'pins.json').chmod(0o444); output.chmod(0o555)
    c.inputs(output)
    return output


def check(bundle,output):
    value,pins,images=c.inputs(bundle)
    c.equal(value['ready_to_time'],False)
    c.equal(str(bundle),value['qualification_path'])
    c.disjoint_attempt(output)
    c.equal(str(output),str(output.resolve()))
    c.equal(str(output),value['qualification_output'])
    output.mkdir(mode=0o700)
    state=dict(schema='leopard-epoch-campaign-checks/v1',bead=c.BEAD,bundle_path=str(bundle),output_path=str(output),
               plan_sha256=pins['files'][c.PLAN],pins=pins,host=c.HOST,child_environment=c.ENV,
               preflight=[],conditions=[],complete=False,real_clocks_read=False)
    handles=[]
    try:
        handles=collector.acquire_locks(); state['scope']=collector.scope_and_host()
        with (output/'launches.jsonl').open('x') as intents:
            def condition(label):
                record=collector.launch(output,intents,label,['/bin/bash',str(bundle/'paired_epoch_condition.sh')],
                                        c.condition_env(),False)
                state['conditions'].append(record); collector.successful(output,record,False)
                c.equal(c.bounded_text(output/record['stdout']).splitlines(),c.condition_lines())
            condition('condition-before')
            for item in c.preflights():
                c.inputs(bundle)
                args=c.command(bundle,item['cell'],item['order'],False)
                c.require('--measure' not in args,'clock-free qualification')
                record=collector.launch(output,intents,c.name(item),args,c.ENV)
                state['preflight'].append(dict(item,**record)); collector.successful(output,record,False)
                profile='native' if item['order']=='NNNN' else 'release'
                c.check_record(c.bounded_text(output/record['stdout']),item['cell'],item['order'],images[profile])
            condition('condition-after')
        c.inputs(bundle)
        c.require(len(handles)==2,'both locks held')
        for fd in handles: os.fstat(fd)
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
    if sys.argv[1:]==['freeze']: print(freeze())
    else:
        c.require(len(sys.argv)==4 and sys.argv[1]=='check',
                  'usage: freeze_paired_epoch_campaign.py freeze | check QUALIFICATION_BUNDLE FRESH_OUTPUT')
        check(Path(sys.argv[2]).resolve(strict=True),Path(sys.argv[3]).absolute())

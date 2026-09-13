#!/usr/bin/env python3
"""Full journal/raw-file validation; does not import the campaign collector."""
import json
from pathlib import Path
import re
import subprocess
import sys

import paired_epoch_campaign as c
import replay_paired_epoch_analysis as analysis


def publication(bundle,commit):
    c.require(type(commit) is str and re.fullmatch('[0-9a-f]{40}',commit) is not None,'preregistration commit')
    subprocess.run(['git','merge-base','--is-ancestor',commit,'origin/codex/claude-fable-5-1-audit'],
                   check=True,timeout=15)
    for name in sorted(c.SOURCES | c.ASSETS | c.EXTRA):
        actual=subprocess.check_output(['git','show',commit+':experiments/leopard2/gf16_high_encode/'+name],timeout=15)
        c.equal(actual.decode(),c.bounded_text(bundle/name,4*1024**2))


def verify(root,bundle,resource_log,*,qualification=False):
    value,pins,images=c.inputs(bundle)
    c.require(type(qualification) is bool,'qualification mode')
    c.equal(value['ready_to_time'],not qualification)
    state=c.read(root/'attempt.json')
    keys={'schema','bead','bundle_path','plan_sha256','pins','host','child_environment','preflight',
          'conditions','complete','scope'}
    keys |= {'real_clocks_read','output_path'} if qualification else {'preregistration','invocations','passive','analysis','qualification'}
    c.equal(sorted(state),sorted(keys))
    c.equal(state['schema'],'leopard-epoch-campaign-'+('checks' if qualification else 'attempt')+'/v1')
    c.equal([state['bead'],state['complete'],state['pins'],state['plan_sha256'],state['host'],state['child_environment']],
            [c.BEAD,True,pins,pins['files'][c.PLAN],c.HOST,c.ENV])
    logical=Path(value['qualification_path' if qualification else 'timing_path'])
    c.equal(state['bundle_path'],str(logical))
    scope=state['scope']; c.equal(sorted(scope),sorted(['host','topology','controller_affinity','memory_max','swap_max',
                                            'initial_memory_events','cgroup']))
    c.equal({k:v for k,v in scope.items() if k!='cgroup'},dict(host=c.HOST,topology='26,90',
            controller_affinity=[0],memory_max=268435456,swap_max=0,initial_memory_events=[0]*6))
    c.require(type(scope['cgroup']) is str and scope['cgroup'].startswith('/') and
              '..' not in Path(scope['cgroup']).parts and scope['cgroup'].endswith('.scope'),'resource cgroup')
    original_output=Path(value['attempt_root'])
    if qualification:
        c.require(type(state['output_path']) is str,'qualification output path')
        original_output=Path(state['output_path'])
        c.require(original_output.is_absolute() and '..' not in original_output.parts and
                  str(original_output)==state['output_path'],'canonical qualification path')
        c.disjoint_attempt(original_output)
        c.equal(str(original_output),value['qualification_output'])
    peak=c.resource(resource_log,scope['cgroup'],
        c.controller_command(logical,original_output,qualification,state.get('preregistration')),
        logical/'paired_epoch_scope.sh')
    if qualification: c.equal(state['real_clocks_read'],False)
    else:
        c.equal(state['qualification'],qualification_gate(pins))
        publication(bundle,state['preregistration'])
        passive=state['passive']; c.equal(sorted(passive),sorted(['before','after','elapsed_ns']))
        c.require(all(type(v) is int and v>=0 for v in passive.values()),'passive integer fields')
        c.equal(passive['before'],passive['after']); c.require(passive['elapsed_ns']>=10**10,'full passive observation')
    expected_files={'attempt.json','launches.jsonl'}; intents=[]
    outcome_keys={'label','command','environment','returncode','timed_out','stdout','stderr','stdout_sha256',
                  'stderr_sha256','sibling_delta','failure'}
    def outcome(row,label,command,env,quiet,item=None):
        c.equal(sorted(row),sorted(outcome_keys | (set(item) if item is not None else set())))
        if item is not None: c.equal({k:row[k] for k in item},item)
        c.equal([row['label'],row['command'],row['environment'],row['returncode'],row['timed_out'],row['failure']],
                [label,command,env,0,False,None])
        if quiet: c.equal(row['sibling_delta'],0)
        elif item is None: c.equal(row['sibling_delta'],None)
        else: c.require(type(row['sibling_delta']) is int and row['sibling_delta']>=0,'preflight sibling count')
        texts={}
        for field in ('stdout','stderr'):
            filename=label+'.'+field; c.equal(row[field],filename); expected_files.add(filename)
            path=root/filename; texts[field]=c.bounded_text(path); c.equal(c.sha(path),row[field+'_sha256'])
        c.equal(texts['stderr'],'')
        filename=label+'.outcome.json'; expected_files.add(filename)
        c.equal(c.read(root/filename),{k:row[k] for k in outcome_keys})
        intents.append(dict(label=label,command=command,environment=env))
        return texts['stdout']
    c.equal(len(state['conditions']),2); c.equal(len(state['preflight']),27)
    def condition(index):
        label='condition-'+('before' if index==0 else 'after')
        text=outcome(state['conditions'][index],label,['/bin/bash',str(logical/'paired_epoch_condition.sh')],
                     c.condition_env(),False)
        c.equal(text.splitlines(),c.condition_lines())
    condition(0)
    for index,(cell,order) in enumerate((cell,order) for cell in range(9) for order in ('NNNN','0000','1111')):
        item=dict(cell=cell,order=order); label=f'check-{cell}-{order}'
        text=outcome(state['preflight'][index],label,c.command(logical,cell,order,False),c.ENV,False,item)
        c.check_record(text,cell,order,images['native' if order=='NNNN' else 'release'])
    entries=[]
    if not qualification:
        c.equal(len(state['invocations']),318)
        for row,item in zip(state['invocations'],analysis.schedule()):
            label=f"cell-{item['cell']}-round-{item['round']}-{item['comparison']}-slot-{item['slot']}"
            outcome(row,label,c.command(logical,item['cell'],item['order'],True),c.ENV,True,item)
            entries.append({k:row[k] for k in (*item,'sibling_delta','stdout','stdout_sha256')})
    condition(1)
    actual=[c.parse(line) for line in c.bounded_text(root/'launches.jsonl').splitlines()]
    c.equal(actual,intents)
    c.equal(sorted(p.name for p in root.iterdir()),sorted(expected_files))
    for name in expected_files:
        c.require((root/name).is_file() and not (root/name).is_symlink(),'regular campaign inventory')
    result=dict(bead=c.BEAD,qualification=qualification,preflights=27,timed_processes=len(entries),
                diagnostic_only=True,production_promotion=False,default_enabled=False,memory_peak=peak)
    if not qualification:
        result['analysis']=analysis.derive(root,entries,images)
        analysis.compare(state['analysis'],result['analysis'])
    c.inputs(bundle)
    return result


def qualification_gate(pins):
    """Authenticate the actual immutable 27-check qualification, not a flag."""
    bundle=Path(c.plan()['qualification_path'])
    value,qualified_pins,_=c.inputs(bundle)
    c.equal(value['ready_to_time'],False)
    c.equal({k:v for k,v in pins['files'].items() if k!=c.PLAN},
            {k:v for k,v in qualified_pins['files'].items() if k!=c.PLAN})
    root=Path(value['qualification_output']); resource=Path(value['qualification_resource'])
    c.require(root.is_dir() and not root.is_symlink() and not root.stat().st_mode&0o222,'readonly qualification')
    for path in [resource,*root.iterdir()]:
        c.require(path.is_file() and not path.is_symlink() and not path.stat().st_mode&0o222,'readonly qualification record')
    return verify(root,bundle,resource,qualification=True)


if __name__=='__main__':
    c.require(len(sys.argv) in (4,5) and (len(sys.argv)==4 or sys.argv[4]=='--qualification'),
              'usage: replay_paired_epoch_campaign.py ATTEMPT BUNDLE RESOURCE_LOG [--qualification]')
    print(json.dumps(verify(Path(sys.argv[1]),Path(sys.argv[2]),Path(sys.argv[3]),
                            qualification=len(sys.argv)==5),sort_keys=True,allow_nan=False))

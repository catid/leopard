"""Independent three-epoch projection/estimation; imports no collector."""
import math
from pathlib import Path

import replay_paired_r19932_screen as prior
import verify_paired_epoch as epoch
from verify_paired_metadata import sha, require, equal
from verify_paired_r19932 import expected

BEAD = 'leopard-79h.38.5.4.19.1.4.4'
SCHEMA = 'leopard-paired-epoch-diagnostic/v1'


def schedule():
    for cell in range(9):
        for rnd in range(3):
            comparisons = [('paired_0110',['0110']),('paired_1001',['1001']),
                           ('same_off',['0000']*4),('same_on',['1111']*4)]
            if cell<2: comparisons += [('native_on',['NNNN','1111','1111','NNNN']),('same_native',['NNNN']*4)]
            for comp,orders in comparisons:
                for slot,order in enumerate(orders):
                    yield dict(cell=cell,round=rnd,comparison=comp,slot=slot,order=order)


def projections(records,cell,order,image):
    equal(sorted(records),sorted([SCHEMA,epoch.ACCOUNT,epoch.META,epoch.PROGRESS]))
    profile = 'native' if order=='NNNN' else 'release'
    count = 6401 if cell==8 else 26
    base = expected(profile,cell,order,256 if cell==8 else 1,True)
    full = dict(base,schema=SCHEMA,clock_source='steady',timed=True,encode_calls=count*12,
                selections=312,per_slot_calls=[count*3]*4,warmup_passes_per_epoch=4,exercise_passes_per_epoch=21)
    full.pop('warmup_passes'); full.pop('exercise_passes')
    actual = records[SCHEMA]
    equal({k:v for k,v in actual.items() if k!='samples'},full)
    samples = actual['samples']
    require(type(samples) is list and len(samples)==252,'three complete sample epochs')
    group = 256 if cell==8 else 1
    for sample in samples:
        require(type(sample) is list and len(sample)==2,'sample shape')
        duration,average = sample
        require(type(duration) is int and duration in range(1,2**53),'exact duration')
        require(type(average) is float and math.isfinite(average) and average*group==duration,'normalized duration')
    phases = []
    for e in (0,1,2):
        phases.append(dict(epoch=e,calls=count*4,selections=104,per_slot_calls=[count]*4,
                           probes=base['probes'],sample_begin=84*e,sample_count=84))
    equal(records[epoch.ACCOUNT],dict(schema=epoch.ACCOUNT,timed=False,epochs=phases))
    equal(records[epoch.PROGRESS],dict(schema=epoch.PROGRESS,selection_count=312,snapshot_count=6,timed=False))
    epoch.validate_metadata(records[epoch.META],profile,cell,order,True,image)
    rows = []
    for e in (0,1,2):
        row = dict(actual)
        row.update(schema='leopard-paired-timer-r19932/v1',encode_calls=count*4,selections=104,
                   per_slot_calls=[count]*4,warmup_passes=4,exercise_passes=21,samples=samples[84*e:84*e+84])
        row.pop('warmup_passes_per_epoch'); row.pop('exercise_passes_per_epoch')
        prior.record(row,cell,order,True)
        rows.append(row)
    return rows


def derive_projected(epochs):
    require(type(epochs) is list and len(epochs)==3,'three epochs required')
    results = []
    for index,rows in enumerate(epochs):
        result = prior.derive(rows)
        if result['decision']=='qualify_default_on_artifact': result['decision']='diagnostic_all_gates_pass'
        result.update(epoch=index,diagnostic_only=True)
        results.append(result)
    homogeneous = []
    for i,item in enumerate(schedule()):
        if item['order'] not in ('0000','1111','NNNN'): continue
        costs = [prior.median([v[1] for v in records[i]['record']['samples']]) for records in epochs]
        homogeneous.append(dict(item,median_ns_per_call=costs,epoch1_over_epoch0=costs[1]/costs[0],
                                epoch2_over_epoch0=costs[2]/costs[0]))
    return dict(bead=BEAD,diagnostic_only=True,production_promotion=False,default_enabled=False,
                epochs=results,homogeneous_processes=homogeneous,
                all_epochs_pass=all(r['decision']=='diagnostic_all_gates_pass' for r in results),
                historical_shift_cause_established=False)


def derive(root,entries,images):
    require(type(entries) is list and len(entries)==318,'complete process inventory')
    epochs = [[],[],[]]
    for entry,item in zip(entries,schedule()):
        equal(sorted(entry),sorted([*item,'sibling_delta','stdout','stdout_sha256']))
        equal({k:entry[k] for k in item},item); equal(entry['sibling_delta'],0)
        name = f"cell-{item['cell']}-round-{item['round']}-{item['comparison']}-slot-{item['slot']}.stdout"
        equal(entry['stdout'],name)
        path = Path(root)/name
        require(path.is_file() and not path.is_symlink() and path.stat().st_size<1024**2,'bounded stdout')
        equal(sha(path),entry['stdout_sha256'])
        records = epoch.records(path.read_text())
        equal(sha(path),entry['stdout_sha256'])
        rows = projections(records,item['cell'],item['order'],images['native' if item['order']=='NNNN' else 'release'])
        for i,row in enumerate(rows): epochs[i].append(dict(item,sibling_delta=0,record=row))
    return derive_projected(epochs)


def compare(actual,derived):
    equal({k:v for k,v in actual.items() if k!='epochs'},{k:v for k,v in derived.items() if k!='epochs'})
    equal(len(actual['epochs']),3)
    for left,right in zip(actual['epochs'],derived['epochs']): prior.compare_analysis(left,right)

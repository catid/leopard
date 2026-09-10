"""Three-epoch diagnostic analysis. No process launch or clock acquisition.

Read one full process record at a time, validate all epochs before projection,
and retain only the small timing projections in memory. Never average epochs.
"""
import math
from pathlib import Path
import statistics

import run_paired_r19932_screen as inherited
import verify_paired_epoch as epoch
from paired_epoch_timing_overlay import BEAD, SCHEMA
from verify_paired_metadata import sha, equal, require
from verify_paired_r19932 import expected


def projections(records,cell,order,image):
    require(type(records) is dict and set(records)=={SCHEMA,epoch.META,epoch.ACCOUNT,epoch.PROGRESS},
            'plain diagnostic record inventory')
    profile = 'native' if order=='NNNN' else 'release'
    group = 256 if cell==8 else 1
    base = expected(profile,cell,order,group,True)
    fixed = dict(base,schema=SCHEMA,clock_source='steady',timed=True,
                 encode_calls=base['encode_calls']*3,selections=312,
                 per_slot_calls=[v*3 for v in base['per_slot_calls']])
    for key in ('warmup_passes','exercise_passes'): fixed[key+'_per_epoch'] = fixed.pop(key)
    row = records[SCHEMA]
    equal({k:v for k,v in row.items() if k!='samples'},fixed)
    samples = row['samples']
    require(type(samples) is list and len(samples)==252,'complete epoch sample inventory')
    for sample in samples:
        require(type(sample) is list and len(sample)==2,'sample pair')
        value,average = sample
        require(type(value) is int and 0<value<2**53,'sample integer interval')
        require(type(average) is float and math.isfinite(average) and average==value/group,
                'sample normalization')
    equal(records[epoch.ACCOUNT],dict(schema=epoch.ACCOUNT,timed=False,epochs=[dict(epoch=e,
        calls=base['encode_calls'],selections=104,per_slot_calls=base['per_slot_calls'],probes=base['probes'],
        sample_begin=e*84,sample_count=84) for e in range(3)]))
    equal(records[epoch.PROGRESS],dict(schema=epoch.PROGRESS,selection_count=312,snapshot_count=6,timed=False))
    epoch.validate_metadata(records[epoch.META],profile,cell,order,True,image)
    # No projection is produced until the complete raw shape has passed.
    return [dict(inherited.expected_record(cell,order,True),samples=samples[e*84:(e+1)*84]) for e in range(3)]


def read_process(root,entry):
    name = inherited.label(entry)+'.stdout'
    equal(entry['stdout'],name)
    path = Path(root)/name
    require(path.is_file() and not path.is_symlink() and path.stat().st_size<1024**2,'bounded raw stdout')
    equal(sha(path),entry['stdout_sha256'])
    records = epoch.records(path.read_text())
    equal(sha(path),entry['stdout_sha256'])
    return records


def analyze_projected(epochs):
    require(type(epochs) is list and len(epochs)==3,'three fixed epochs')
    results = []
    for index,rows in enumerate(epochs):
        result = inherited.analyze(rows)
        if result['decision']=='qualify_default_on_artifact': result['decision']='diagnostic_all_gates_pass'
        result.update(epoch=index,diagnostic_only=True)
        results.append(result)
    homogeneous = []
    for i,item in enumerate(inherited.schedule()):
        if len(set(item['order'])) != 1: continue
        costs = [statistics.median(s[1] for s in rows[i]['record']['samples']) for rows in epochs]
        homogeneous.append(dict(item,median_ns_per_call=costs,epoch1_over_epoch0=costs[1]/costs[0],
                                epoch2_over_epoch0=costs[2]/costs[0]))
    return dict(bead=BEAD,diagnostic_only=True,production_promotion=False,default_enabled=False,
                epochs=results,homogeneous_processes=homogeneous,
                all_epochs_pass=all(r['decision']=='diagnostic_all_gates_pass' for r in results),
                historical_shift_cause_established=False)


def analyze(root,entries,images):
    require(type(entries) is list and len(entries)==318,'complete ordered process inventory')
    epochs = [[],[],[]]
    for entry,item in zip(entries,inherited.schedule()):
        equal(sorted(entry),sorted([*item,'sibling_delta','stdout','stdout_sha256']))
        equal({k:entry[k] for k in item},item); equal(entry['sibling_delta'],0)
        records = read_process(root,entry)
        projected = projections(records,item['cell'],item['order'],images['native' if item['order']=='NNNN' else 'release'])
        for index,row in enumerate(projected): epochs[index].append(dict(item,sibling_delta=0,record=row))
    return analyze_projected(epochs)

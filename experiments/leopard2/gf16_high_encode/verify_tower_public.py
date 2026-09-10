#!/usr/bin/env python3
"""Collector-free tower frontend oracle. No codecs, clocks or subprocesses run."""
import copy
import json
from pathlib import Path
import re
import sys

import verify_avx2_adjacent_public as previous
from build_tower_public import ARCHIVES, BEAD, require, sha
from tower_public_overlay import adapt
from verify_gf16_callback_probe import model_shape
from verify_paired_r19932 import equal, parse, read

CELLS = previous.CELLS
SELECTED = (0,1,2,4,8)
FIELDS = ('selected_passes','source_rows','source_bytes','output_rows','output_bytes',
          'inverse_pairs','forward_pairs','accumulating_pairs')
PRIOR_MANIFEST = 'f3e58727918ded2e06b62c250ff1b19bde5d7860c107dbf9a8166a67fa0f4dfa'


def base_profile(profile):
    require(profile in ARCHIVES,'profile')
    return 'original' if profile == 'original-sanitize' else profile


def schedules(profile):
    return previous.schedules(base_profile(profile))


def valid(profile, cell, schedule, group):
    previous.valid(base_profile(profile),cell,schedule,group)


def structural(cell):
    """Count selected passes and actual tower callback expansion, NOT time shares.

    Unlike canonical range implementations, tower range wrappers invoke all
    four edges, including zero-skew edges. The external traversal is unchanged.
    """
    require(type(cell) is int and 0 <= cell < len(CELLS),'cell')
    if cell not in SELECTED: return [0]*8
    k,r,b,_,_ = CELLS[cell]
    tile,passes = 32768,b//32768
    traversal = model_shape((k,r,b,3,tile,passes))
    inverse = forward = accumulating = 0
    for (op,distance,_mask,_hint,_bytes),calls in traversal.items():
        if op == 'ifft2': inverse += calls
        elif op == 'fft2': forward += calls
        elif op == 'ifft2_xor': accumulating += calls
        elif op == 'ifft4_range': inverse += 4*distance*calls
        elif op == 'fft4_range': forward += 4*distance*calls
    return [passes,k*passes,k*b,r*passes,r*b,inverse,forward,accumulating]


def counts(profile, cell, states, initialized):
    operations = sum(s == '1' for s in states) if profile in ('trace','sanitize') else 0
    return dict(values=[n*operations for n in structural(cell)],initializations=int(initialized))


def expected(profile, cell, schedule, group, mode, clock):
    valid(profile,cell,schedule,group)
    # Reuse qualified API/workload/clock accounting, not the old kernel oracle.
    base = base_profile(profile)
    if base in ('trace','sanitize'): base = 'release'
    row = previous.expected(base,cell,schedule,group,mode,clock)
    row.update(schema='leopard-tower-public/v1',codec=profile+':'+ARCHIVES[profile],
               traced=profile in ('trace','sanitize'))
    del row['pair_probes'],row['pair_totals']
    runtime = profile in ('release','trace','sanitize')
    row['tower_probes'] = [counts(profile,cell,state,
        runtime and cell in SELECTED and '1' in schedule[:i+1]) for i,state in enumerate(schedule)]
    row['tower_totals'] = counts(profile,cell,schedule*(25*group if mode!='--check' else 0),
        runtime and cell in SELECTED and '1' in schedule)
    return row


def witness(profile, *args):
    return previous.witness(base_profile(profile),*args)


def inventory():
    rows = previous.inventory()
    # Retain all nine cases, all four runtime schedules, both GF8 group sizes,
    # abort/synthetic clock tests and malformed requests. Add original sanitizer.
    original = [copy.deepcopy(r) for r in rows if r['profile'] == 'original']
    for row in original:
        row['profile'] = 'original-sanitize'
        row['label'] = row['label'].replace('original-','original-sanitize-',1)
    at = next(i for i,r in enumerate(rows) if r['profile'] == 'release')
    rows[at:at] = original
    for profile,case,code in (('release','canary',0),('sanitize','canary',0),
                              ('sanitize','underflow',1),('sanitize','overflow',1)):
        rows.append(dict(label=profile+'-guard-'+case,profile=profile,binary='guards',
                         args=[case],code=code,parity=False,fault=None))
    return rows


def verify_record(row, output, error):
    if row['binary'] == 'guards':
        if row['code'] == 0:
            equal(error,''); equal(output,'both canaries rejected; zero-size and restored buffers pass\n')
        else:
            equal(output,'')
            require('ERROR: AddressSanitizer: use-after-poison' in error and
                    'READ of size 1' in error,'actual poisoned read refusal')
        return
    lines = [parse(line) for line in output.splitlines()]
    actual = {r['schema']:r for r in lines}
    require(len(lines) == len(actual),'duplicate schema')
    wanted = []
    p,binary,args = row['profile'],row['binary'],row['args']
    if '-bad-' in row['label']:
        require(bool(error),'CLI diagnostic')
    elif binary == 'group-unit':
        equal(error,''); wanted = [dict(schema='paired-group-unit/v1',cases=37,timed=False)]
    else:
        mode,c,s,g = args; c,g = int(c),int(g)
        if row['code'] == 0:
            equal(error,'')
            clock = 'steady' if binary == 'plain' else binary
            wanted = [expected(p,c,s,g,mode,clock)]
            if binary != 'plain': wanted += [witness(p,c,s,g)]
            if mode == '--clock-exercise': wanted += [previous.clocks(g)]
        elif row['code'] == 86:
            equal(error,'unexpected driver benchmark clock\n')
            wanted = [witness(p,c,s,g,4)]
        else:
            equal(error,('group duration exceeds exact binary64 integer range' if row['fault']=='huge'
                         else 'nonpositive or reversed grouped clock interval')+'\n')
            wanted = [witness(p,c,s,g,4,g),previous.clocks(g,1)]
    equal(actual,{r['schema']:r for r in wanted})


def scope_footer(output, code, maximum):
    """Independently inspect actual child exit marker and all resource counters."""
    require(output.count('TOWER_CHILD_EXIT=') == 1,'one scope exit record')
    native,resource = output.split('TOWER_CHILD_EXIT=')
    lines = resource.splitlines()
    equal(lines[0],str(code))
    equal(lines[1],'memory.peak')
    peak = int(lines[2])
    require(0 < peak < maximum,'scope memory peak')
    equal(lines[3:],['memory.max',str(maximum),'memory.events','low 0','high 0','max 0',
        'oom 0','oom_kill 0','oom_group_kill 0','memory.swap.current','0','memory.swap.max','0'])
    return native,peak


def split_scope(output, error, code):
    native,peak = scope_footer(output,code,256*1024**2)
    header,stderr = error.split('\n',1)
    require(re.fullmatch(r'Running as unit: run-[\w-]+\.scope; invocation ID: [0-9a-f]{32}',header)
            is not None,'systemd scope header')
    return native,stderr,peak


def input_snapshot(root, name, digest):
    """Resolve every build input to an exact retained copy, not mutable live paths."""
    path = Path(name)
    if path.name == 'build_tower_public.py':
        candidate = root/'tools'/path.name
    elif path.parent == Path('/home/catid/leopard/experiments/leopard2/gf16_high_encode'):
        candidate = root/'build'/path.name
    elif path.suffix == '.a':
        profiles = [p for p,d in ARCHIVES.items() if d == digest]
        require(len(profiles)==1,'known input archive')
        candidate = root/'build'/profiles[0]/'codec.a'
    elif path.name == 'clock_guard.cpp':
        candidate = root/'build/clock_guard.cpp'
    elif path.suffix == '.h':
        headers = 'native-headers' if '/pure-checks-v2/source/' in name else 'l2-headers'
        candidate = root/'build'/headers/path.name
    else:
        raise ValueError('unknown retained build input: '+name)
    equal(sha(candidate),digest)
    return candidate


def replay(root):
    build,checks = read(root/'build/build.json'),read(root/'checks/checks.json')
    for state in (build,checks): equal([state['bead'],state['completed'],state['timed']],[BEAD,True,False])
    for name,digest in build['inputs'].items(): input_snapshot(root,name,digest)
    for name,digest in build['artifacts'].items(): equal(sha(root/'build'/name),digest)
    for base,kind,target in (('avx2_adjacent_public.cpp','driver','tower_public.cpp'),
                             ('avx2_adjacent_public_witness.cpp','witness','tower_public_witness.cpp')):
        equal((root/'build'/target).read_text(),adapt((root/'build'/base).read_text(),kind))
    original_root = Path(build['root'])
    for profile,digest in ARCHIVES.items():
        equal(sha(root/'build'/profile/'codec.a'),digest)
        common = 'native' if profile=='native' else 'sanitize' if 'sanitize' in profile else 'l2-release'
        for binary in ('plain','synthetic','abort'):
            links = [c for c in build['commands'] if c[-1] == str(original_root/'build'/profile/binary)]
            require(len(links)==1,'unique profile link')
            equal([a for a in links[0] if a.endswith('/driver.o')],
                  [str(original_root/'build'/(common+'-objects')/'driver.o')])
            symbols = (root/'build'/profile/(binary+'-undefined.txt')).read_text()
            equal('_ZNSt6chrono3_V212steady_clock3nowEv' in symbols,binary=='plain')
    rows = inventory()
    equal([r['label'] for r in checks['records']],[r['label'] for r in rows])
    peaks, total, comparisons = [],0,0
    for record,row in zip(checks['records'],rows):
        folder,label = root/'checks',row['label']
        out,err = folder/(label+'.stdout'),folder/(label+'.stderr')
        equal(record['returncode'],row['code']); equal(record['fault'],row['fault'])
        equal(sha(out),record['stdout_sha256']); equal(sha(err),record['stderr_sha256'])
        args = list(row['args'])
        if row['parity']: args.append(str(Path(checks['root'])/'checks'/(label+'.parity')))
        equal(record['args'],[row['profile'],row['binary'],*args])
        equal(record['scope_command'],['systemd-run','--user','--scope','--expand-environment=no',
            '-p','MemoryMax=256M','-p','MemorySwapMax=0','bash',str(original_root/'build/tower_public_scope.sh'),
            'flock','-n','/tmp/leopard-gf8-authoritative.lock','timeout','--signal=TERM','--kill-after=5',
            '120','prlimit','--cpu=60:60','--core=0:0','--',
            str(original_root/'build'/row['profile']/row['binary']),*args])
        native,stderr,peak = split_scope(out.read_text(),err.read_text(),row['code'])
        equal(record['memory_peak'],peak)
        peaks.append(peak); verify_record(row,native,stderr)
        if row['parity']:
            cell = int(row['args'][1]); parity = folder/(label+'.parity')
            equal(parity.stat().st_size,CELLS[cell][1]*CELLS[cell][2])
            equal(sha(parity),record['parity_sha256'])
            baseline = folder/f'native-{cell}-NNNN-1-plain.parity'
            if parity != baseline:
                total += previous.compare(parity,baseline); comparisons += 1
    equal(sorted(p.name for p in (root/'checks').glob('*.stdout')),sorted(r['label']+'.stdout' for r in rows))
    # Fresh native executions must reproduce the earlier independently retained
    # native product outputs too. Do not count reference self-comparisons.
    prior = root/'prior-native'
    if not prior.is_dir():
        prior = Path('/home/catid/leopard/.research/leopard-79h/avx2-adjacent-public.NFkTeG')
    equal(sha(prior/'SHA256SUMS'),PRIOR_MANIFEST)
    pins = dict(line.split('  ',1)[::-1] for line in (prior/'SHA256SUMS').read_text().splitlines())
    prior_bytes = 0
    for cell in range(9):
        name = f'checks/native-{cell}-NNNN-1-plain.parity'
        equal(sha(prior/name),pins[name])
        prior_bytes += previous.compare(root/name,prior/name)
    resources = {phase:scope_footer((root/(phase+'.log')).read_text(),0,limit*1024**2)[1]
                 for phase,limit in (('build',512),('checks',256))}
    return dict(bead=BEAD,positive=sum(r['code']==0 for r in rows),
        clock_aborts=sum(r['code']==86 for r in rows),clock_faults=sum(bool(r['fault']) for r in rows),
        cli_refusals=sum('-bad-' in r['label'] for r in rows),
        poisoned_read_refusals=2,parity_comparisons=comparisons,parity_bytes=total,
        prior_native_bytes=prior_bytes,maximum_native_peak=max(peaks),scope_peaks=resources,
        structural_counts=[dict(zip(FIELDS,structural(c))) for c in range(9)],
        timed=False,default_enabled=False,performance_qualified=False)


if __name__ == '__main__': print(json.dumps(replay(Path(sys.argv[1]).resolve()),sort_keys=True))

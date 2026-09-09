#!/usr/bin/env python3
"""Collector-free, clock-free replay for leopard-79h.38.5.4.19.1.1."""
import hashlib
import json
from pathlib import Path
import sys

from verify_auto_r19932_checks import ARCHIVES, CELLS, INPUTS, OUTPUTS, public_expected, scope

BEAD = 'leopard-79h.38.5.4.19.1.1'
MASK = (1 << 64) - 1


def require(ok, message):
    if not ok:
        raise ValueError(message)


def equal(actual, expected):
    require(json.dumps(actual, sort_keys=True, allow_nan=False) ==
            json.dumps(expected, sort_keys=True, allow_nan=False), 'record mismatch')


def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def parse(text):
    def unique(pairs):
        result = {}
        for key, value in pairs:
            require(key not in result, 'duplicate JSON key')
            result[key] = value
        return result
    return json.loads(text, object_pairs_hook=unique)


def read(path):
    return parse(path.read_text())


def valid(profile, cell, schedule, group):
    require(profile in ARCHIVES and type(cell) is int and 0 <= cell < 9, 'profile/cell')
    require(schedule in (('NNNN',) if profile == 'native' else ('0110','1001','0000','1111')), 'schedule')
    require(type(group) is int and (group == 1 or (cell == 8 and group == 256)), 'group')


def expected(profile, cell, schedule, group, exercise):
    valid(profile, cell, schedule, group)
    require(type(exercise) is bool, 'exercise')
    native = profile == 'native'
    k,r,size = CELLS[cell]
    old = public_expected(profile, 0, cell)
    calls = 1 + (25 * group if exercise else 0)
    probes = [0 if native else int(2 <= cell <= 4 or (cell < 2 and state == '1')) for state in schedule]
    return dict(schema='leopard-paired-r19932/v1', codec=profile+':'+ARCHIVES[profile],
        cell=cell,k=k,r=r,bytes=size,api=old['api'],schedule=schedule,group=group,
        warmup_passes=4 if exercise else 0,exercise_passes=21 if exercise else 0,
        encode_calls=4*calls,selections=104 if exercise else 4,per_slot_calls=[calls]*4,
        probes=probes,scratch_bytes=old['scratch_bytes'],input_hash=INPUTS[cell],
        output_hash=OUTPUTS[cell],timed=False,default_enabled=False)


def witness(profile, cell, schedule, group, kind):
    valid(profile, cell, schedule, group)
    require(kind in ('check','exercise','clock-guard'), 'kind')
    api = 2 if profile == 'native' else 1 if cell == 1 else 0
    states = [2 if s == 'N' else int(s) for s in schedule]
    sequence = list(states)
    for _ in range(0 if kind == 'check' else 4 if kind == 'clock-guard' else 25):
        for state in states:
            sequence.extend([state] * group)
    counts, apis, digest = [0]*3, [0]*3, 14695981039346656037
    for state in sequence:
        counts[state] += 1
        apis[api] += 1
        digest = ((digest ^ (state + 4*api)) * 1099511628211) & MASK
    return dict(schema='leopard-paired-witness/v1',calls=len(sequence),states=counts,
                apis=apis,order_hash=f'{digest:016x}')


def cases():
    rows = []
    for profile in ARCHIVES:
        for cell in range(9):
            for schedule in (('NNNN',) if profile=='native' else ('0110','1001','0000','1111')):
                for group in ((1,256) if cell==8 else (1,)):
                    for kind, binary in (('check','plain'),('exercise','witness')):
                        label = f'{profile}-{cell}-{schedule}-{group}-{kind}'
                        rows.append((label,profile,cell,schedule,group,kind,binary,0))
        for cell,group in ((0,1),(8,256)):
            for schedule in (('NNNN',) if profile=='native' else ('0110','1001','0000','1111')):
                label = f'{profile}-{cell}-{schedule}-{group}-clock-guard'
                rows.append((label,profile,cell,schedule,group,'clock-guard','witness',86))
    return rows


def bad_arguments(profile):
    schedule = 'NNNN' if profile=='native' else '0110'
    base = ['--check','0',schedule,'1']
    return [[],['--check'],['--measure','0',schedule,'1'],base+['a','b'],
            ['--check','9',schedule,'1'],['--check','00',schedule,'1'],
            ['--check','0',schedule,'0'],['--check','0',schedule,'256'],
            ['--check','8',schedule,'0256'],['--check','0','0101','1'],
            ['--check','0','111','1'],['--check','0',('0110' if profile=='native' else 'NNNN'),'1'],
            ['--clock-guard','0',schedule,'1','forbidden-parity']]


def compare(left, right):
    total = 0
    with left.open('rb') as a, right.open('rb') as b:
        while True:
            chunk = a.read(65536)
            require(chunk == b.read(65536), 'full native parity differs')
            if not chunk:
                return total
            total += len(chunk)


def guards(root):
    folder = root/'guards'
    build, checks = read(folder/'build.json'), read(folder/'checks.json')
    for state in (build,checks):
        equal([state['bead'],state['completed'],state['timed']], [BEAD,True,False])
    for name,digest in build['artifacts'].items():
        equal(sha(folder/name),digest)
    wanted = [('release-canary',0),('sanitize-canary',0),('sanitize-underflow',1),('sanitize-overflow',1)]
    equal([(r['label'],r['returncode']) for r in checks['records']],wanted)
    for record,(label,code) in zip(checks['records'],wanted):
        stdout,stderr = folder/(label+'.stdout'),folder/(label+'.stderr')
        equal(sha(stdout),record['stdout_sha256']); equal(sha(stderr),record['stderr_sha256'])
        if code:
            equal(stdout.read_text(),'')
            require('ERROR: AddressSanitizer: use-after-poison' in stderr.read_text(), 'ASan poison detection')
            require('READ of size 1' in stderr.read_text(), 'ASan read boundary detection')
        else:
            equal(stderr.read_text(),'')
            equal(stdout.read_text(),'both canaries rejected; zero-size and restored buffers pass\n')
    return dict(canary_checks=2,expected_asan_rejections=2)


def replay(root):
    build = read(root/'build/build.json')
    equal([build['bead'],build['completed'],build['timed']], [BEAD,True,False])
    for name,digest in build['inputs'].items():
        equal(sha(Path(name)),digest)
    for name,digest in build['artifacts'].items():
        equal(sha(root/'build'/name),digest)
    for profile,digest in ARCHIVES.items():
        equal(sha(root/'build'/profile/'codec.a'),digest)
    state = read(root/'checks/checks.json')
    equal([state['bead'],state['completed'],state['timed']], [BEAD,True,False])
    rows = cases()
    equal([row['label'] for row in state['records']],
          [row[0] for row in rows]+[f'{p}-bad-{i}' for p in ARCHIVES for i in range(len(bad_arguments(p)))])
    parity_bytes, parity_files = 0, 0
    for record, (label,profile,cell,schedule,group,kind,binary,code) in zip(state['records'],rows):
        args = ['--'+kind,str(cell),schedule,str(group)]
        parity = root/'checks'/(label+'.parity')
        if kind=='exercise':
            args.append(str(Path(state['root'])/'checks'/parity.name))
        equal(record['args'],[profile,binary,*args])
        equal(record['returncode'],code)
        stdout = root/'checks'/(label+'.stdout')
        stderr = root/'checks'/(label+'.stderr')
        equal(sha(stdout),record['stdout_sha256']); equal(sha(stderr),record['stderr_sha256'])
        equal(stderr.read_text(), 'unexpected driver benchmark clock\n' if code==86 else '')
        lines = stdout.read_text().splitlines()
        wanted = [] if code else [expected(profile,cell,schedule,group,kind=='exercise')]
        if binary=='witness': wanted.append(witness(profile,cell,schedule,group,kind))
        equal([parse(line) for line in lines],wanted)
        if kind=='exercise':
            equal(parity.stat().st_size,CELLS[cell][1]*CELLS[cell][2])
            equal(sha(parity),record['parity_sha256'])
            baseline = root/'checks'/f'native-{cell}-NNNN-{group}-exercise.parity'
            if profile!='native':
                parity_bytes += compare(parity,baseline); parity_files += 1
    for record in state['records'][len(rows):]:
        profile,_,index = record['label'].split('-')
        equal(record['args'],[profile,'plain',*bad_arguments(profile)[int(index)]])
        equal(record['returncode'],1)
        stdout = root/'checks'/(record['label']+'.stdout')
        stderr = root/'checks'/(record['label']+'.stderr')
        equal(sha(stdout),record['stdout_sha256']); equal(sha(stderr),record['stderr_sha256'])
        equal(stdout.read_text(),''); require(bool(stderr.read_text()),'missing refusal diagnostic')
    equal(sorted(p.name for p in (root/'checks').glob('*.stdout')),
          sorted(r['label']+'.stdout' for r in state['records']))
    # Retain the actual max-event count: checks exited 0 with no OOM or swap,
    # but reached the 256 MiB ceiling. This is deliberately NOT an all-zero claim.
    resources = {phase:scope((root/(phase+'.log')).read_text(),
                            (512 if phase.endswith('build') else 256)*1024**2,
                            max_events=1519 if phase=='checks' else 0)
                 for phase in ('build','checks','guard-build','guard-checks')}
    return dict(bead=BEAD, timed=False, default_enabled=False, positive=sum(r[-1]==0 for r in rows),
                clock_aborts=sum(r[-1]==86 for r in rows),cli_refusals=len(state['records'])-len(rows),
                parity_files=parity_files,parity_bytes=parity_bytes,guards=guards(root),resources=resources)


if __name__=='__main__':
    print(json.dumps(replay(Path(sys.argv[1]).resolve()),sort_keys=True))

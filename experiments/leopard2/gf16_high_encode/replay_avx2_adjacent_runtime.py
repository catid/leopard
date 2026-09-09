#!/usr/bin/env python3
"""Collector-free focused-runtime replay. Full parity/timing gates remain open."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys
from audit_avx2_adjacent_runtime import audit
from audit_avx2_pair_schedule import resource
from avx2_adjacent_runtime import overlay
from verify_paired_r19932 import equal,parse,require,sha
from verify_auto_r19932_checks import ARCHIVES,guard_expected

BEAD = 'leopard-79h.38.5.4.18.3.1'
STATIC = Path('/home/catid/leopard/.research/leopard-79h/avx2-adjacent-qualified.wxkx9v2b')
SELECTORS = ['adjacent','forward-ranges','pairs','split',*map(str,range(8)),'roundtrip','concurrent']
BAD = [[],['off'],['bad','--pairs'],['on','--pairs','extra'],['off','8'],['on','--measure']]


def focused(text,profile,mode,selector):
    require(profile in ('release','trace','sanitize') and mode in ('off','on') and selector in SELECTORS,
            'focused selector identity')
    lines = text.splitlines(); require(bool(lines),'missing focused output')
    record = parse(lines.pop())
    equal(sorted(record),['blocks','calls','mode','schema','timed','trace'])
    enabled = int(mode=='on'); traced = profile!='release'
    equal([record['schema'],record['mode'],record['trace'],record['timed']],
          ['adjacent-runtime-focused/v1',enabled,traced,False])
    for key in ('calls','blocks'):
        require(type(record[key]) is list and len(record[key])==2,'family dimension')
        for family in record[key]:
            require(type(family) is list and len(family)==2 and
                    all(type(v) is int and 0<=v<2**64 for v in family),'state/counter dimension')
            equal(family[1-enabled],0)
            if not traced: equal(family,[0,0])
    if traced:
        require(all(family[enabled]>0 for family in record['calls']),'missing traced families')
        if selector=='adjacent':
            require(record['calls'][0][enabled]>=65535 and record['calls'][1][enabled]>=65535+918,
                    'missing exercised pair calls')
    if selector in tuple(map(str,range(8))):
        require(len(lines)==1,'guard output inventory')
        equal(parse(lines[0]),guard_expected(True,int(selector)))
    elif selector=='concurrent':
        equal([len(lines),lines.count('roundtrip cell 1 passed'),lines.count('roundtrip cell 2 passed'),
               lines.count('four-thread both-field roundtrips passed')],[18,9,8,1])
    else:
        wanted = {
            'adjacent':['adjacent pairs: forward=65535 accumulating=65535 boundary_accumulations=918'],
            'forward-ranges':['forward range cases: 384'],'pairs':['pair kernel cases: 66147'],
            'split':['split range cases: 256'],
            'roundtrip':['roundtrip cell '+str(i)+' passed' for i in range(3)]}
        equal(lines,wanted[selector])
    return {key:[family[enabled] for family in record[key]] for key in ('calls','blocks')}


def replay(root):
    build = json.loads((root/'build/build.json').read_text())
    checks = json.loads((root/'checks/checks.json').read_text())
    for state in (build,checks): equal([state['bead'],state['completed'],state['timed']],[BEAD,True,False])
    equal(sorted(build['profiles']),['release','sanitize','trace'])
    oldroot = Path(next(iter(build['inputs']))).parents[2]
    def locate(name):
        path = Path(name)
        return root/path.relative_to(oldroot) if path.is_relative_to(oldroot) else path
    for pins in (build['inputs'],build['external_inputs'],checks['pins']):
        for name,digest in pins.items(): equal(sha(locate(name)),digest)
    source,drivers = root/'build/source',root/'build/drivers'
    equal((source/'Leopard2BackendAVX2.cpp').read_text(),overlay((STATIC/'source/Leopard2BackendAVX2.cpp').read_text()))
    original = (drivers/'test_avx2_adjacent_schedule.cpp').read_text()
    equal((drivers/'adjacent_focused.inc').read_text(),original.replace(
        'int main(int argc, char** argv)','int LeoAdjacentFocusedMain(int argc, char** argv)',1))
    equal(sha(drivers/'test_avx2_adjacent_schedule.cpp'),sha(STATIC/'checks-build/drivers/test_avx2_adjacent_schedule.cpp'))
    for profile,entry in build['profiles'].items():
        directory = root/'build'/profile; archive = directory/'candidate.a'; original = Path(entry['original'])
        equal(entry['original_sha256'],ARCHIVES['sanitize' if profile=='sanitize' else 'release'])
        equal(sha(original),entry['original_sha256']); equal(sha(archive),entry['archive_sha256'])
        members = subprocess.check_output(['ar','t',str(original)],text=True).splitlines()
        require(len(members)==len(set(members))==24,'original member count')
        equal(subprocess.check_output(['ar','t',str(archive)],text=True).splitlines(),members+['avx2_adjacent_control.cpp.o'])
        equal(sorted(entry['unchanged_members']),sorted(set(members)-{'Leopard2BackendAVX2.cpp.o'}))
        for member in members+['avx2_adjacent_control.cpp.o']:
            actual = subprocess.check_output(['ar','p',str(archive),member])
            if member=='Leopard2BackendAVX2.cpp.o':
                equal(hashlib.sha256(actual).hexdigest(),entry['object_sha256'])
                equal(hashlib.sha256(actual).hexdigest(),sha(directory/member))
            elif member=='avx2_adjacent_control.cpp.o': equal(hashlib.sha256(actual).hexdigest(),entry['files'][member])
            else:
                require(actual==subprocess.check_output(['ar','p',str(original),member]),'unrelated archive member')
                equal(hashlib.sha256(actual).hexdigest(),entry['unchanged_members'][member])
        for name,digest in entry['files'].items(): equal(sha(directory/name),digest)
    inventory = []
    for profile in ('release','trace','sanitize'):
        inventory.append((profile+'-unit',profile,'control-unit',[],0,None,None))
        for mode in ('off','on'):
            for selector in SELECTORS:
                arg = selector if selector in tuple(map(str,range(8))) else '--'+selector
                inventory.append((profile+'-'+mode+'-'+selector,profile,'focused',[mode,arg],0,mode,selector))
        for i,args in enumerate(BAD): inventory.append((profile+'-bad-'+str(i),profile,'focused',args,1,None,None))
    equal([r['label'] for r in checks['records']],[x[0] for x in inventory])
    observations = {}
    for row,(label,profile,binary,args,code,mode,selector) in zip(checks['records'],inventory):
        equal([row['returncode'],row['expected']],[code,code])
        equal(row['argv'],['/usr/bin/prlimit','--cpu=30:30','--fsize=1048576:1048576','--',
                           str(oldroot/'build'/profile/binary),*args])
        for suffix in ('stdout','stderr'): equal(sha(root/'checks'/(label+'.'+suffix)),row[suffix+'_sha256'])
        text = (root/'checks'/(label+'.stdout')).read_text()
        error = (root/'checks'/(label+'.stderr')).read_text()
        if code==0:
            equal(error,'')
            if binary=='control-unit':
                equal(parse(text),dict(schema='adjacent-control-unit/v1',trace=profile!='release',timed=False))
            else: observations[label] = focused(text,profile,mode,selector)
        else: require(bool(error),'CLI refusal missing')
    equal(sorted(p.name for p in (root/'checks').iterdir()),sorted(['checks.json']+
          [x[0]+'.'+s for x in inventory for s in ('stdout','stderr')]))
    for selector in SELECTORS:
        expected = observations['trace-off-'+selector]
        for profile,mode in (('trace','on'),('sanitize','off'),('sanitize','on')):
            equal(observations[profile+'-'+mode+'-'+selector],expected)
    codegen = audit(root)
    equal(codegen,json.loads((root/'codegen.json').read_text()))
    return dict(bead=BEAD,positive=87,cli_refusals=18,timed=False,production_changed=False,
        focused_matrix_pass=True,full_native_parity_qualified=False,clock_boundary_qualified=False,
        performance_qualified=False,full_runtime_qualification=False,
        trace_counts_include_initialization=True,worker_trace_counts_aggregated=False,
        trace_counts={s:observations['trace-off-'+s] for s in SELECTORS},codegen=codegen,
        resources={name:resource(root/(name+'.log'),limit) for name,limit in
                   (('build',536870912),('checks',268435456),('codegen',268435456))})


if __name__=='__main__': print(json.dumps(replay(Path(sys.argv[1]).resolve()),indent=2,sort_keys=True))

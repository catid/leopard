#!/usr/bin/env python3
"""Read-only full public qualification replay; no codec or benchmark execution."""
import json
from pathlib import Path
import sys

from verify_avx2_adjacent_public import replay, BEAD, read, equal, require, sha, compare, scope


def qualification(root):
    result = replay(root)
    folder = root/'guards'
    build,checks = read(folder/'build.json'),read(folder/'checks.json')
    for state in (build,checks): equal([state['bead'],state['completed'],state['timed']],[BEAD,True,False])
    for name,digest in build['artifacts'].items(): equal(sha(folder/name),digest)
    equal(sha(folder/'test_avx2_adjacent_public_guards.cpp'),
          sha(Path(__file__).with_name('test_avx2_adjacent_public_guards.cpp')))
    public_build = read(root/'build/build.json')
    original_root = Path(read(root/'checks/checks.json')['root'])
    equal(len(build['commands']),2)
    for profile,command in zip(('release','sanitize'),build['commands']):
        common = 'l2-release-objects' if profile=='release' else 'sanitize-objects'
        recipe, = [c for c in public_build['commands'] if '-c' in c and
                   c[-1]==str(original_root/'build'/common/'driver.o')]
        flags = recipe[recipe.index('c++')+1:recipe.index('-c')]
        equal(command,['prlimit','--cpu=120:120','--','c++',*flags,'-I'+str(original_root/'build'),
            str(original_root/'guards/test_avx2_adjacent_public_guards.cpp'),
            str(original_root/'build'/profile/'link.o'),str(original_root/'build'/common/'abort-clock.o'),
            str(original_root/'build'/profile/'codec.a'),str(original_root/'build/clock_guard.cpp'),
            '-Wl,--wrap=_ZNSt6chrono3_V212steady_clock3nowEv','-o',str(original_root/'guards'/profile)])
    wanted = [('release-canary',0),('sanitize-canary',0),('sanitize-underflow',1),('sanitize-overflow',1)]
    equal([(r['label'],r['returncode']) for r in checks['records']],wanted)
    for row,(label,code) in zip(checks['records'],wanted):
        profile,case = label.split('-')
        equal(row['args'],['prlimit','--cpu=10:10','--',str(original_root/'guards'/profile),case])
        out,err = folder/(label+'.stdout'),folder/(label+'.stderr')
        equal(sha(out),row['stdout_sha256']); equal(sha(err),row['stderr_sha256'])
        if code:
            equal(out.read_text(),'')
            require('ERROR: AddressSanitizer: use-after-poison' in err.read_text() and
                    'READ of size 1' in err.read_text(), 'ASan poisoned boundary rejection')
        else:
            equal(err.read_text(),'')
            equal(out.read_text(),'both canaries rejected; zero-size and restored buffers pass\n')
    result['guard_checks'] = dict(canaries=2,expected_asan_rejections=2)
    for phase,limit in (('guard-build',512*1024**2),('guard-checks',256*1024**2)):
        result['resources'][phase] = scope((root/(phase+'.log')).read_text(),limit,max_events=0)
    # New native executions must also reproduce the original retained product
    # parity, not just each other's output. Batch cell8 deliberately maps to0.
    reference = Path('/home/catid/leopard/.research/leopard-79h/avx2-pair-screen.m_vukpir/preparation')
    result['prior_native_bytes'] = 0
    result['prior_native_sha256'] = []
    for cell in range(9):
        path = reference/f'native-native-{0 if cell==8 else cell}.parity'
        actual = root/'checks'/f'native-{cell}-NNNN-1-plain.parity'
        result['prior_native_bytes'] += compare(actual,path)
        result['prior_native_sha256'].append(sha(path))
    # Carry the already-passed focused matrix forward by exact record identity;
    # do not rerun that matrix or infer that the new frontend repeats its tests.
    focused = Path('/home/catid/leopard/.research/leopard-79h/avx2-adjacent-runtime-focused.yMDquC')
    pins = {'build/build.json':'5c81b185be5dfa33c7677a2536b35df6e9c5df03483ba1cc2cbeb47573809e17',
            'checks/checks.json':'764de17680b7b050b4a48b58b75df7535ad80ffc6e820e0518e6fbb319041731',
            'replay.json':'8bee89cfc31835576230e96fde9e8d0d899394a4dfaafdf9b00cb9302944771c'}
    for name,digest in pins.items(): equal(sha(focused/name),digest)
    result['prior_focused_record_pins'] = pins
    result['full_runtime_qualification'] = True
    result['performance_qualified'] = False
    return result


if __name__=='__main__': print(json.dumps(qualification(Path(sys.argv[1]).resolve()),sort_keys=True))

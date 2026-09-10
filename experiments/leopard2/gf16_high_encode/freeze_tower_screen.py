#!/usr/bin/env python3
"""Freeze qualified inputs and run only --check; never starts a timing attempt."""
import json
from pathlib import Path
import shutil
import subprocess
import sys

from run_tower_screen import (PLAN,PROTOCOL,PINS,SOURCE_FILES,read,verify,validate,validate_plan)
from verify_paired_r19932 import equal,require,sha
from replay_tower_screen import inputs, record

ROOT = Path('/tmp/leopard-tower-screen.0mjTlC')
QUALIFIED = Path('/home/catid/leopard/.research/leopard-79h/tower-public-qualified.KPkPef')
OLD = Path('/tmp/leopard-paired-integration.la42Cz/frozen')


def freeze():
    require(not (ROOT/'attempt1').exists(),'attempt already exists')
    source = Path(__file__).resolve().parent
    validate_plan(read(source/PLAN))
    frozen = ROOT/'frozen'; frozen.mkdir(mode=0o700)
    paths = {name:source/name for name in SOURCE_FILES|{PLAN}}
    paths.update({
        'native':QUALIFIED/'build/native/plain', 'original':QUALIFIED/'build/original/plain',
        'current':QUALIFIED/'build/release/plain',
        'native.a':QUALIFIED/'build/native/codec.a', 'original.a':QUALIFIED/'build/original/codec.a',
        'current.a':QUALIFIED/'build/release/codec.a',
        'l2-driver.o':QUALIFIED/'build/l2-release-objects/driver.o',
        'native-driver.o':QUALIFIED/'build/native-objects/driver.o',
        'build.json':QUALIFIED/'build/build.json', 'checks.json':QUALIFIED/'checks/checks.json',
        'qualification.log':QUALIFIED/'replay-normal.log','qualification.manifest':QUALIFIED/'SHA256SUMS',
        'leopard2.cpp':OLD/'leopard2.cpp','Leopard2Direct.h':OLD/'Leopard2Direct.h',
        'check.sh':OLD/'check.sh','check-condition.sh':OLD/'check-condition.sh'})
    for name,path in paths.items():
        require(path.is_file() and not path.is_symlink(),'input file')
        if name in PINS: equal(sha(path),PINS[name])
        # Actual private copies, not hard links to another lane's executables.
        shutil.copyfile(path,frozen/name)
        equal(sha(path),sha(frozen/name))
        (frozen/name).chmod(0o555 if name in ('native','original','current') else 0o444)
    pins = dict(schema='leopard-tower-pins/v1',files={n:sha(frozen/n) for n in sorted(paths)})
    with (frozen/'pins.json').open('x') as stream: json.dump(pins,stream,indent=2); stream.write('\n')
    (frozen/'pins.json').chmod(0o444); frozen.chmod(0o555)
    verify(frozen,pins); inputs(ROOT)
    checks = ROOT/'preflight-qualification'; checks.mkdir(mode=0o700)
    state = dict(bead=PROTOCOL['bead'],timed=False,completed=False,records=[])
    env = dict(PATH='/usr/bin:/bin',LANG='C',LC_ALL='C',OMP_NUM_THREADS='1',OMP_DYNAMIC='FALSE',OMP_THREAD_LIMIT='1')
    try:
        for cell in range(9):
            for order in ('NNNN','PPPP','0000','1111'):
                name = f'check-{cell}-{order}'
                argv = ['/usr/bin/taskset','-c','26','/usr/bin/prlimit','--cpu=30:30','--core=0:0',
                    '--fsize=1048576:1048576','--',str(frozen/('native' if order=='NNNN' else 'original' if order=='PPPP' else 'current')),
                    '--check',str(cell),order,str(PROTOCOL['groups'][cell])]
                outpath,errpath = checks/(name+'.stdout'),checks/(name+'.stderr')
                with outpath.open('xb') as out,errpath.open('xb') as err:
                    child = subprocess.run(argv,stdout=out,stderr=err,env=env,timeout=60)
                require(child.returncode==0 and errpath.stat().st_size==0,'check failed: '+name)
                row = read(outpath); validate(row,cell,order,False)
                record(row,cell,order,False)
                state['records'].append(dict(name=name,argv=argv,returncode=child.returncode,record=row))
            print(f'cell {cell}: native/original/OFF/ON four-call preflights validated without clocks',flush=True)
        verify(frozen,pins); inputs(ROOT)
        state['completed'] = True
    finally:
        with (checks/'checks.json').open('x') as stream: json.dump(state,stream,indent=2); stream.write('\n')


if __name__=='__main__':
    require(len(sys.argv)==1,'no arguments; one fixed freeze path')
    freeze()

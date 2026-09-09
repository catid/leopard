#!/usr/bin/env python3
"""Untimed audit of actual runtime object loops; not a performance claim."""
import json
from pathlib import Path
import re
import subprocess
import sys
from audit_avx2_adjacent_schedule import FUNCTIONS,isa_ceiling
from audit_avx2_pair_schedule import instructions,require,sha


def loops(text,function):
    parts = re.split(r'^([0-9a-f]+) <(.+)>:\n',text,flags=re.M)
    bodies = [parts[i+2] for i in range(1,len(parts),3) if parts[i+1]==function]
    require(len(bodies)==1,'missing/duplicate function')
    ops = instructions(bodies[0]); candidates = []
    for end,_,assembly in ops:
        branch = re.match(r'j(?!mp\b)[a-z]+\s+([0-9a-f]+)\s',assembly)
        if not branch or int(branch[1],16)>=end: continue
        start = int(branch[1],16); body = [o for o in ops if start<=o[0]<=end]
        if sum(o[2].split()[0]=='vpshufb' for o in body)==8:
            require(body[0][0]==start,'invalid branch target')
            candidates.append((start,end,body))
    inner = [(a,b,body) for a,b,body in candidates if not any(
        a<=c and d<=b and (a,b)!=(c,d) for c,d,_ in candidates)]
    require(bool(inner),'missing eight-shuffle loop')
    result = []
    for start,end,body in sorted(inner):
        require(not any('leo_adjacent_' in o[2] or re.match(r'call\w*\b',o[2]) for o in body),
                'runtime control/observer call in vector loop')
        result.append(dict(start=hex(start),end=hex(end),instructions=len(body),byte_shuffles=8,
            stack_references=[o[2] for o in body if re.search(r'\([^)]*%(?:rsp|rbp)\b',o[2])]))
    return result


def audit(root):
    build = json.loads((root/'build/build.json').read_text())
    require(build['completed'] is True and build['timed'] is False,'complete untimed build')
    directory = root/'build/release'; obj = directory/'Leopard2BackendAVX2.cpp.o'
    require(sha(obj)==build['profiles']['release']['object_sha256'],'object drift')
    actual = subprocess.check_output(['objdump','-drwC',str(obj)],text=True)
    require(instructions(actual)==instructions((directory/'disassembly.txt').read_text()),'disassembly drift')
    isa_ceiling(actual)
    result = {key:loops(actual,name) for key,name in FUNCTIONS.items()}
    for key,values in result.items(): require(len(values)==(1 if key=='inverse' else 2),'runtime loop inventory')
    return dict(bead=build['bead'],timed=False,codec_executed=False,
                object_sha256=sha(obj),loops=result,isa_ceiling_pass=True,
                runtime_control_outside_vector_loops=True)


if __name__=='__main__': print(json.dumps(audit(Path(sys.argv[1]).resolve()),indent=2,sort_keys=True))

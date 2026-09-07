#!/usr/bin/env python3
"""Independent stdlib replay of the fixed .16 attempt; never executes a codec."""
import hashlib
import json
import math
from pathlib import Path
import statistics
import sys

PLAN_SHA = "029b62d7370664f56f950846b6061302d6e6f8cb46127d5f5fc72ee43f0b6f81"
PINS_SHA = "de90ed1c119b0fd58613aa61930ce374b961e76044e6990f365989c9c964add8"
ORDERS = ((0,1,2,3,3,2,1,0), (1,2,3,0,0,3,2,1), (2,3,0,1,1,0,3,2))


def require(value, message):
    if not value: raise ValueError(message)


def sha(path):
    with path.open("rb") as stream: return hashlib.file_digest(stream,"sha256").hexdigest()


def read(path):
    require(path.stat().st_size <= 1048576, "bounded record")
    return json.loads(path.read_text())


def equal(actual, expected):
    require(json.dumps(actual,sort_keys=True) == json.dumps(expected,sort_keys=True), "record mismatch")


def close(actual, expected):
    if type(expected) is float:
        require(type(actual) is float and math.isfinite(actual) and
                math.isclose(actual,expected,rel_tol=1e-12,abs_tol=1e-12), "derived floating value")
    elif type(expected) is dict:
        require(type(actual) is dict and set(actual)==set(expected), "derived object schema")
        for key in expected: close(actual[key],expected[key])
    elif type(expected) is list:
        require(type(actual) is list and len(actual)==len(expected), "derived list schema")
        for a,b in zip(actual,expected): close(a,b)
    else:
        equal(actual,expected)


def checked_pair(root,label,cell,mode,measured,expected):
    record,trace=read(root/(label+".stdout")),read(root/(label+".stderr"))
    equal({key:value for key,value in record.items() if key!="samples_ns"},expected[cell])
    values=record.get("samples_ns")
    require(type(values) is list and len(values)==(21 if measured else 0) and
            all(type(value) is int and value>0 for value in values), "raw samples")
    encodes=26 if measured else 1
    matches=encodes*2 if cell==0 else 0
    equal(trace,dict(schema="gfni-combined-timing/v1",cell=cell,mode=mode,encodes=encodes,
        calls=encodes*(2,2,1,1,2,1)[cell],matches=matches,
        first=matches if mode in (1,3) else 0,terminal=matches if mode in (2,3) else 0,
        timed=measured,exercise=False))
    return record,trace


def resource_peak(path,complete):
    lines=path.read_text().splitlines()
    require(lines.count("memory.peak")==1,"resource block")
    index=lines.index("memory.peak"); peak=int(lines[index+1])
    require(0<peak<=268435456 and lines[index+2:]==["memory.max","268435456","memory.events",
        "low 0","high 0","max 0","oom 0","oom_kill 0","oom_group_kill 0",
        "memory.swap.current","0","memory.swap.max","0"],"resource envelope")
    equal([line for line in lines if line.startswith("\tExit status:")],
          ["\tExit status: %d" % (0 if complete else 1)])
    return peak


def derive(medians):
    cells=[]
    for cell in range(6):
        rounds={comparison:{str(mode):[] for mode in (1,2,3)} for comparison in ("factorial","same_off")}
        interaction=[]
        for round_id,order in enumerate(ORDERS):
            for comparison in ("factorial","same_off"):
                logs=[]
                for mode in range(4):
                    slots=[slot for slot,value in enumerate(order) if value==mode]
                    logs.append(sum(math.log(medians[cell,round_id,comparison,slot]) for slot in slots)/2)
                for mode in (1,2,3): rounds[comparison][str(mode)].append(logs[0]-logs[mode])
                if comparison=="factorial": interaction.append(logs[1]+logs[2]-logs[0]-logs[3])
        cells.append(dict(cell=cell,
            round_ratios={comparison:{mode:[math.exp(v) for v in values] for mode,values in modes.items()}
                          for comparison,modes in rounds.items()},
            ratios={comparison:{mode:math.exp(sum(values)/3) for mode,values in modes.items()}
                    for comparison,modes in rounds.items()},
            interaction_rounds=[math.exp(v) for v in interaction],interaction_factor=math.exp(sum(interaction)/3)))
    controls=[v for cell in cells for v in cell["ratios"]["same_off"].values()]
    controls += [v for cell in cells[1:] for v in cell["ratios"]["factorial"].values()]
    valid=all(1/1.02<=v<=1.02 for v in controls)
    target=cells[0]["ratios"]["factorial"]["3"]>=1.05 and min(cells[0]["round_ratios"]["factorial"]["3"])>1
    decision="inconclusive_controls" if not valid else "continue_to_future_qualification" if target else "reject_for_this_screen"
    return dict(decision=decision,cells=cells,aggregate_controls=33,confidence_intervals=False,
                production_promotion=False,exact_leopard1_claim=False,authoritative_v19=False)


def replay(root):
    frozen,attempt=root/"frozen",root/"attempt"
    require(sha(frozen/"gfni_combined_screen_plan.json")==PLAN_SHA,"preregistered plan")
    require(sha(frozen/"pins.json")==PINS_SHA,"frozen inventory")
    pins=read(frozen/"pins.json")
    require(len(pins["files"])==16,"pin count")
    for name,expected_sha in pins["files"].items():
        path=frozen/name
        require(Path(name).name==name and path.is_file() and not path.is_symlink() and
                not path.stat().st_mode & 0o222 and sha(path)==expected_sha,"frozen artifact")
    plan=read(frozen/"gfni_combined_screen_plan.json")
    equal(pins["files"],dict(plan["artifact_sha256"],**{"gfni_combined_screen_plan.json":PLAN_SHA}))
    expected=read(frozen/"expected.json"); state=read(attempt/"attempt.json")
    identity=state["executable_identity"]
    require(type(identity) is list and len(identity)==4 and all(type(v) is int and v>=0 for v in identity) and
            identity[1]>0 and identity[2]==(frozen/"timing").stat().st_size and not identity[3]&0o222,
            "recorded stable executable identity")
    equal(state["pins"],pins); equal(state["host"],plan["host"])
    require(state["schema"]=="gfni-combined-screen-attempt/v1" and state["plan_sha256"]==PLAN_SHA,"attempt identity")
    require(len(state["preflight"])==24,"preflight count")
    for cell in range(6):
        for mode in range(4):
            record,trace=checked_pair(attempt,"check-%d-%d"%(cell,mode),cell,mode,False,expected)
            equal(state["preflight"][cell*4+mode],dict(record=record,trace=trace))
    passive=state["passive"]
    require(type(passive["elapsed_ns"]) is int and passive["elapsed_ns"]>=10000000000,"passive duration")
    rows=state["invocations"]
    require(type(rows) is list and len(rows)<=288 and type(state["complete"]) is bool,"attempt bounds")
    cursor=0; medians={}
    for cell in range(6):
        for round_id,order in enumerate(ORDERS):
            for comparison in ("factorial","same_off"):
                for slot,variant in enumerate(order):
                    if cursor>=len(rows): continue
                    row=rows[cursor]; mode=variant if comparison=="factorial" else 0
                    equal([row[key] for key in ("cell","round","comparison","slot","variant","mode")],
                          [cell,round_id,comparison,slot,variant,mode])
                    require(type(row["sibling_delta"]) is int and row["sibling_delta"]>=0 and
                            (row["sibling_delta"]==0 or cursor==len(rows)-1),"sibling history")
                    label="cell-%d-round-%d-%s-slot-%d"%(cell,round_id,comparison,slot)
                    record,trace=checked_pair(attempt,label,cell,mode,True,expected)
                    equal(row["record"],record); equal(row["trace"],trace)
                    medians[cell,round_id,comparison,slot]=statistics.median(record["samples_ns"])
                    cursor+=1
    peak=resource_peak(root/"attempt.scope.log",state["complete"])
    if not state["complete"]:
        require("analysis" not in state,"partial attempt analyzed")
        passive_failure=passive["after"]>passive["before"] and not rows
        timed_failure=passive["after"]==passive["before"] and rows and rows[-1]["sibling_delta"]>0
        require(passive_failure or timed_failure,"unsupported failure requires separate investigation")
        equal(state["failure"],"ValueError: %ssibling activity; attempt stopped"%("passive " if passive_failure else ""))
        return dict(complete=False,timed_invocations=len(rows),performance_inference=False,
                    failure=state["failure"],peak_bytes=peak)
    require(len(rows)==288 and passive["after"]==passive["before"] and
            all(row["sibling_delta"]==0 for row in rows) and "failure" not in state,"complete isolation")
    analysis=derive(medians); close(state["analysis"],analysis)
    return dict(complete=True,timed_invocations=288,preflights=24,decision=analysis["decision"],
        cells=analysis["cells"],aggregate_controls=33,peak_bytes=peak,
        production_promotion=False,exact_leopard1_claim=False,
        performance_inference=analysis["decision"]!="inconclusive_controls")


if __name__=="__main__":
    require(len(sys.argv)==2,"usage: replay_gfni_combined_screen.py ROOT")
    print(json.dumps(replay(Path(sys.argv[1])),sort_keys=True))

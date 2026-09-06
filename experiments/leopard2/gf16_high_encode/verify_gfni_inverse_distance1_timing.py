#!/usr/bin/env python3
"""Retained-only timing-driver qualification: no benchmark execution."""
import hashlib
import json
from pathlib import Path
import sys
from verify_gfni_source_stage import equal, require, resource_peak


def digest(path):
    with path.open("rb") as stream: return hashlib.file_digest(stream,"sha256").hexdigest()


def replay(root, baseline, original):
    peaks=[]
    total=0
    require(digest(root / "bin/release-default") == digest(original / "bin/release") ==
            "64d0767a3e611c6cf3b7fa20856ac44d10800449ca7272a67a9833cf73e438d9", "default16 binary unchanged")
    for flavor in ("release","sanitize"):
        for capacity in (16,64):
            path=root / "checks" / ("%s-capacity-%d.log" % (flavor,capacity))
            peaks.append(resource_peak(path))
            line="distance-one capacity %d: both modes, overflow, exact records, reset, neighbor passed" % capacity
            require(path.read_text().splitlines().count(line)==1,"capacity proof")
        for mode in (0,1):
            for cell in range(6):
                for operation in (("check","exercise") if cell==0 else ("check",)):
                    path=root / "checks" / ("%s-%d-%d-%s" % (flavor,mode,cell,operation))
                    peaks.append(resource_peak(path.with_suffix(".log")))
                    encodes=26 if operation=="exercise" else 1
                    expected=(baseline / ("%d-current.json" % cell)).read_bytes()
                    require(path.with_suffix(".jsonl").read_bytes()==expected*encodes,"public workload records")
                    matches=encodes*2 if cell==0 else 0
                    equal(json.loads(path.with_suffix(".trace.json").read_text()),{
                        "schema":"gfni-inverse-distance1-timing/v1","cell":cell,"enabled":bool(mode),
                        "encodes":encodes,"calls":encodes*(2,2,1,1,2,1)[cell],"matches":matches,
                        "changed":matches if mode else 0,"timed":False,"exercise":operation=="exercise"})
                    if flavor=="release" and operation=="check":
                        parity=path.with_suffix(".parity")
                        reference=baseline / ("%d-main.parity" % cell)
                        require(parity.stat().st_size==reference.stat().st_size,"parity length")
                        with parity.open("rb") as actual, reference.open("rb") as original_bytes:
                            while True:
                                block=actual.read(65536)
                                require(block==original_bytes.read(65536),"exact Leopard1 parity")
                                if not block: break
                        total+=parity.stat().st_size
    peaks.append(resource_peak(root / "guards.log"))
    for flavor in ("release","sanitize"):
        for i in range(8):
            path=root / "guards" / ("%s-%d" % (flavor,i))
            require(path.with_suffix(".stdout").read_bytes()==b"","guard stdout")
            lines=path.with_suffix(".stderr").read_text().splitlines()
            require(len(lines)==1 and lines[0].startswith("distance-one timing driver: "),"guard stderr")
    require(total==122028032,"parity total")
    return dict(fixed_records=24,clock_free_exercise_records=104,capacity_tests=4,
        default16_release_byte_identical=True,malformed_requests=16,full_parity_bytes=total,
        validated_scopes=len(peaks),peak_bytes=max(peaks),timed_invocations=0,performance_inference=False)


if __name__=="__main__":
    require(len(sys.argv)==4,"usage: verify_gfni_inverse_distance1_timing.py ROOT BASELINE ORIGINAL")
    print(json.dumps(replay(*(Path(value) for value in sys.argv[1:])),sort_keys=True))

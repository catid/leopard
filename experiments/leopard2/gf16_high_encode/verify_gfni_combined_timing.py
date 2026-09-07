#!/usr/bin/env python3
"""Retained-only four-mode front-end qualification; executes no benchmark."""
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
            "46c9b86ced9386173cbb1aa22ebc0e2d13066984ab4ac6798bb3ac2b38dada61", "default16 binary unchanged")
    for flavor in ("release","sanitize"):
        for capacity in (16,64):
            path=root / "checks" / ("%s-capacity-%d.log" % (flavor,capacity))
            peaks.append(resource_peak(path))
            line="combined capacity %d: four modes, overflow, exact records, reset, neighbor passed" % capacity
            require(path.read_text().splitlines().count(line)==1,"capacity proof")
        for mode in range(4):
            for cell in range(6):
                for operation in (("check","exercise") if cell==0 else ("check",)):
                    path=root / "checks" / ("%s-%d-%d-%s" % (flavor,mode,cell,operation))
                    peaks.append(resource_peak(path.with_suffix(".log")))
                    encodes=26 if operation=="exercise" else 1
                    expected=(baseline / ("%d-current.json" % cell)).read_bytes()
                    require(path.with_suffix(".jsonl").read_bytes()==expected*encodes,"public workload records")
                    matches=encodes*2 if cell==0 else 0
                    equal(json.loads(path.with_suffix(".trace.json").read_text()),{
                        "schema":"gfni-combined-timing/v1","cell":cell,"mode":mode,
                        "encodes":encodes,"calls":encodes*(2,2,1,1,2,1)[cell],"matches":matches,
                        "first":matches if mode&1 else 0,"terminal":matches if mode&2 else 0,
                        "timed":False,"exercise":operation=="exercise"})
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
            require(len(lines)==1 and lines[0].startswith("combined timing driver: "),"guard stderr")
    require(total==244056064,"parity total")
    return dict(fixed_records=48,clock_free_exercise_records=208,capacity_tests=4,
        default16_release_byte_identical=True,malformed_requests=16,full_parity_bytes=total,
        validated_scopes=len(peaks),peak_bytes=max(peaks),timed_invocations=0,performance_inference=False)


if __name__=="__main__":
    require(len(sys.argv)==4,"usage: verify_gfni_combined_timing.py ROOT BASELINE ORIGINAL")
    print(json.dumps(replay(*(Path(value) for value in sys.argv[1:])),sort_keys=True))

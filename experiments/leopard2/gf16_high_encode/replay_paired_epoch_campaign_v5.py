#!/usr/bin/env python3
"""Replay v5 journal and resource footer without launching codec binaries."""
import json
from pathlib import Path
import sys

import paired_epoch_campaign_v5 as campaign
import replay_paired_epoch_campaign as replay

if __name__ == '__main__':
    if len(sys.argv) != 3:
        raise SystemExit('usage: replay_paired_epoch_campaign_v5.py ATTEMPT RESOURCE_LOG')
    replay.c = campaign._v4._base
    replay.publication = campaign.preregistration
    result = replay.verify(Path(sys.argv[1]),
                           Path(campaign.TIMING_ROOT),
                           Path(sys.argv[2]))
    print(json.dumps(result, sort_keys=True, allow_nan=False))

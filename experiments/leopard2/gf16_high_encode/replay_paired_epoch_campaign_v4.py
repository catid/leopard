#!/usr/bin/env python3
"""Replay and validate the v4 successor without launching codec binaries."""
import json
from pathlib import Path
import sys

import paired_epoch_campaign_v4 as campaign
import replay_paired_epoch_campaign as replay


if __name__ == '__main__':
    if len(sys.argv) != 3:
        raise SystemExit('usage: replay_paired_epoch_campaign_v4.py ATTEMPT RESOURCE_LOG')
    # replay's verifier resolves c.* through the inherited module object.
    replay.c.BEAD = campaign.BEAD
    replay.c.ROOT = campaign._base.ROOT
    replay.c.plan = campaign.plan
    replay.c.SOURCES = campaign._base.SOURCES
    replay.c.FILES = campaign._base.FILES
    replay.publication = campaign.preregistration
    result = replay.verify(Path(sys.argv[1]),
                           Path(campaign.TIMING_ROOT),
                           Path(sys.argv[2]))
    print(json.dumps(result, sort_keys=True, allow_nan=False))

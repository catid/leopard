#!/usr/bin/env python3
"""Run only the separately preregistered v4 paired AUTO successor."""
import sys
from pathlib import Path

import paired_epoch_campaign_v4 as campaign


if __name__ == '__main__':
    campaign._base.require(len(sys.argv) == 3,
                           'usage: run_paired_epoch_campaign_v4.py PREREGISTERED_FROZEN PUSHED_COMMIT')
    campaign.run(Path(sys.argv[1]).resolve(strict=True), sys.argv[2])

#!/usr/bin/env python3
"""Run the one-time resource-captured v5 AUTO successor."""
import sys
from pathlib import Path

import paired_epoch_campaign_v5 as campaign

if __name__ == '__main__':
    campaign.require(len(sys.argv) == 3,
                     'usage: run_paired_epoch_campaign_v5.py PREREGISTERED_FROZEN PUSHED_COMMIT')
    campaign.run(Path(sys.argv[1]).resolve(strict=True), sys.argv[2])

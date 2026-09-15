"""Resource-captured successor configuration for AUTO R19932 integration."""
import subprocess

import paired_epoch_campaign_v4 as _v4

BEAD = 'leopard-79h.38.5.4.19.1.6'
PLAN = 'paired_epoch_campaign_v5_plan.json'
ROOT = '/home/catid/leopard/.research/leopard-79h/paired-epoch-campaign-v5'
ATTEMPT_ROOT = '/tmp/leopard-paired-epoch-successor-20260915-v5/attempt1'
TIMING_ROOT = ROOT + '/frozen'
QUALIFIED_ROOT = _v4.QUALIFIED_ROOT

# Reuse only the already-qualified implementation; all v5 paths and source
# identities are distinct and are pinned in the new bundle.
_v4.BEAD = BEAD
_v4.PLAN = PLAN
_v4.ROOT = ROOT
_v4.ATTEMPT_ROOT = ATTEMPT_ROOT
_v4.TIMING_ROOT = TIMING_ROOT
_v4.QUALIFIED_ROOT = QUALIFIED_ROOT
_v4._base.BEAD = BEAD
_v4._base.PLAN = PLAN
_v4._base.ROOT = _v4._base.Path(ROOT)
_v4._base.SOURCES = set(_v4._base.SOURCES) | {
    'paired_epoch_campaign_v5.py', 'run_paired_epoch_campaign_v5.py',
    'replay_paired_epoch_campaign_v5.py'}
_v4._base.SOURCE_ROOTS = _v4._base.qualified.FINAL_ROOTS + tuple(sorted(
    set(_v4._base.NEW_SOURCES) | {
        'paired_epoch_campaign_v4.py', 'run_paired_epoch_campaign_v4.py',
        'replay_paired_epoch_campaign_v4.py', 'paired_epoch_campaign_v5.py',
        'run_paired_epoch_campaign_v5.py', 'replay_paired_epoch_campaign_v5.py'}))
_v4._base.EXTRA = set(_v4._base.EXTRA) - {_v4._original['PLAN']} | {PLAN}
_v4._base.FILES = (_v4._base.SOURCES | _v4._base.ASSETS |
                   set(_v4._base.PROOFS) | _v4._base.EXTRA |
                   {'native', 'current'})


def plan(ready=False):
    value = _v4.plan(ready)
    value['frontend_commit'] = 'paired-epoch-v5-resource-capture'
    return value


def preregistration(bundle, commit):
    return _v4.preregistration(bundle, commit)


def inputs(bundle, **kwargs):
    return _v4.inputs(bundle, **kwargs)


def qualification_gate(pins):
    return _v4.qualification_gate(pins)


_v4._base.plan = plan


def controller_command(bundle, output, qualification, commit=None):
    script = 'freeze_paired_epoch_campaign.py' if qualification else \
        'run_paired_epoch_campaign_v5.py'
    args = ['check', str(bundle), str(output)] if qualification else \
        [str(bundle), commit]
    return ['/usr/bin/timeout', '--signal=TERM', '--kill-after=5', '1200',
            '/usr/bin/prlimit', '--cpu=600:600', '--core=0:0', '--',
            '/usr/bin/python3', '-B', str(bundle / script), *args]


_v4._base.controller_command = controller_command
_v4._runner.c = _v4._base
_v4._runner.preregistration = preregistration
_v4._runner.qualification_gate = qualification_gate
_v4._runner.analysis.BEAD = BEAD

run = _v4.run
require = _v4._base.require

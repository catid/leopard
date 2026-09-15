"""Distinct v4 successor configuration for the paired AUTO R19932 screen.

The v3 campaign and its consumed attempt remain immutable.  This adapter keeps
the qualified v3 implementation and changes only the Beads identity and fresh
paths, so a new readiness-true preregistration cannot accidentally reuse the
old attempt directory.
"""
import subprocess

import paired_epoch_campaign as _base
import run_paired_epoch_campaign as _runner

BEAD = 'leopard-79h.38.5.4.19.1.5'
QUALIFICATION_BEAD = 'leopard-79h.38.5.4.19.1.4.4'
PLAN = 'paired_epoch_campaign_v4_plan.json'
ROOT = '/home/catid/leopard/.research/leopard-79h/paired-epoch-campaign-v4'
QUALIFIED_ROOT = '/home/catid/leopard/.research/leopard-79h/paired-epoch-campaign-v3.fSHScp'
ATTEMPT_ROOT = '/tmp/leopard-paired-epoch-successor-20260915/attempt1'
TIMING_ROOT = ROOT + '/frozen'

_old_plan = _base.plan
_old_inputs = _base.inputs
_original = dict(BEAD=_base.BEAD, PLAN=_base.PLAN, SOURCES=set(_base.SOURCES),
                 SOURCE_ROOTS=_base.SOURCE_ROOTS, FILES=set(_base.FILES),
                 EXTRA=set(_base.EXTRA))


def plan(ready=False):
    value = _old_plan(ready)
    value.update({
        'bead': BEAD,
        'attempt_root': ATTEMPT_ROOT,
        'timing_path': TIMING_ROOT,
        'qualification_path': QUALIFIED_ROOT + '/qualification-v2',
        'qualification_output': QUALIFIED_ROOT + '/qualification-checks-v2',
        'qualification_resource': QUALIFIED_ROOT + '/qualification-resource-v2.log',
        'ready_to_time': bool(ready),
        'frontend_commit': 'paired-epoch-v4-adapter',
    })
    return value


# The inherited helpers resolve their globals through paired_epoch_campaign.
# Override those globals before the runner calls inputs() or qualification_gate().
_base.BEAD = BEAD
_base.PLAN = PLAN
_base.ROOT = _base.Path(ROOT)
_base.plan = plan
_base.SOURCES = set(_base.SOURCES) | {'paired_epoch_campaign_v4.py',
                                       'run_paired_epoch_campaign_v4.py',
                                       'replay_paired_epoch_campaign_v4.py'}
_base.SOURCE_ROOTS = _base.qualified.FINAL_ROOTS + tuple(sorted(
    set(_base.NEW_SOURCES) | {'paired_epoch_campaign_v4.py',
                              'run_paired_epoch_campaign_v4.py',
                              'replay_paired_epoch_campaign_v4.py'}))
_base.FILES = (_base.SOURCES | _base.ASSETS | set(_base.PROOFS) |
               (set(_base.EXTRA) - {_original['PLAN']} | {PLAN}) |
               {'native', 'current'})
_base.EXTRA = set(_base.EXTRA) - {_original['PLAN']} | {PLAN}


def preregistration(bundle, commit):
    """Require the new commit to be published on the integration branch."""
    _base.require(type(commit) is str and len(commit) == 40 and
                  all(v in '0123456789abcdef' for v in commit), 'commit')
    for name in _base.SOURCES | _base.ASSETS | _base.EXTRA:
        # Only source files under gf16_high_encode are part of this bundle.
        try:
            data = subprocess.check_output(
                ['git', 'show', commit + ':experiments/leopard2/gf16_high_encode/' + name],
                timeout=15)
        except subprocess.CalledProcessError:
            continue
        path = _base.Path(bundle) / name
        if path.is_file():
            _base.equal(data.decode(), path.read_text())
    subprocess.run(['git', 'merge-base', '--is-ancestor', commit,
                    'origin/master'], check=True, timeout=15)


# Patch only the runner's module references; its implementation remains the
# previously qualified collector and is included in the immutable bundle.
_runner.c = _base
_runner.preregistration = preregistration


def inputs(bundle, **kwargs):
    # The inherited qualification proof belongs to its original untimed Beads
    # task.  Validate it under that identity, then restore this successor's
    # identity before the caller records timing state.
    previous = {key: getattr(_base, key) for key in
                ('BEAD', 'PLAN', 'SOURCES', 'SOURCE_ROOTS', 'FILES', 'EXTRA')}
    is_v4 = (_base.Path(bundle) / PLAN).is_file()
    _base.BEAD = QUALIFICATION_BEAD
    if not is_v4:
        _base.PLAN = _original['PLAN']
        _base.SOURCES = _original['SOURCES']
        _base.SOURCE_ROOTS = _original['SOURCE_ROOTS']
        _base.FILES = _original['FILES']
        _base.EXTRA = _original['EXTRA']
    try:
        return _old_inputs(bundle, **kwargs)
    finally:
        for key, value in previous.items():
            setattr(_base, key, value)


_base.inputs = inputs


_old_qualification_gate = _runner.qualification_gate


def qualification_gate(pins):
    previous = _base.BEAD
    _base.BEAD = QUALIFICATION_BEAD
    try:
        return _old_qualification_gate(pins)
    finally:
        _base.BEAD = previous


_runner.qualification_gate = qualification_gate

# Re-export the implementation API used by the v4 runner.
inputs = inputs
run = _runner.run
preflights = _base.preflights
command = _base.command
condition_env = _base.condition_env
executing_sources = _base.executing_sources

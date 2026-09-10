"""Adversarial replay of real retained records, never codec execution."""
import copy
import os
from pathlib import Path
import unittest
from unittest import mock

import verify_paired_metadata as verifier


@unittest.skipUnless(os.environ.get('PAIRED_METADATA_TEST_ROOT'),'requires explicitly selected retained evidence')
class RetainedReplayTests(unittest.TestCase):
    def corrupt(self,name,mutation):
        root=Path(os.environ['PAIRED_METADATA_TEST_ROOT']).resolve(strict=True)
        text=(root/name).read_text(); original=verifier.parse
        changed=copy.deepcopy(original(text)); mutation(changed)
        def parse(value): return copy.deepcopy(changed) if value==text else original(value)
        with mock.patch.object(verifier,'parse',side_effect=parse), self.assertRaises(ValueError):
            verifier.replay(root)

    def test_changed_command_lock(self):
        def change(state):
            args=state['records'][0]['scope_command']
            args[args.index('/tmp/leopard-gf8-authoritative.lock')]='/tmp/wrong-lock'
        self.corrupt('checks/checks.json',change)

    def test_removed_cpu_bound(self):
        self.corrupt('checks/checks.json',lambda s:s['records'][0]['scope_command'].remove('--cpu=60:60'))

    def test_forged_memory_peak(self):
        self.corrupt('checks/checks.json',lambda s:s['records'][0].__setitem__('memory_peak',1))

    def test_missing_executable_digest(self):
        self.corrupt('build/build.json',lambda s:s['artifacts'].pop('native/abort'))

    def test_empty_build_inventory(self):
        self.corrupt('build/build.json',lambda s:s.__setitem__('artifacts',{}))

    def test_missing_source_input(self):
        def change(state):
            state['inputs'].pop(next(n for n in state['inputs'] if n.endswith('/PairedRuntimeMetadata.h')))
        self.corrupt('build/build.json',change)

    def test_empty_tool_inventory(self):
        self.corrupt('checks/checks.json',lambda s:s.__setitem__('tools',{}))

    def test_unit_wrong_source(self):
        def change(state):
            command=state['commands'][0]; command[command.index('-c')+1]='wrong.cpp'
        self.corrupt('units/build.json',change)

    def test_unit_missing_sanitizer(self):
        self.corrupt('units/build.json',lambda s:s['commands'][4].remove('-fsanitize=address,undefined'))

    def test_unit_wrong_object(self):
        def change(state):
            command=state['commands'][1]
            index=next(i for i,a in enumerate(command) if a.endswith('/units/native.o'))
            command[index]=command[index].replace('/units/native.o','/build/native/driver.o')
        self.corrupt('units/build.json',change)


if __name__=='__main__': unittest.main()

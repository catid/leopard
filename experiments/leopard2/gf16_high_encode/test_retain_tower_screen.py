"""Small retention-only tests; no codecs, benchmark clocks or source artifacts."""
import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import retain_tower_screen as retain


class RetentionTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory(prefix='tower-retention-unit-')
        self.root=Path(self.temp.name)
        self.source=self.root/'source'; self.source.mkdir()
        self.output=self.root/'output'; self.output.mkdir(mode=0o700)
        (self.source/'attempt1').mkdir()
        (self.source/'attempt1/attempt.json').write_text('{"complete":false,"failure":"test"}\n')
        (self.source/'executable').write_bytes(b'qualified bytes\0')
        (self.source/'executable').chmod(0o555)

    def tearDown(self):
        for path in self.root.rglob('*'):
            if path.is_dir() and not path.is_symlink(): path.chmod(0o700)
        self.temp.cleanup()

    def test_private_readonly_copies_and_manifest(self):
        out=io.StringIO()
        with contextlib.redirect_stdout(out): retain.retain(self.source,self.output)
        data=json.loads(out.getvalue()); self.assertEqual(data['files'],3)
        for original in (self.source/'executable',self.source/'attempt1/attempt.json'):
            target=self.output/original.relative_to(self.source)
            self.assertEqual(target.read_bytes(),original.read_bytes())
            self.assertNotEqual((target.stat().st_dev,target.stat().st_ino),
                                (original.stat().st_dev,original.stat().st_ino))
            self.assertEqual(target.stat().st_mode & 0o777,0o444)
        self.assertEqual(self.output.stat().st_mode & 0o777,0o555)
        self.assertEqual(data['manifest_sha256'],retain.sha(self.output/'SHA256SUMS'))

    def test_existing_destination(self):
        (self.output/'owned').write_text('preserve')
        with self.assertRaises(ValueError): retain.retain(self.source,self.output)
        self.assertEqual((self.output/'owned').read_text(),'preserve')

    def test_reserved_root_manifest_before_copy(self):
        (self.source/'SHA256SUMS').write_text('existing manifest')
        with self.assertRaisesRegex(ValueError,'reserved root manifest'):
            retain.retain(self.source,self.output)
        self.assertFalse(any(self.output.iterdir()))
        self.assertEqual((self.source/'SHA256SUMS').read_text(),'existing manifest')

    def test_root_and_member_links(self):
        link=self.root/'alias'; link.symlink_to(self.source,target_is_directory=True)
        with self.assertRaises(ValueError): retain.retain(link,self.output)
        (self.source/'alias').symlink_to(self.source/'executable')
        with self.assertRaises(ValueError): retain.retain(self.source,self.output)
        self.assertFalse(any(self.output.iterdir()))

    def test_overlapping_trees(self):
        nested=self.source/'nested'; nested.mkdir(mode=0o700)
        with self.assertRaises(ValueError): retain.retain(self.source,nested)
        with self.assertRaises(ValueError): retain.retain(self.source,self.source)

    def test_journal_and_destination_mode(self):
        missing=self.root/'missing'; missing.mkdir()
        with self.assertRaises(ValueError): retain.retain(missing,self.output)
        self.output.chmod(0o755)
        with self.assertRaises(ValueError): retain.retain(self.source,self.output)

    def test_source_changed_after_copy(self):
        original_sha=retain.sha; calls=0
        target=self.source/'attempt1/attempt.json'
        def altered(path):
            nonlocal calls
            if path==target:
                calls+=1
                if calls==3: return '0'*64
            return original_sha(path)
        with mock.patch.object(retain,'sha',side_effect=altered),self.assertRaises(ValueError):
            retain.retain(self.source,self.output)
        self.assertFalse((self.output/'SHA256SUMS').exists())


if __name__=='__main__': unittest.main()

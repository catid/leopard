"""Small private-copy and refusal tests; no codec or benchmark clocks."""
import json
from pathlib import Path
import tempfile
import unittest
import retain_paired_epoch as retain


class RetentionTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory(); self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name); self.source=self.root/'source'; self.source.mkdir(mode=0o700)
        self.out=self.root/'out'; self.out.mkdir(mode=0o700)
        for name in ('checks','units'):
            folder=self.source/name; folder.mkdir()
            (folder/'checks.json').write_text(json.dumps(dict(bead=retain.BEAD,completed=True,timed=False)))
        for mode in ('normal','optimized'): (self.source/('replay-'+mode+'.json')).write_text('{}\n')

    def writable(self,root):
        # Test-owned temporary trees only; permits TemporaryDirectory cleanup.
        root.chmod(0o700)
        for path in root.rglob('*'):
            if path.is_dir(): path.chmod(0o700)

    def test_private_readonly_copy(self):
        report=retain.retain(self.source,self.out)
        self.addCleanup(self.writable,self.out)
        self.assertEqual(report['files'],5)
        self.assertEqual(self.out.stat().st_mode&0o222,0)
        for path in self.source.rglob('*'):
            if path.is_file():
                copied=self.out/path.relative_to(self.source)
                self.assertEqual(path.read_bytes(),copied.read_bytes())
                self.assertNotEqual((path.stat().st_dev,path.stat().st_ino),(copied.stat().st_dev,copied.stat().st_ino))
                self.assertEqual(copied.stat().st_mode&0o222,0)
        self.assertEqual(retain.sha(self.out/'SHA256SUMS'),report['manifest_sha256'])

    def test_stopped_history(self):
        stopped=self.root/'old'; (stopped/'checks').mkdir(parents=True)
        (stopped/'checks/checks.json').write_text(json.dumps(dict(bead=retain.BEAD,completed=False,timed=False,
                                                               records=[dict(returncode=143)])))
        report=retain.retain(self.source,self.out,stopped)
        self.addCleanup(self.writable,self.out)
        self.assertTrue(report['stopped_attempt_included'])
        self.assertTrue((self.out/'stopped_attempt/checks/checks.json').is_file())

    def test_reserved_manifest(self):
        (self.source/'SHA256SUMS').write_text('collision')
        with self.assertRaises(ValueError): retain.retain(self.source,self.out)
        self.assertFalse(any(self.out.iterdir()))

    def test_existing_freeze_manifest(self):
        manifest=self.source/'final-tools.json'; manifest.write_text('preserve')
        with self.assertRaises(ValueError): retain.freeze(self.source)
        self.assertEqual(manifest.read_text(),'preserve')
        self.assertFalse((self.source/'final_tools').exists())

    def test_reserved_history(self):
        (self.source/'stopped_attempt').mkdir()
        with self.assertRaises(ValueError): retain.retain(self.source,self.out)

    def test_overlapping_histories(self):
        with self.assertRaises(ValueError): retain.retain(self.source,self.out,self.source)
        with self.assertRaises(ValueError): retain.retain(self.source,self.out,self.source/'nested')
        with self.assertRaises(ValueError): retain.retain(self.source,self.out,self.root)

    def test_valid_nested_stopped_history_rejected_before_copy(self):
        # Both journals are valid independently and the destination is outside
        # both sources. Without the mutual-disjointness guard, the stopped
        # evidence would be copied twice, under nested/ and stopped_attempt/.
        stopped=self.source/'nested'; (stopped/'checks').mkdir(parents=True)
        (stopped/'checks/checks.json').write_text(json.dumps(dict(
            bead=retain.BEAD,completed=False,timed=False,records=[dict(returncode=143)])))
        with self.assertRaisesRegex(ValueError,'disjoint complete and stopped histories'):
            retain.retain(self.source,self.out,stopped)
        self.assertFalse(any(self.out.iterdir()))

    def test_source_link(self):
        (self.source/'alias').symlink_to(self.source/'replay-normal.json')
        with self.assertRaises(ValueError): retain.retain(self.source,self.out)

    def test_occupied_or_nonprivate_destination(self):
        (self.out/'owned').write_text('keep')
        with self.assertRaises(ValueError): retain.retain(self.source,self.out)
        self.assertEqual((self.out/'owned').read_text(),'keep')
        (self.out/'owned').unlink(); self.out.chmod(0o755)
        with self.assertRaises(ValueError): retain.retain(self.source,self.out)

    def test_disagreeing_replays(self):
        (self.source/'replay-optimized.json').write_text('{"different":true}')
        with self.assertRaises(ValueError): retain.retain(self.source,self.out)

    def test_incomplete_checks(self):
        (self.source/'checks/checks.json').write_text(json.dumps(dict(bead=retain.BEAD,completed=False,timed=False)))
        with self.assertRaises(ValueError): retain.retain(self.source,self.out)

    def test_newline_name(self):
        (self.source/'unsafe\nname').write_text('x')
        with self.assertRaises(ValueError): retain.retain(self.source,self.out)


if __name__=='__main__': unittest.main()

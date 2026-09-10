"""Private-copy and refusal checks on test-owned evidence trees."""
import json
from pathlib import Path
import tempfile
import unittest
import retain_paired_epoch_timing as retain


class RetentionTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory(prefix='leopard-epoch-retain-'); self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name); self.source=self.root/'source'; self.failed=self.root/'failed'
        (self.source/'checks').mkdir(parents=True); (self.failed/'build').mkdir(parents=True)
        self.out=self.root/'out'; self.out.mkdir(mode=0o700)
        (self.source/'checks/checks.json').write_text(json.dumps(dict(bead=retain.BEAD,completed=True,real_clocks_read=False)))
        (self.failed/'build/build.json').write_text(json.dumps(dict(bead=retain.BEAD,completed=False,real_clocks_read=False)))
        for mode in ('normal','optimized'): (self.source/('replay-'+mode+'.json')).write_text('{}\n')

    def copy(self): return retain.retain(self.source,self.out,self.failed)

    def test_private_readonly_copy_with_failure(self):
        try:
            result=self.copy()
            self.assertEqual(result['files'],5); self.assertTrue(result['failed_build_included'])
            self.assertFalse(result['collector_qualified'])
            self.assertEqual(self.out.stat().st_mode&0o222,0)
            for folder,prefix in ((self.source,''),(self.failed,'failed_build')):
                for path in folder.rglob('*'):
                    if path.is_file():
                        other=self.out/prefix/path.relative_to(folder)
                        self.assertEqual(path.read_bytes(),other.read_bytes())
                        self.assertEqual(other.stat().st_mode&0o222,0)
                        self.assertNotEqual((path.stat().st_dev,path.stat().st_ino),(other.stat().st_dev,other.stat().st_ino))
            self.assertEqual(retain.sha(self.out/'SHA256SUMS'),result['manifest_sha256'])
        finally:
            self.out.chmod(0o700)
            for path in self.out.rglob('*'):
                if path.is_dir(): path.chmod(0o700)

    def test_reserved_names(self):
        for name in ('SHA256SUMS','failed_build'):
            path=self.source/name; path.write_text('keep')
            with self.assertRaises(ValueError): self.copy()
            self.assertFalse(any(self.out.iterdir())); self.assertEqual(path.read_text(),'keep'); path.unlink()

    def test_links_and_newlines(self):
        path=self.failed/'alias'; path.symlink_to(self.failed/'build/build.json')
        with self.assertRaises(ValueError): self.copy()
        path.unlink(); path=self.source/'bad\nname'; path.write_text('x')
        with self.assertRaises(ValueError): self.copy()
        self.assertFalse(any(self.out.iterdir()))

    def test_disjoint(self):
        nested=self.source/'nested'; (nested/'build').mkdir(parents=True)
        (nested/'build/build.json').write_text(json.dumps(dict(bead=retain.BEAD,completed=False,real_clocks_read=False)))
        # Both journals independently satisfy their respective success/failure
        # contracts; only the mutual-overlap guard rejects this double copy.
        with self.assertRaisesRegex(ValueError,'disjoint evidence trees'):
            retain.retain(self.source,self.out,nested)
        self.assertFalse(any(self.out.iterdir()))
        for source,failed in ((self.source,self.source),(self.source,self.source/'nested'),(self.root,self.failed)):
            with self.assertRaises(ValueError): retain.retain(source,self.out,failed)
        self.assertFalse(any(self.out.iterdir()))

    def test_nonprivate_and_occupied_destination(self):
        self.out.chmod(0o755)
        with self.assertRaises(ValueError): self.copy()
        self.out.chmod(0o700); (self.out/'keep').write_text('owned')
        with self.assertRaises(ValueError): self.copy()
        self.assertEqual((self.out/'keep').read_text(),'owned')

    def test_disagreeing_replays(self):
        (self.source/'replay-optimized.json').write_text('{"different":true}')
        with self.assertRaises(ValueError): self.copy()
        self.assertFalse(any(self.out.iterdir()))

    def test_incomplete_or_wrong_failure(self):
        checks=self.source/'checks/checks.json'; original=checks.read_text()
        checks.write_text(json.dumps(dict(bead=retain.BEAD,completed=False,real_clocks_read=False)))
        with self.assertRaises(ValueError): self.copy()
        checks.write_text(original); (self.failed/'checks').mkdir()
        with self.assertRaises(ValueError): self.copy()
        self.assertFalse(any(self.out.iterdir()))


if __name__=='__main__': unittest.main()

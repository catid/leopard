"""Final-source identity regressions; no native codec or benchmark clocks."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

import freeze_paired_epoch_timing as freezer
import verify_paired_epoch_timing as verify

ROOT=Path(__file__).resolve().parent


class FinalToolTests(unittest.TestCase):
    def setUp(self):
        raw=os.environ.get('LEO_PAIRED_EPOCH_TIMING_EVIDENCE')
        if not raw: self.skipTest('new native fixtures require LEO_PAIRED_EPOCH_TIMING_EVIDENCE')
        self.temp=tempfile.TemporaryDirectory(prefix='leopard-epoch-tools-'); self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name); self.folder=self.root/'final_tools'; self.folder.mkdir()
        (self.root/'build').mkdir()
        shutil.copyfile(Path(raw)/'build/reference_build.json',self.root/'build/reference_build.json')
        self.pins={}
        for name in verify.FINAL_PYTHON | set(verify.FINAL_ASSETS):
            shutil.copyfile(ROOT/name,self.folder/name); self.pins[name]=verify.sha(self.folder/name)
        self.manifest()

    def manifest(self):
        (self.root/'final-tools.json').write_text(json.dumps(dict(bead=verify.BEAD,real_clocks_read=False,files=self.pins)))

    def changed(self,name):
        path=self.folder/name; path.write_text(path.read_text()+'\n')
        self.pins[name]=verify.sha(path); self.manifest()

    def test_complete_inventory(self):
        verify.verify_final_tools(self.root)

    def test_missing_dependency_and_digest(self):
        name='verify_auto_gfni_boundary_checks.py'
        (self.folder/name).unlink(); self.pins.pop(name); self.manifest()
        with self.assertRaises(ValueError): verify.verify_final_tools(self.root)

    def test_changed_asset_with_new_digest(self):
        self.changed('paired_timer_clock.cpp')
        with self.assertRaises(ValueError): verify.verify_final_tools(self.root)

    def test_executing_verifier_identity(self):
        self.changed('verify_paired_epoch_timing.py')
        with self.assertRaises(ValueError): verify.verify_final_tools(self.root)

    def executing_check(self,source,extra=''):
        program=('import sys; from pathlib import Path; sys.path.insert(0,sys.argv[1]); '
                 'import verify_paired_epoch_timing as verify; import paired_epoch_analysis; '
                 'import replay_paired_epoch_analysis; '+extra+'verify.verify_final_tools(Path(sys.argv[2]))')
        return subprocess.run([sys.executable,'-I','-B',*([] if __debug__ else ['-O']),'-c',program,
                               str(source),str(self.root)],capture_output=True,text=True,timeout=10)

    def test_drifted_executing_dependency(self):
        source=self.root/'executing'; shutil.copytree(self.folder,source)
        self.assertEqual(self.executing_check(source).returncode,0)
        for name in ('verify_paired_metadata.py','paired_epoch_timing_elf.py',
                     'paired_epoch_analysis.py','replay_paired_epoch_analysis.py'):
            with self.subTest(name=name):
                original=(source/name).read_text(); (source/name).write_text(original+'\nDRIFTED = True\n')
                result=self.executing_check(source)
                self.assertNotEqual(result.returncode,0)
                self.assertIn('executing dependency hash: '+name,result.stderr)
                (source/name).write_text(original)

    def test_drifted_duplicate_main(self):
        source=self.root/'executing'; shutil.copytree(self.folder,source)
        alias=source/'alias'; alias.mkdir()
        (alias/'verify_paired_epoch_timing.py').write_text((source/'verify_paired_epoch_timing.py').read_text()+'\nDRIFTED=True\n')
        extra=("import importlib.util; spec=importlib.util.spec_from_file_location('epoch_alias', "
               "Path(sys.argv[1])/'alias/verify_paired_epoch_timing.py'); module=importlib.util.module_from_spec(spec); "
               "sys.modules['epoch_alias']=module; spec.loader.exec_module(module); ")
        result=self.executing_check(source,extra)
        self.assertNotEqual(result.returncode,0)
        self.assertIn('executing dependency hash: verify_paired_epoch_timing.py',result.stderr)

    def test_unlisted_directory(self):
        (self.folder/'unlisted').mkdir()
        with self.assertRaises(ValueError): verify.verify_final_tools(self.root)

    def test_freeze_occupied_preserves_files(self):
        before=(self.root/'final-tools.json').read_bytes()
        with self.assertRaises(ValueError): freezer.freeze(self.root)
        self.assertEqual((self.root/'final-tools.json').read_bytes(),before)

    def test_freeze_terminal_gate_and_success(self):
        fresh=self.root/'fresh'; (fresh/'build').mkdir(parents=True); (fresh/'checks').mkdir()
        shutil.copyfile(self.root/'build/reference_build.json',fresh/'build/reference_build.json')
        checks=fresh/'checks/checks.json'
        checks.write_text(json.dumps(dict(bead=verify.BEAD,completed=False,real_clocks_read=False)))
        with self.assertRaises(ValueError): freezer.freeze(fresh)
        self.assertFalse((fresh/'final_tools').exists())
        checks.write_text(json.dumps(dict(bead=verify.BEAD,completed=True,real_clocks_read=False)))
        try:
            freezer.freeze(fresh); verify.verify_final_tools(fresh)
            self.assertEqual((fresh/'final_tools').stat().st_mode&0o222,0)
        finally:
            if (fresh/'final_tools').exists(): (fresh/'final_tools').chmod(0o700)


if __name__=='__main__': unittest.main()

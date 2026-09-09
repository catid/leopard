"""Pure overlay/parser/native-record regression tests. No codec or clocks."""
import copy
import json
from pathlib import Path
import unittest
from avx2_adjacent_runtime import overlay,once
from audit_avx2_adjacent_runtime import loops
from audit_avx2_adjacent_schedule import isa_ceiling
from replay_avx2_adjacent_runtime import focused

STATIC = Path('/home/catid/leopard/.research/leopard-79h/avx2-adjacent-qualified.wxkx9v2b/source/Leopard2BackendAVX2.cpp')


def record(mode=0):
    return dict(schema='adjacent-runtime-focused/v1',mode=mode,trace=True,timed=False,
                calls=[[66413,0],[66585,0]] if not mode else [[0,66413],[0,66585]],
                blocks=[[66178,0],[259516,0]] if not mode else [[0,66178],[0,259516]])


def output(row):
    return 'adjacent pairs: forward=65535 accumulating=65535 boundary_accumulations=918\n'+json.dumps(row)+'\n'


def assembly():
    text = '00000000 <fixture>:\n'
    for start in (0,16):
        for offset in range(8): text += f'{start+offset:x}:\t90\tvpshufb %ymm0,%ymm1,%ymm2\n'
        text += f'{start+8:x}:\t90\tjne {start:x} <fixture>\n'
    return text


class Tests(unittest.TestCase):
    def test_exact_overlay_and_original_preserved(self):
        original = STATIC.read_text(); actual = overlay(original)
        self.assertEqual(STATIC.read_text(),original)
        self.assertEqual(actual.count('const bool enabled ='),2)
        self.assertIn('const bool enabled = !Inverse && leo_adjacent_schedule_enabled;',actual)
        self.assertEqual(actual.count('LeoAdjacentRecord('),2)
        self.assertEqual(actual.count('AVX2FF16AdjacentProductAdd('),2)
        self.assertIn('AVX2FF16Butterfly2PreparedImpl<Inverse, false>',actual)
        self.assertIn('AVX2FF16IFFTButterfly2XorImpl<false>',actual)
        # Both scalar-tail bodies survive verbatim; only vector schedules differ.
        a = original.index('    const uint64_t symbols = (byte_count - offset) / 2;',original.index('static void AVX2FF16IFFTButterfly2Xor('))
        tail = original[a:original.index('\n#endif // LEO_HAS_FF16',a)]
        self.assertIn(tail,actual)

    def test_overlay_rejects_wrong_or_duplicate_input(self):
        original = STATIC.read_text()
        for bad in (original+'\n',original.replace('byte_count','other'),overlay(original),''):
            with self.assertRaises(ValueError): overlay(bad)
        for text in ('','xx'):
            with self.assertRaises(ValueError): once(text,'x','y')

    def test_two_disjoint_loops(self):
        found = loops(assembly(),'fixture')
        self.assertEqual(len(found),2)
        self.assertEqual([v['instructions'] for v in found],[9,9])
        self.assertEqual([v['stack_references'] for v in found],[[],[]])

    def test_loops_refuse_missing_duplicate_or_control_inside(self):
        text = assembly()
        for bad in (text.replace('<fixture>:','<missing>:'),text+text,
                    text.replace('vpshufb %ymm0,%ymm1,%ymm2','vpxor %ymm0,%ymm1,%ymm2'),
                    text.replace('jne 0 <fixture>','jne 0 <leo_adjacent_schedule_enabled>')):
            with self.assertRaises(ValueError): loops(bad,'fixture')

    def test_full_isa_rejections(self):
        for raw,op in [('62 01','vpxor %ymm0,%ymm1,%ymm2'),('90','vpxor %ymm16,%ymm1,%ymm2'),
                       ('90','vmovdqa %zmm0,%zmm1'),('90','vpternlogd $1,%ymm0,%ymm1,%ymm2'),
                       ('90','vgf2p8affineqb $0,%ymm0,%ymm1,%ymm2')]:
            with self.assertRaises(ValueError): isa_ceiling('0:\t'+raw+'\t'+op+'\n')

    def test_mode_and_profile_records(self):
        expected = dict(calls=[66413,66585],blocks=[66178,259516])
        for p in ('trace','sanitize'):
            for m in ('off','on'): self.assertEqual(focused(output(record(int(m=='on'))),p,m,'adjacent'),expected)
        row = dict(record(),trace=False,calls=[[0,0],[0,0]],blocks=[[0,0],[0,0]])
        self.assertEqual(focused(output(row),'release','off','adjacent'),dict(calls=[0,0],blocks=[0,0]))

    def test_record_identity_failures(self):
        row = record()
        for key,value in [('schema','other'),('mode',1),('mode',False),('trace',False),('timed',True),('extra',1)]:
            with self.assertRaises(ValueError): focused(output(dict(row,**{key:value})),'trace','off','adjacent')
        for args in [('wrong','off','adjacent'),('trace','bad','adjacent'),('trace','off','bad')]:
            with self.assertRaises(ValueError): focused(output(row),*args)
        with self.assertRaises(ValueError): focused(output(row).replace('65535','65534',1),'trace','off','adjacent')

    def test_counter_types_inventory_state_and_missing_calls(self):
        for value in ([],[[1,0]],[[1,0,0],[1,0]],[[True,0],[1,0]],[[1.0,0],[1,0]],
                      [[-1,0],[1,0]],[[2**64,0],[1,0]],[[66413,1],[66585,0]],[[0,0],[0,0]],
                      [[66413,0],[1,0]]):
            with self.assertRaises(ValueError): focused(output(dict(record(),calls=value)),'trace','off','adjacent')
        bad = copy.deepcopy(record()); bad['blocks'][1][1]=1
        with self.assertRaises(ValueError): focused(output(bad),'trace','off','adjacent')


if __name__=='__main__': unittest.main()

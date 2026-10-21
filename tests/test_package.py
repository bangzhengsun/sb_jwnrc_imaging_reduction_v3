import tempfile
import unittest
from pathlib import Path
from sbnrc.cli import grid, main, verify

class PackageTests(unittest.TestCase):
    def test_non_square_grid(self):
        c={'mosaic':{'shape_60':[8000,10000],'pixel_scale_60_arcsec':0.06,'center_ra_deg':12.,'center_dec_deg':-5.,'rotation_deg':19.}}
        g=grid(c,'30mas')
        self.assertEqual(g['shape'],[16000,20000])
        self.assertEqual(g['crpix_step'],[7999,9999])
        self.assertEqual(g['pixel_scale_arcsec'],0.03)
        self.assertEqual(g['crval'],[12.,-5.])
        self.assertEqual(g['rotation'],19.)
    def test_snapshot_integrity(self):
        self.assertGreater(verify(),50)
    def test_export_preserves_existing_directory(self):
        with tempfile.TemporaryDirectory() as tmp:
            dest=Path(tmp)/'record'
            main(['export',str(dest)])
            self.assertTrue((dest/'nexus/scripts/SBNRC_01_stage1.py').exists())
            self.assertFalse((dest/'skill').exists())
            marker=dest/'keep';marker.write_text('original')
            with self.assertRaises(SystemExit):main(['export',str(dest)])
            self.assertEqual(marker.read_text(),'original')

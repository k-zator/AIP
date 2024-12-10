import unittest
from filereader.aip_reader import AipReader
import os
FIXTURE_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'test_files', '')

aip_water = f"{FIXTURE_DIR}/water_ssip.xml"

class TestAipReader(unittest.TestCase):
    def setUp(self):
        self.aip = AipReader(aip_water)

    def test_list_atoms_length(self):
        self.assertEqual(len(self.aip.list_atoms), 3)
   
    def test_atom_names(self):
        list_atoms = ['a1', "a2", "a3"]
        for i in self.aip.list_atoms:
            self.assertTrue(i.aname in list_atoms)

    def test_list_bonds_length(self):
        self.assertEqual(len(self.aip.list_bonds), 2)

    def test_list_aips_length(self):
        self.assertEqual(len(self.aip.list_aips), 4)

    def test_aip_info(self):
        dict_values = {
                i.atom_neigh.aname: i.value for i in self.aip.list_aips
                }
        self.assertAlmostEqual(dict_values["a2"], 2.7305060481279995)
        aip_point = self.aip.list_aips[0]
        self.assertAlmostEqual(aip_point.atom_type, "O.3.water" )
        self.assertEqual(aip_point.atom_type, "O.3.water" )
        self.assertAlmostEqual(aip_point.mepsvalue, -0.088262 )
        self.assertAlmostEqual(aip_point.fraction, 1.0 )

    def test_surface(self):
        self.assertAlmostEqual(self.aip.surface.total, 37.249475717728636)
        self.assertAlmostEqual(self.aip.surface.positive, 19.19858898634103)
        self.assertAlmostEqual(self.aip.surface.negative, 18.05088673138761)
        self.assertAlmostEqual(self.aip.surface.positive_polar, 9.472175572526938)
        self.assertAlmostEqual(self.aip.surface.negative_polar, 18.05088673138761)
        self.assertAlmostEqual(self.aip.surface.total_nonpolar, 9.726413413814088)
        self.assertAlmostEqual(self.aip.surface.positive_nonpolar, 9.726413413814088)
        self.assertAlmostEqual(self.aip.surface.negative_nonpolar, 0.0)
        self.assertEqual(self.aip.surface.numberOFMEPSPoints, 5128)
        self.assertAlmostEqual(self.aip.surface.electrostaticPotentialMax, 0.084445)
        self.assertAlmostEqual(self.aip.surface.electrostaticPotentialMin, -0.070811)


if __name__ == '__main__':
    unittest.main()

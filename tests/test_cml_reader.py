import unittest
from filereader.cml_reader import CmlReader
import os
FIXTURE_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'test_files', '')

cml_water = f"{FIXTURE_DIR}/water_ns.cml"
cml_ethanol = f"{FIXTURE_DIR}/ethanolnonamespace.cml"
cml_nitromethane = f"{FIXTURE_DIR}/nitromethane.cml"


class TestCmlReader(unittest.TestCase):
    def setUp(self):
        self.cml = CmlReader(cml_water)
        self.no_namespace_cml = CmlReader(cml_ethanol, ns=None)
        self.cml_nitromethane = CmlReader(cml_nitromethane, ns=None)

    def tearDown(self):
        del self.cml
        del self.no_namespace_cml

    def test_list_atoms_length(self):
        self.assertEqual(len(self.cml.list_atoms), 3)
   
    def test_atom_names(self):
        list_atoms = ['a1', "a2", "a3"]
        for i in self.cml.list_atoms:
            self.assertTrue(i.aname in list_atoms)

    def test_list_bonds_length(self):
        self.assertEqual(len(self.cml.list_bonds), 2)

    def test_no_namespace_cml(self):
        atom = self.no_namespace_cml.list_atoms[0]
        self.assertEqual(atom.aname, "a1")
        self.assertEqual(atom.element, "C")
        self.assertAlmostEqual(atom.xyz[0], -0.888300)

    def test_no_namespace_cml(self):
        atom = self.cml_nitromethane.list_atoms[0]
        self.assertEqual(atom.aname, "a1")
        self.assertEqual(atom.element, "N")
        self.assertAlmostEqual(atom.xyz[0], 0.060000)

if __name__ == '__main__':
    unittest.main()

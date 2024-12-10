import unittest
from filereader.cml_reader import CmlReader
from filereader.mol2reader import Mol2Reader
import os
from lxml import etree
import networkx as nx
FIXTURE_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'test_files', '')

cml_ethanol = f"{FIXTURE_DIR}/ethanolnonamespace.cml"
mol2_ethanol = f"{FIXTURE_DIR}/ethanol.mol2"
cml_sybyl = f"{FIXTURE_DIR}/ethanolnonamespacewithsybyl.cml"
cml_namespace = f"{FIXTURE_DIR}/ethanolnamespace.cml"
cml_namespace_sybyl = f"{FIXTURE_DIR}/ethanolnamespacesybyl.cml"
cml_nitromethane = f"{FIXTURE_DIR}/nitromethane.cml"
mol2_nitromethane = f"{FIXTURE_DIR}/nitromethane.mol2"
cml_water = f"{FIXTURE_DIR}/water_no_ns.cml"
mol2_water = f"{FIXTURE_DIR}/water.mol2"


class TestCmlMol2Comparison(unittest.TestCase):
    def setUp(self):
        self.cml = CmlReader(cml_ethanol, ns=None)
        self.mol2 = Mol2Reader(mol2_ethanol, from_file=True)
        self.cml_sybyl = CmlReader(cml_sybyl, ns=None)
        self.cml_ns = CmlReader(cml_namespace)
        self.cml_ns_sybyl = CmlReader(cml_namespace_sybyl)
        self.cml_nitromethane = CmlReader(cml_nitromethane, ns=None)
        self.mol2_nitromethane = Mol2Reader(mol2_nitromethane, from_file=True)
        self.cml_water = CmlReader(cml_water, ns=None)
        self.mol2_water = Mol2Reader(mol2_water, from_file=True)
    def tearDown(self):
        del self.cml
        del self.mol2
        del self.cml_sybyl
        del self.cml_ns
        del self.cml_ns_sybyl
        del self.cml_nitromethane
        del self.mol2_nitromethane
        del self.cml_water
        del self.mol2_water

    def test_list_atoms_length(self):
        atom0 = self.cml.list_atoms[0]
        atom1 = self.mol2.list_atoms[0]
        self.assertTrue(self.cml._compare_atoms_(atom0, atom1))

    def test_sybyl_atoms(self):
        sybyl_dict = self.cml._get_sybyl_dict_(self.mol2)
        atom_type_0 = sybyl_dict["a1"]
        self.assertEqual(len(sybyl_dict.keys()), 9)
        self.assertEqual(atom_type_0, "C.3")

    def test_sybyl_nonamespace(self):
        tree_new_sybyl = self.cml.get_cml_with_sybyl_info(self.mol2)
        tree_sybyl = self.cml_sybyl.tree
        self.assertEqual(
                etree.tostring(tree_sybyl), 
                etree.tostring(tree_new_sybyl)
                )

    def test_sybyl_namespace(self):
        tree_new_sybyl = self.cml_ns.get_cml_with_sybyl_info(self.mol2)
        tree_sybyl = self.cml_ns_sybyl.tree
        self.assertEqual(
                etree.tostring(tree_sybyl), 
                etree.tostring(tree_new_sybyl)
                )
   
    def test_list_atoms_length_nitromethane(self):
        atom0 = self.cml_nitromethane.list_atoms[0]
        atom1 = self.mol2_nitromethane.list_atoms[0]
        self.assertTrue(self.cml._compare_atoms_(atom0, atom1))

    def test_water_assignment(self):
        network = self.cml_water.get_cml_network_with_sybyl(self.mol2_water)
        dict_element = nx.get_node_attributes(network,'elementType')
        dict_sybyl =  nx.get_node_attributes(network,'sybyl')
        self.assertEqual(dict_element["a1"], "O")
        self.assertEqual(dict_sybyl["a1"], "O.3")

if __name__ == '__main__':
    unittest.main()

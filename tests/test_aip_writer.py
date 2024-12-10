import unittest
from filewriter.AIP_writer import AIPWriter, CML_NAMESPACE_DICT
from aip_footprinting.AIP_class import AIP
from aip_footprinting.surface_class import Surface
from aip_footprinting.read_MEPS import MEPS
from aip_footprinting.footprinting import Footprinting
from aip_footprinting.atom_class import AtomSet
import os
from lxml import etree
import copy
import tempfile

FIXTURE_DIR = os.path.join(os.path.dirname(
    os.path.realpath(__file__)), 'test_files', '')

cml_nw = f"{FIXTURE_DIR}/methane_nwchem_out.cml"
cube_nw = f"{FIXTURE_DIR}/methane_merged.cube"
cube_middle = f"{FIXTURE_DIR}/methane_0104_merged.cube"


class TestAIPWriter(unittest.TestCase):
    def setUp(self):
        self.aip_writer = AIPWriter(cml_nw)
        aip1 = AIP()
        self.value = 1.0
        self.x = 0
        self.atom = "a"
        self.type = "H.O"
        self.MEP_value = "0.001"
        self.isosurface = "0.0300"
        self.aip_area_fraction = "0.5"
        aip1.set_value(self.value)
        aip1.set_xyz(self.x, self.x, self.x)
        aip1.set_nearest_atom(self.atom)
        aip1.set_atom_type(self.type)
        aip1.set_mepsvalue(self.MEP_value)
        aip1.set_isosurface(self.isosurface)
        aip1.set_AIP_area_fraction(self.aip_area_fraction)
        aip2 = copy.deepcopy(aip1)
        aip2.value = 2.0
        self.list_aips = [aip1, aip2]
        meps = MEPS(cube_nw, cml_nw)
        meps_middle = MEPS(cube_middle, cml_nw)
        atom = AtomSet()
        atom._create_from_MEPS(meps)
        aip = Footprinting(meps, 0, meps_middle, atom)
        self.surface = Surface(meps, aip, atom)

    def tearDown(self):
        del self.aip_writer
        del self.surface
        del self.value
        del self.atom
        del self.type
        del self.MEP_value
        del self.isosurface
        del self.list_aips

    def test_add_aip_information(self):
        self.new_tree = self.aip_writer.add_aip_information(
            self.aip_writer.tree, self.list_aips)
        element = self.new_tree.xpath(
            "//ssip:SSIP", namespaces=CML_NAMESPACE_DICT)[0]
        self.assertEqual(element.attrib["{{{}}}nearestAtomID".format(
            CML_NAMESPACE_DICT["ssip"])], str(self.atom))
        self.assertEqual(element.attrib["{{{}}}x3".format(
            CML_NAMESPACE_DICT["cml"])], str(self.x))
        self.assertEqual(element.attrib["{{{}}}y3".format(
            CML_NAMESPACE_DICT["cml"])], str(self.x))
        self.assertEqual(element.attrib["{{{}}}z3".format(
            CML_NAMESPACE_DICT["cml"])], str(self.x))
        self.assertEqual(element.attrib["{{{}}}MEPvalue".format(
            CML_NAMESPACE_DICT["ssip"])], str(self.MEP_value))
        self.assertEqual(element.attrib["{{{}}}value".format(
            CML_NAMESPACE_DICT["ssip"])], str(self.value))
        self.assertEqual(element.attrib["{{{}}}aipAtomType".format(
            CML_NAMESPACE_DICT["ssip"])], str(self.type))
        self.assertEqual(element.attrib["{{{}}}isosurface".format(
            CML_NAMESPACE_DICT["ssip"])], str(self.isosurface))
        self.assertEqual(element.attrib["{{{}}}aipAreaFraction".format(
            CML_NAMESPACE_DICT["ssip"])], str(self.aip_area_fraction))

    def test_aip_string(self):
        tp = tempfile.NamedTemporaryFile(
            mode="w", suffix=".xml", prefix="myname")
        self.aip_writer.write_file(self.surface, self.list_aips,  tp.name)
        #self.aip_writer.write_file(self.surface, self.list_aips,  "example_ssips.xml")
        parser = etree.XMLParser(remove_blank_text=True)
        parsed_tree = copy.deepcopy(etree.parse(tp.name, parser))
        element = parsed_tree.xpath(
            "//ssip:SSIP", namespaces=CML_NAMESPACE_DICT)[0]
        self.assertEqual(element.attrib["{{{}}}nearestAtomID".format(
            CML_NAMESPACE_DICT["ssip"])], str(self.atom))
        self.assertEqual(element.attrib["{{{}}}x3".format(
            CML_NAMESPACE_DICT["cml"])], str(self.x))
        self.assertEqual(element.attrib["{{{}}}y3".format(
            CML_NAMESPACE_DICT["cml"])], str(self.x))
        self.assertEqual(element.attrib["{{{}}}z3".format(
            CML_NAMESPACE_DICT["cml"])], str(self.x))
        self.assertEqual(element.attrib["{{{}}}MEPvalue".format(
            CML_NAMESPACE_DICT["ssip"])], str(self.MEP_value))
        self.assertEqual(element.attrib["{{{}}}value".format(
            CML_NAMESPACE_DICT["ssip"])], str(self.value))
        self.assertEqual(element.attrib["{{{}}}aipAtomType".format(
            CML_NAMESPACE_DICT["ssip"])], str(self.type))
        self.assertEqual(element.attrib["{{{}}}isosurface".format(
            CML_NAMESPACE_DICT["ssip"])], str(self.isosurface))
        self.assertEqual(element.attrib["{{{}}}aipAreaFraction".format(
            CML_NAMESPACE_DICT["ssip"])], str(self.aip_area_fraction))
        tp.close()


if __name__ == '__main__':
    unittest.main()

#    nwchemcmlutils can generate MEPS files with NWChem and SSIP descriptions.
#    Copyright (C) 2019  Mark D. Driver
#
#    nwchemcmlutils is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Affero General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Affero General Public License for more details.
#
#    You should have received a copy of the GNU Affero General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 09:00:10 2016

@author: mark
"""

import unittest
import logging
import copy
import pathlib
from lxml import etree
import nwchemcmlutils.nwUtils.cml_reading_util as cml_reading_util

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.INFO)


class CMLReadingUtilTestCase(unittest.TestCase):
    """Test class for the cml_reading_util script.
    """

    def setUp(self):
        """setUp for the class.
        """
        self.maxDiff = None
        parent_directory = pathlib.Path(__file__).parents[1]
        water_cml_file = etree.parse(
            (parent_directory / "resources/examplecml.xml").absolute().as_posix()
        )
        self.cml_namespace = {"cml": "http://www.xml-cml.org/schema"}
        water_cml = water_cml_file.findall(
            "cml:molecule", namespaces=self.cml_namespace
        )[0]
        self.water_cml = water_cml
        self.bad_water_1 = water_cml_file.findall(
            "cml:molecule", namespaces=self.cml_namespace
        )[1]
        self.bad_water_2 = water_cml_file.findall(
            "cml:molecule", namespaces=self.cml_namespace
        )[2]
        bad_atom_array = cml_reading_util.extractAtomArray(self.bad_water_2)
        self.bad_oxygen = cml_reading_util.createAtomList(bad_atom_array)[0]
        self.atom_array = cml_reading_util.extractAtomArray(self.water_cml)
        self.bond_array = cml_reading_util.extractBondArray(self.water_cml)
        self.atom_list = cml_reading_util.createAtomList(self.atom_array)
        self.oxygen_atom = copy.deepcopy(self.atom_list[0])
        self.atom_attrib = cml_reading_util.getAtomListAttributes(self.atom_list)

    def tearDown(self):
        """tearDown for the class.
        """
        del self.cml_namespace
        del self.water_cml
        del self.atom_array
        del self.bond_array
        del self.atom_list
        del self.atom_attrib
        del self.bad_water_1
        del self.bad_water_2

    def test_molecule_cml(self):
        """test to check that the correct molecule is being used in the tests.
        """
        expected_molecule = (
            u'<cml:molecule xmlns:cml="http://www.xml-cml.o'
            + 'rg/schema" xmlns:ssip="http://www-hunter.ch.cam.ac.uk/SSIP" ssip:stdInChIKey="XLYOFNOQVPJJNP-UHFFFAO'
            + 'YSA-N" StuartId="13_01_20121_new" cml:id="XLYOFNOQ'
            + 'VPJJNP-UHFFFAOYSA-N">\n <cml:atomArray>\n  <cm'
            + 'l:atom cml:elementType="O" cml:id="a1" cml:x3="-0.198000" cml:y'
            + '3="0.000000" cml:z3="0.347000"/>\n  <cml:atom cml:elem'
            + 'entType="H" cml:id="a2" cml:x3="0.760000" cml:y3="0.000000'
            + '" cml:z3="0.204000"/>\n  <cml:atom cml:elementType="H"'
            + ' cml:id="a3" cml:x3="-0.561000" cml:y3="0.000000" cml:z3="-0.5'
            + '51000"/>\n </cml:atomArray>\n <cml:bondArray>\n'
            + '  <cml:bond cml:atomRefs2="a1 a2" cml:order="1"/>\n  '
            + '<cml:bond cml:atomRefs2="a1 a3" cml:order="1"/>\n </cml'
            + ":bondArray>\n</cml:molecule>\n"
        )
        actual_molecule = etree.tounicode(self.water_cml)
        self.assertEqual(actual_molecule, expected_molecule)

    def test_extract_atom_array(self):
        """Test to see that the correct atom array is extracted.
        """
        expected_atom_array = (
            '<cml:atomArray xmlns:cml="http://www.xml-cml'
            + '.org/schema" xmlns:ssip="http://www-hunter.ch.cam.ac.uk/SSIP">\n  <cml:atom cml:elementType="O" cml:i'
            + 'd="a1" cml:x3="-0.198000" cml:y3="0.000000" cml:z3="0.34'
            + '7000"/>\n  <cml:atom cml:elementType="H" cml:id="a2"'
            + ' cml:x3="0.760000" cml:y3="0.000000" cml:z3="0.204000"/>'
            + '\n  <cml:atom cml:elementType="H" cml:id="a3" cml:x3="-0'
            + '.561000" cml:y3="0.000000" cml:z3="-0.551000"/>\n </'
            + "cml:atomArray>\n "
        )
        actual_atom_array = etree.tounicode(self.atom_array)
        self.assertMultiLineEqual(actual_atom_array, expected_atom_array)

    def test_extract_atom_array_error(self):
        """Test to see that the expected error is raised if the molecule
        contains a number of atom arrays not equal to 1.
        """
        with self.assertRaises(ValueError) as err:
            cml_reading_util.extractAtomArray(self.bad_water_1)
        expected_args = "Incorrect number of atom arrays."
        actual_args = err.exception.args[0]
        self.assertEqual(actual_args, expected_args)

    def test_extract_bond_array(self):
        """Test to see that the correct bond array is extracted.
        """
        expected_bond_array = (
            '<cml:bondArray xmlns:cml="http://www.xml-cm'
            + 'l.org/schema" xmlns:ssip="http://www-hunter.ch.cam.ac.uk/SSIP">\n  <cml:bond cml:atomRefs2="a1 a2'
            + '" cml:order="1"/>\n  <cml:bond cml:atomRefs2="a1 a3"'
            + ' cml:order="1"/>\n </cml:bondArray>\n'
        )
        actual_bond_array = etree.tounicode(self.bond_array)
        self.assertEqual(actual_bond_array, expected_bond_array)

    def test_extract_bond_array_error(self):
        """Test to see that the expected error is raised if the molecule
        contains a number of atom arrays not equal to 1.
        """
        with self.assertRaises(ValueError) as err:
            cml_reading_util.extractBondArray(self.bad_water_1)
        expected_args = "Incorrect number of bond arrays."
        actual_args = err.exception.args[0]
        self.assertEqual(actual_args, expected_args)

    def test_atom_cml(self):
        """Test to see that the expected atom CML is being used for future
        tests.
        """
        expected_atom = (
            u'<cml:atom xmlns:cml="http://www.xml-cml.org/schem'
            + 'a" cml:elementType="O" cml:id="a1" cml:x3="-0.198000" cml:y3="0.00'
            + '0000" cml:z3="0.347000"/>\n  '
        )
        actual_atom = etree.tounicode(self.oxygen_atom)
        self.assertEqual(actual_atom, expected_atom)

    def test_get_atom_element_type(self):
        """Test to see that the expected elementType is being returned.
        Also contains tests of error handling, using the bad cml.
        """
        expected_atom_element_type = "O"
        actual_atom_element_type = cml_reading_util.getAtomElementType(self.oxygen_atom)
        self.assertEqual(actual_atom_element_type, expected_atom_element_type)

        # now to test error handling
        with self.assertRaises(KeyError) as err:
            cml_reading_util.getAtomElementType(self.bad_oxygen)
        expected_err_args = "no elementType"
        actual_err_args = err.exception.args[0]
        self.assertEqual(actual_err_args, expected_err_args)

    def test_get_atom_x3(self):
        """Test to see that the expected x3 coordinate is being returned.
        """
        expected_atom_x3 = "-0.198000"
        actual_atom_x3 = cml_reading_util.getAtomX3(self.oxygen_atom)
        self.assertEqual(actual_atom_x3, expected_atom_x3)
        # now to test error handling
        with self.assertRaises(KeyError) as err:
            cml_reading_util.getAtomX3(self.bad_oxygen)
        expected_err_args = "no x3"
        actual_err_args = err.exception.args[0]
        self.assertEqual(actual_err_args, expected_err_args)

    def test_get_atom_y3(self):
        """Test to see that the expected y3 coordinate is being returned.
        """
        expected_atom_y3 = "0.000000"
        actual_atom_y3 = cml_reading_util.getAtomY3(self.oxygen_atom)
        self.assertEqual(actual_atom_y3, expected_atom_y3)
        # now to test error handling
        with self.assertRaises(KeyError) as err:
            cml_reading_util.getAtomY3(self.bad_oxygen)
        expected_err_args = "no y3"
        actual_err_args = err.exception.args[0]
        self.assertEqual(actual_err_args, expected_err_args)

    def test_get_atom_z3(self):
        """Test to see that the expected z3 coordinate is being returned.
        """
        expected_atom_z3 = "0.347000"
        actual_atom_z3 = cml_reading_util.getAtomZ3(self.oxygen_atom)
        self.assertEqual(actual_atom_z3, expected_atom_z3)
        # now to test error handling
        with self.assertRaises(KeyError) as err:
            cml_reading_util.getAtomZ3(self.bad_oxygen)
        expected_err_args = "no z3"
        actual_err_args = err.exception.args[0]
        self.assertEqual(actual_err_args, expected_err_args)

    def test_get_atom_isotope_number(self):
        """Test to see that the expected isotopeNumber is returned.
        """
        expected_atom_isotope_number = "16"
        self.oxygen_atom.set("{http://www.xml-cml.org/schema}isotopeNumber", "16")
        actual_atom_isotope_number = cml_reading_util.getAtomIsotopeNumber(
            self.oxygen_atom
        )
        self.assertEqual(actual_atom_isotope_number, expected_atom_isotope_number)
        # now to test error handling
        with self.assertRaises(KeyError) as err:
            cml_reading_util.getAtomIsotopeNumber(self.bad_oxygen)
        expected_err_args = "no isotopeNumber"
        actual_err_args = err.exception.args[0]
        self.assertEqual(actual_err_args, expected_err_args)

    def test_get_atom_attributes(self):
        """Test to see that the correct selection of atom attributes is
        returned.
        """
        expected_attributes = {
            "element": "O",
            "x3": "-0.198000",
            "y3": "0.000000",
            "z3": "0.347000",
        }
        actual_attributes = cml_reading_util.getAtomAttributes(self.oxygen_atom)
        self.assertDictEqual(actual_attributes, expected_attributes)
        # test with isotope number
        expected_attributes["isotope_number"] = "16"
        self.oxygen_atom.set("{http://www.xml-cml.org/schema}isotopeNumber", "16")
        actual_attributes = cml_reading_util.getAtomAttributes(self.oxygen_atom)
        self.assertDictEqual(actual_attributes, expected_attributes)
        # now to test error handling
        with self.assertRaises(KeyError) as err:
            cml_reading_util.getAtomAttributes(self.bad_oxygen)
        expected_err_args = "no elementType"
        actual_err_args = err.exception.args[0]
        self.assertEqual(actual_err_args, expected_err_args)

    def test_get_atom_list(self):
        """Test to see the expected atom list is returned.
        """
        expected_atom_list = [
            u'<cml:atom xmlns:cml="http://www.xml-cml.org/schem'
            + 'a" xmlns:ssip="http://www-hunter.ch.cam.ac.uk/SSIP" cml:elementType="O" cml:id="a1" cml:x3="-0.198000" cml:y3="0.00'
            + '0000" cml:z3="0.347000"/>\n  ',
            u'<cml:atom xmlns:cml="http://www.xml-cml.org/s'
            + 'chema" xmlns:ssip="http://www-hunter.ch.cam.ac.uk/SSIP" cml:elementType="H" cml:id="a2" cml:x3="0.760000" '
            + 'cml:y3="0.000000" cml:z3="0.204000"/>\n  ',
            u'<cml:atom xmlns:cml="http://www.xml-cml.org/'
            + 'schema" xmlns:ssip="http://www-hunter.ch.cam.ac.uk/SSIP" cml:elementType="H" cml:id="a3" cml:x3="-0.561000'
            + '" cml:y3="0.000000" cml:z3="-0.551000"/>\n ',
        ]
        actual_atom_list = [etree.tounicode(atom_cml) for atom_cml in self.atom_list]
        self.assertEqual(actual_atom_list, expected_atom_list)

    def test_get_atom_list_attributes(self):
        """Test to see if the expected dictionary of atom attribute
        dictionaries is returned.
        """
        expected_atom_list_attributes = {
            0: {"element": "O", "x3": "-0.198000", "y3": "0.000000", "z3": "0.347000"},
            1: {"element": "H", "x3": "0.760000", "y3": "0.000000", "z3": "0.204000"},
            2: {"element": "H", "x3": "-0.561000", "y3": "0.000000", "z3": "-0.551000"},
        }
        actual_atom_list_attributes = cml_reading_util.getAtomListAttributes(
            self.atom_list
        )
        self.assertEqual(actual_atom_list_attributes, expected_atom_list_attributes)

    def test_conv_att_to_atm_line_list(self):
        """test to see if expected atom line list is generated form the
        attribute dict.
        """
        expected_list = [
            "    O       -0.198000        0.000000        0.347000\n",
            "    H        0.760000        0.000000        0.204000\n",
            "    H       -0.561000        0.000000       -0.551000\n",
        ]
        actual_list = cml_reading_util.convertAtomListAttribToAtomLineList(
            self.atom_attrib
        )
        self.assertEqual(actual_list, expected_list)

    def test_atom_line_list_from_array(self):
        """test to see if expected atom line list is generated from the array.
        """
        expected_list = [
            "    O       -0.198000        0.000000        0.347000\n",
            "    H        0.760000        0.000000        0.204000\n",
            "    H       -0.561000        0.000000       -0.551000\n",
        ]
        actual_list = cml_reading_util.createAtomLineListFromArray(self.atom_array)
        self.assertEqual(actual_list, expected_list)

    def test_get_element_set_from_dict(self):
        """test to see if expected element set is generated from the attribute
        dictionary.
        """
        expected_set = {"O", "H"}
        actual_set = cml_reading_util.getElementSetFromAttributeDict(self.atom_attrib)
        self.assertEqual(actual_set, expected_set)

    def test_get_element_set_from_array(self):
        """test to see if expected element set is generated from the atom array.
        """
        expected_set = {"O", "H"}
        actual_set = cml_reading_util.createElementSetFromArray(self.atom_array)
        self.assertEqual(actual_set, expected_set)

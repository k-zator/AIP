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
Created on Sat Mar 26 16:28:58 2016

@author: mark
"""


import unittest
import logging
import os
import copy
import pathlib
from lxml import etree
import numpy as np
import nwchemcmlutils.nwUtils.cml_reading_util as cml_reading_util
import nwchemcmlutils.nwUtils.cml_writing_util as cml_writing_util


logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.WARN)


class CMLWritingUtilTestCase(unittest.TestCase):
    """Test class for the cml_reading_util script.
    """

    def setUp(self):
        """setUp for the class.
        """
        self.maxDiff = None
        # we read in the file with the example molecule in.
        self.parent_directory = pathlib.Path(__file__).parents[1]
        water_cml_file = etree.parse(
            (self.parent_directory / "resources/examplecml.xml").absolute().as_posix()
        )
        self.cml_namespace = {"cml": "http://www.xml-cml.org/schema"}
        water_cml = water_cml_file.findall(
            "cml:molecule", namespaces=self.cml_namespace
        )[0]
        self.water_cml = water_cml
        self.atom_array = cml_reading_util.extractAtomArray(self.water_cml)
        self.bond_array = cml_reading_util.extractBondArray(self.water_cml)
        self.atom_list = cml_reading_util.createAtomList(self.atom_array)
        self.oxygen_atom = copy.deepcopy(self.atom_list[0])
        self.new_coordinates = {
            0: {
                "coordinates": np.array(
                    [-1.59116481e-02, -1.23028388e-29, 3.70600788e-01]
                ),
                "element": "O",
            },
            1: {
                "coordinates": np.array(
                    [9.13545825e-01, 1.38240500e-29, 9.75447324e-02]
                ),
                "element": "H",
            },
            2: {
                "coordinates": np.array(
                    [-5.00634177e-01, 1.05052504e-30, -4.68145521e-01]
                ),
                "element": "H",
            },
        }

    def tearDown(self):
        """tearDown for the class.
        """
        del self.cml_namespace
        del self.water_cml
        del self.atom_array
        del self.bond_array
        del self.atom_list
        del self.oxygen_atom
        del self.new_coordinates
        if os.path.isfile("water_test_output.cml"):
            os.remove("water_test_output.cml")

    def test_set_atom_x3(self):
        """test to see if x3 attribute is correctly set for the oxygen atom.
        """
        expected_atom = (
            u'<cml:atom xmlns:cml="http://www.xml-cml.org/schem'
            + 'a" cml:elementType="O" cml:id="a1" cml:x3="-1.59116481e-02" cml:y3'
            + '="0.000000" cml:z3="0.347000"/>\n  '
        )
        new_x3_coord = "-1.59116481e-02"
        cml_writing_util.setAtomX3(new_x3_coord, self.oxygen_atom)
        self.assertEqual(etree.tounicode(self.oxygen_atom), expected_atom)

    def test_set_atom_y3(self):
        """test to see if y3 attribute is correctly set for the oxygen atom.
        """
        expected_atom = (
            u'<cml:atom xmlns:cml="http://www.xml-cml.org/schem'
            + 'a" cml:elementType="O" cml:id="a1" cml:x3="-0.198000" cml:y3'
            + '="-1.23028388e-29" cml:z3="0.347000"/>\n  '
        )
        new_y3_coord = "-1.23028388e-29"
        cml_writing_util.setAtomY3(new_y3_coord, self.oxygen_atom)
        self.assertEqual(etree.tounicode(self.oxygen_atom), expected_atom)

    def test_set_atom_z3(self):
        """test to see if z3 attribute is correctly set for the oxygen atom.
        """
        expected_atom = (
            u'<cml:atom xmlns:cml="http://www.xml-cml.org/schem'
            + 'a" cml:elementType="O" cml:id="a1" cml:x3="-0.198000" cml:y3'
            + '="0.000000" cml:z3="3.70600788e-01"/>\n  '
        )
        new_z3_coord = "3.70600788e-01"
        cml_writing_util.setAtomZ3(new_z3_coord, self.oxygen_atom)
        self.assertEqual(etree.tounicode(self.oxygen_atom), expected_atom)

    def test_set_atom_coordinates(self):
        """test to see if all 3 coordinates are changed.
        """
        expected_atom = (
            u'<cml:atom xmlns:cml="http://www.xml-cml.org/schem'
            + 'a" cml:elementType="O" cml:id="a1" cml:x3="-0.0159116481" cml:y3'
            + '="-1.23028388e-29" cml:z3="0.370600788"/>\n  '
        )
        new_coord_array = np.array([-1.59116481e-02, -1.23028388e-29, 3.70600788e-01])
        cml_writing_util.setAtomCoordinates(new_coord_array, self.oxygen_atom)
        self.assertEqual(etree.tounicode(self.oxygen_atom), expected_atom)

    def test_set_coords_for_array(self):
        """test to see if all the coordinates in an atom array are correctly
        changed.
        """
        expected_atom_array = (
            u'<cml:atomArray xmlns:cml="http://www.xml-cm'
            + 'l.org/schema" xmlns:ssip="http://www-hunter.ch.cam.ac.uk/SSIP">\n  <cml:atom cml:elementType="O" '
            + 'cml:id="a1" cml:x3="-0.0159116481" cml:y3="-1.23028388e-'
            + '29" cml:z3="0.370600788"/>\n  <cml:atom cml:elementT'
            + 'ype="H" cml:id="a2" cml:x3="0.913545825" cml:y3="1.38240'
            + '5e-29" cml:z3="0.0975447324"/>\n  <cml:atom cml:elem'
            + 'entType="H" cml:id="a3" cml:x3="-0.500634177" cml:y3="1.'
            + '05052504e-30" cml:z3="-0.468145521"/>\n </cml:at'
            + "omArray>\n "
        )
        cml_writing_util.setAtomCoordinatesForAtomArray(
            self.new_coordinates, self.atom_array
        )
        self.assertEqual(etree.tounicode(self.atom_array), expected_atom_array)

    def test_molecule_cml_after_coord_change(self):
        """test to see if coordinate change is mirrored in the molecule itself.
        """
        expected_water_cml = u"""<cml:molecule xmlns:cml="http://www.xml-cml.org/schema" xmlns:ssip="http://www-hunter.ch.cam.ac.uk/SSIP" ssip:stdInChIKey="XLYOFNOQVPJJNP-UHFFFAOYSA-N" StuartId="13_01_20121_new" cml:id="XLYOFNOQVPJJNP-UHFFFAOYSA-N">
 <cml:atomArray>
  <cml:atom cml:elementType="O" cml:id="a1" cml:x3="-0.0159116481" cml:y3="-1.23028388e-29" cml:z3="0.370600788"/>
  <cml:atom cml:elementType="H" cml:id="a2" cml:x3="0.913545825" cml:y3="1.382405e-29" cml:z3="0.0975447324"/>
  <cml:atom cml:elementType="H" cml:id="a3" cml:x3="-0.500634177" cml:y3="1.05052504e-30" cml:z3="-0.468145521"/>
 </cml:atomArray>
 <cml:bondArray>
  <cml:bond cml:atomRefs2="a1 a2" cml:order="1"/>
  <cml:bond cml:atomRefs2="a1 a3" cml:order="1"/>
 </cml:bondArray>
</cml:molecule>
"""
        cml_writing_util.setAtomCoordinatesForAtomArray(
            self.new_coordinates, self.atom_array
        )
        self.assertEqual(etree.tounicode(self.water_cml), expected_water_cml)

    def test_write_cml_to_file(self):
        """test to see if correct output file is produced.
        """
        expected_filename = (
            (self.parent_directory / "resources/water_out_example.cml")
            .absolute()
            .as_posix()
        )
        actual_filename = "water_test_output.cml"
        cml_writing_util.writeMoleculeCMLToFile(self.water_cml, actual_filename)
        with open(expected_filename, "r") as expected_file:
            with open(actual_filename, "r") as actual_file:
                actual_file_read_in = actual_file.read()
                expected_file_read_in = expected_file.read()
                self.assertEqual(actual_file_read_in, expected_file_read_in)

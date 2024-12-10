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
Script for tests of CMLFileReader functions.

@author: mark
"""

import unittest
import logging
import pathlib
import nwchemcmlutils.calcCreator.CMLFileReader as CMLFileReader

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.WARN)


class CMLFileReaderTestCase(unittest.TestCase):
    """Test case for the CMLFileReader module.
    """

    def setUp(self):
        """Set up for tests.
        """
        self.parent_directory = pathlib.Path(__file__).parents[1]
        self.example_cml_filename = (
            (self.parent_directory / "resources/examplecml.cml").absolute().as_posix()
        )
        self.example_cml_file = CMLFileReader.readCMLFile(self.example_cml_filename)
        self.molecule_cml_list = CMLFileReader.extractMoleculeList(
            self.example_cml_file
        )
        self.molecule_cml = self.molecule_cml_list[0]
        self.molecule_cml_string = CMLFileReader.writeMoleculeCMLForNwinFile(
            self.molecule_cml
        )

    def tearDown(self):
        """Tear down afterr test.
        """
        del self.example_cml_file
        del self.molecule_cml_list
        del self.molecule_cml
        del self.molecule_cml_string

    def test_len_molecule_list(self):
        """Test to check how many molecules are in list.
        """
        expected_length = 1
        actual_length = len(self.molecule_cml_list)
        self.assertEqual(actual_length, expected_length)

    def test_extract_stdinchikey(self):
        """Test to see StdInChiKey extracted matches expected value.
        """
        expected_value = "XLYOFNOQVPJJNP-UHFFFAOYSA-N"
        actual_value = CMLFileReader.extractStdInChiKey(self.molecule_cml)
        self.assertEqual(actual_value, expected_value)

    def test_write_molecule_cml_for_nwin(self):
        """test to see if expected string is written.
        """
        self.maxDiff = None
        expected_string = u'<cml:molecule xmlns:cml="http://www.xml-cml.org/schema" xmlns:ssip="http://www-hunter.ch.cam.ac.uk/SSIP" ssip:stdInChIKey="XLYOFNOQVPJJNP-UHFFFAOYSA-N" ssip:StuartId="13_01_20121_new" cml:id="XLYOFNOQVPJJNP-UHFFFAOYSA-N">\n   <cml:atomArray>\n    <cml:atom cml:elementType="O" cml:id="a1" cml:x3="-0.198000" cml:y3="0.000000" cml:z3="0.347000"/>\n    <cml:atom cml:elementType="H" cml:id="a2" cml:x3="0.760000" cml:y3="0.000000" cml:z3="0.204000"/>\n    <cml:atom cml:elementType="H" cml:id="a3" cml:x3="-0.561000" cml:y3="0.000000" cml:z3="-0.551000"/>\n   </cml:atomArray>\n   <cml:bondArray>\n    <cml:bond cml:atomRefs2="a1 a2" cml:order="1"/>\n    <cml:bond cml:atomRefs2="a1 a3" cml:order="1"/>\n   </cml:bondArray>\n  </cml:molecule>\n  '
        actual_string = self.molecule_cml_string
        self.assertEqual(actual_string, expected_string)

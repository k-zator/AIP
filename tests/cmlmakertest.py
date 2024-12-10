# -*- coding: utf-8 -*-
#    cmlgenerator creates fully namespaced CML for molecules from input structures.
#    Copyright (C) 2019  Mark D. Driver
#
#    cmlgenerator is free software: you can redistribute it and/or modify
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
"""
Script for testing the cmlmaker methods.

@author: mark
"""

import logging
import unittest
import re
import os
import pathlib
from lxml import etree
import cmlgenerator.cmlmaker as cmlmaker

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.WARN)


class CMLMakerTestCase(unittest.TestCase):
    """Test case for cmlmaker.
    """

    def setUp(self):
        """Test set up.
        """
        self.maxDiff = None
        parent_directory = pathlib.Path(__file__).parents[0]
        self.expected_filename = (
            (parent_directory / "test_files" / "ethanolnamespaced.cml")
            .absolute()
            .as_posix()
        )
        self.mol_filename = (
            (parent_directory / "test_files" / "ethanol.mol").absolute().as_posix()
        )
        self.inchikey = "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"
        self.expected_cml = """<cml:molecule xmlns:ssip="http://www-hunter.ch.cam.ac.uk/SSIP" xmlns:cml="http://www.xml-cml.org/schema" cml:id="LFQSCWFLJHTTHZ-UHFFFAOYSA-N" ssip:stdInChIKey="LFQSCWFLJHTTHZ-UHFFFAOYSA-N">
  <cml:atomArray>
    <cml:atom cml:id="a1" cml:elementType="C" cml:hydrogenCount="3" cml:x3="-0.888300" cml:y3="0.167000" cml:z3="-0.027300"/>
    <cml:atom cml:id="a2" cml:elementType="C" cml:hydrogenCount="2" cml:x3="0.465800" cml:y3="-0.511600" cml:z3="-0.036800"/>
    <cml:atom cml:id="a3" cml:elementType="O" cml:hydrogenCount="1" cml:x3="1.431100" cml:y3="0.322900" cml:z3="0.586700"/>
    <cml:atom cml:id="a4" cml:elementType="H" cml:hydrogenCount="0" cml:x3="-0.848700" cml:y3="1.117500" cml:z3="-0.569500"/>
    <cml:atom cml:id="a5" cml:elementType="H" cml:hydrogenCount="0" cml:x3="-1.647100" cml:y3="-0.470400" cml:z3="-0.489600"/>
    <cml:atom cml:id="a6" cml:elementType="H" cml:hydrogenCount="0" cml:x3="-1.196400" cml:y3="0.397800" cml:z3="0.997700"/>
    <cml:atom cml:id="a7" cml:elementType="H" cml:hydrogenCount="0" cml:x3="0.792000" cml:y3="-0.722400" cml:z3="-1.059700"/>
    <cml:atom cml:id="a8" cml:elementType="H" cml:hydrogenCount="0" cml:x3="0.424600" cml:y3="-1.455900" cml:z3="0.513800"/>
    <cml:atom cml:id="a9" cml:elementType="H" cml:hydrogenCount="0" cml:x3="1.467100" cml:y3="1.155000" cml:z3="0.084800"/>
  </cml:atomArray>
  <cml:bondArray>
    <cml:bond cml:atomRefs2="a1 a2" cml:order="1"/>
    <cml:bond cml:atomRefs2="a2 a3" cml:order="1"/>
    <cml:bond cml:atomRefs2="a1 a4" cml:order="1"/>
    <cml:bond cml:atomRefs2="a1 a5" cml:order="1"/>
    <cml:bond cml:atomRefs2="a1 a6" cml:order="1"/>
    <cml:bond cml:atomRefs2="a2 a7" cml:order="1"/>
    <cml:bond cml:atomRefs2="a2 a8" cml:order="1"/>
    <cml:bond cml:atomRefs2="a3 a9" cml:order="1"/>
  </cml:bondArray>
</cml:molecule>
"""
        self.expected_cml_regex = """<molecule>
 <atomArray>
  <atom id="a1" elementType="C" x3="(-)?([0-9][.])([0-9]){6}" y3="(-)?([0-9][.])([0-9]){6}" z3="(-)?([0-9][.])([0-9]){6}"/>
  <atom id="a2" elementType="C" x3="(-)?([0-9][.])([0-9]){6}" y3="(-)?([0-9][.])([0-9]){6}" z3="(-)?([0-9][.])([0-9]){6}"/>
  <atom id="a3" elementType="O" x3="(-)?([0-9][.])([0-9]){6}" y3="(-)?([0-9][.])([0-9]){6}" z3="(-)?([0-9][.])([0-9]){6}"/>
  <atom id="a4" elementType="H" x3="(-)?([0-9][.])([0-9]){6}" y3="(-)?([0-9][.])([0-9]){6}" z3="(-)?([0-9][.])([0-9]){6}"/>
  <atom id="a5" elementType="H" x3="(-)?([0-9][.])([0-9]){6}" y3="(-)?([0-9][.])([0-9]){6}" z3="(-)?([0-9][.])([0-9]){6}"/>
  <atom id="a6" elementType="H" x3="(-)?([0-9][.])([0-9]){6}" y3="(-)?([0-9][.])([0-9]){6}" z3="(-)?([0-9][.])([0-9]){6}"/>
  <atom id="a7" elementType="H" x3="(-)?([0-9][.])([0-9]){6}" y3="(-)?([0-9][.])([0-9]){6}" z3="(-)?([0-9][.])([0-9]){6}"/>
  <atom id="a8" elementType="H" x3="(-)?([0-9][.])([0-9]){6}" y3="(-)?([0-9][.])([0-9]){6}" z3="(-)?([0-9][.])([0-9]){6}"/>
  <atom id="a9" elementType="H" x3="(-)?([0-9][.])([0-9]){6}" y3="(-)?([0-9][.])([0-9]){6}" z3="(-)?([0-9][.])([0-9]){6}"/>
 </atomArray>
 <bondArray>
  <bond atomRefs2="a1 a2" order="1"/>
  <bond atomRefs2="a2 a3" order="1"/>
  <bond atomRefs2="a1 a4" order="1"/>
  <bond atomRefs2="a1 a5" order="1"/>
  <bond atomRefs2="a1 a6" order="1"/>
  <bond atomRefs2="a2 a7" order="1"/>
  <bond atomRefs2="a2 a8" order="1"/>
  <bond atomRefs2="a3 a9" order="1"/>
 </bondArray>
</molecule>
"""

    def tearDown(self):
        """clean up after tests.
        """
        del self.expected_cml
        del self.expected_cml_regex
        del self.expected_filename
        if os.path.isfile(self.inchikey + ".cml"):
            os.remove(self.inchikey + ".cml")
        del self.inchikey

    def test_generate_cml_file_from_smiles(self):
        """Test to see if expected file is created from SMILES.
        """
        actual_filename = self.inchikey + ".cml"
        cmlmaker.generate_cml_file_from_smiles("CCO")
        with open(actual_filename, "r") as actual_file:
            actual_contents = actual_file.read()
            self.assertTrue(
                type(re.match(self.expected_cml_regex, actual_contents)) is not None
            )

    def test_generate_cml_file_from_file(self):
        """Test to see if cml file is generated from input file.
        """
        actual_filename = self.inchikey + ".cml"
        cmlmaker.generate_cml_file_from_file(self.mol_filename, "mol")
        with open(actual_filename, "r") as actual_file:
            actual_contents = actual_file.read()
            with open(self.expected_filename, "r") as expected_file:
                expected_contents = expected_file.read()
                self.assertMultiLineEqual(expected_contents, actual_contents)

    def test_create_cml_from_smiles(self):
        """Test to see if expected molecule element is produced.
        """
        actual_inchikey, actual_cml = cmlmaker.create_cml_from_smiles("CCO")
        self.assertEqual("LFQSCWFLJHTTHZ-UHFFFAOYSA-N", actual_inchikey)
        self.assertTrue(
            type(
                re.match(
                    self.expected_cml_regex,
                    etree.tounicode(actual_cml, pretty_print=True),
                )
            )
            is not None
        )

    def test_create_cml_from_file(self):
        """Test to see if cml is generated from file.
        """
        actual_inchikey, actual_cml = cmlmaker.create_cml_from_file(
            self.mol_filename, "mol"
        )
        self.assertMultiLineEqual(
            self.expected_cml, etree.tounicode(actual_cml, pretty_print=True)
        )

    def test_create_filename_frominchikey(self):
        """Test to see if expected file name is returned.
        """
        self.assertEqual(
            "LFQSCWFLJHTTHZ-UHFFFAOYSA-N.cml",
            cmlmaker.create_filename_frominchikey(self.inchikey),
        )

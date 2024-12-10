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
Created on Tue Jul 31 14:17:42 2018

@author: mark
"""

import logging
import unittest
import pathlib
import os
from lxml import etree
import cmlgenerator.cmlnamespacing as cmlname

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.WARN)


class CMLNamespacingTestCase(unittest.TestCase):
    """Test case for smilestocml methods
    """

    def setUp(self):
        """Set up before tests.
        """
        self.maxDiff = None
        parent_directory = pathlib.Path(__file__).parents[0]
        self.expected_filename = (
            (parent_directory / "test_files" / "ethanolnamespaced.cml")
            .absolute()
            .as_posix()
        )
        self.unnamespaced_filename = (
            (parent_directory / "test_files" / "ethanolnonamespace.cml")
            .absolute()
            .as_posix()
        )
        self.unnamespaced_filename2 = (
            (parent_directory / "test_files" / "LFQSCWFLJHTTHZ-UHFFFAOYSA-N.cml")
            .absolute()
            .as_posix()
        )
        self.inchikey = "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"
        self.input_cml = """<molecule>
 <atomArray>
  <atom id="a1" elementType="C" hydrogenCount="3" x3="-0.888300" y3="0.167000" z3="-0.027300"/>
  <atom id="a2" elementType="C" hydrogenCount="2" x3="0.465800" y3="-0.511600" z3="-0.036800"/>
  <atom id="a3" elementType="O" hydrogenCount="1" x3="1.431100" y3="0.322900" z3="0.586700"/>
  <atom id="a4" elementType="H" hydrogenCount="0" x3="-0.848700" y3="1.117500" z3="-0.569500"/>
  <atom id="a5" elementType="H" hydrogenCount="0" x3="-1.647100" y3="-0.470400" z3="-0.489600"/>
  <atom id="a6" elementType="H" hydrogenCount="0" x3="-1.196400" y3="0.397800" z3="0.997700"/>
  <atom id="a7" elementType="H" hydrogenCount="0" x3="0.792000" y3="-0.722400" z3="-1.059700"/>
  <atom id="a8" elementType="H" hydrogenCount="0" x3="0.424600" y3="-1.455900" z3="0.513800"/>
  <atom id="a9" elementType="H" hydrogenCount="0" x3="1.467100" y3="1.155000" z3="0.084800"/>
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
        self.molecule_element = cmlname.read_molecule_string(self.input_cml)
        self.atom_array_element = self.molecule_element.getchildren()[0]
        self.atom_element = self.atom_array_element.getchildren()[0]
        self.bond_array_element = self.molecule_element.getchildren()[1]
        self.bond_element = self.bond_array_element.getchildren()[0]
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
    def tearDown(self):
        """Clean up after tests
        """
        if os.path.isfile(self.inchikey + ".cml"):
            os.remove(self.inchikey + ".cml")
        del self.inchikey
        del self.input_cml
        del self.molecule_element
        del self.atom_array_element
        del self.atom_element
        del self.bond_array_element
        del self.bond_element

    def test_convert_unnamespaced_file(self):
        """Test expected conversion is carried out.

        Returns
        -------
        None.

        """
        actual_filename = self.inchikey + ".cml"
        cmlname.convert_unnamespaced_file(self.unnamespaced_filename2, "")
        with open(actual_filename, "r") as actual_file:
            actual_contents = actual_file.read()
            with open(self.expected_filename, "r") as expected_file:
                expected_contents = expected_file.read()
                self.assertMultiLineEqual(expected_contents, actual_contents)
        
    def test_convert_molecule_from_file(self):
        """Test expected CML is produced.

        Returns
        -------
        None.

        """
        actual_cml = cmlname.convert_molecule_from_file(self.inchikey, self.unnamespaced_filename)
        self.assertMultiLineEqual(
            self.expected_cml, etree.tounicode(actual_cml, pretty_print=True)
        )
    
    def test_write_cml_to_file(self):
        """Test expected contents are written.

        Returns
        -------
        None.

        """
        actual_filename = self.inchikey + ".cml"
        cmlname.write_cml_to_file(etree.XML(self.expected_cml), actual_filename)
        with open(actual_filename, "r") as actual_file:
            actual_contents = actual_file.read()
            with open(self.expected_filename, "r") as expected_file:
                expected_contents = expected_file.read()
                self.assertMultiLineEqual(expected_contents, actual_contents)
    
    def test_convert_molecule_from_string(self):
        """Test to see if molecule is converted correctly from string.
        """
        
        actual_element = cmlname.convert_molecule_from_string(
            self.inchikey, self.input_cml
        )
        self.assertMultiLineEqual(
            self.expected_cml, etree.tounicode(actual_element, pretty_print=True)
        )

    def test_convert_molecule(self):
        """Test to see if molecule is converted
        """
        expected_element = """<cml:molecule xmlns:ssip="http://www-hunter.ch.cam.ac.uk/SSIP" xmlns:cml="http://www.xml-cml.org/schema" cml:id="LFQSCWFLJHTTHZ-UHFFFAOYSA-N" ssip:stdInChIKey="LFQSCWFLJHTTHZ-UHFFFAOYSA-N">
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
        actual_element = cmlname.convert_molecule(self.inchikey, self.molecule_element)
        self.assertMultiLineEqual(
            expected_element, etree.tounicode(actual_element, pretty_print=True)
        )

    def test_append_inchikey(self):
        """Test to see if inchikey is appended
        """
        expected_element = """<molecule xmlns:ns0="http://www.xml-cml.org/schema" xmlns:ns1="http://www-hunter.ch.cam.ac.uk/SSIP" ns0:id="LFQSCWFLJHTTHZ-UHFFFAOYSA-N" ns1:stdInChIKey="LFQSCWFLJHTTHZ-UHFFFAOYSA-N">
 <atomArray>
  <atom id="a1" elementType="C" hydrogenCount="3" x3="-0.888300" y3="0.167000" z3="-0.027300"/>
  <atom id="a2" elementType="C" hydrogenCount="2" x3="0.465800" y3="-0.511600" z3="-0.036800"/>
  <atom id="a3" elementType="O" hydrogenCount="1" x3="1.431100" y3="0.322900" z3="0.586700"/>
  <atom id="a4" elementType="H" hydrogenCount="0" x3="-0.848700" y3="1.117500" z3="-0.569500"/>
  <atom id="a5" elementType="H" hydrogenCount="0" x3="-1.647100" y3="-0.470400" z3="-0.489600"/>
  <atom id="a6" elementType="H" hydrogenCount="0" x3="-1.196400" y3="0.397800" z3="0.997700"/>
  <atom id="a7" elementType="H" hydrogenCount="0" x3="0.792000" y3="-0.722400" z3="-1.059700"/>
  <atom id="a8" elementType="H" hydrogenCount="0" x3="0.424600" y3="-1.455900" z3="0.513800"/>
  <atom id="a9" elementType="H" hydrogenCount="0" x3="1.467100" y3="1.155000" z3="0.084800"/>
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
        cmlname.append_inchikey(self.molecule_element, self.inchikey)
        self.assertMultiLineEqual(
            expected_element, etree.tounicode(self.molecule_element, pretty_print=True)
        )

    def test_namespace_molecule(self):
        """Test to see if molecule is namespaced.
        """
        expected_element = """<cml:molecule xmlns:ssip="http://www-hunter.ch.cam.ac.uk/SSIP" xmlns:cml="http://www.xml-cml.org/schema">
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
        actual_element = cmlname.namespace_molecule(self.molecule_element)
        self.assertMultiLineEqual(
            expected_element, etree.tounicode(actual_element, pretty_print=True)
        )

    def test_namespace_bondarray(self):
        """Test to see if bond array is namespaced.
        """
        expected_element = """<cml:bondArray xmlns:ssip="http://www-hunter.ch.cam.ac.uk/SSIP" xmlns:cml="http://www.xml-cml.org/schema">
  <cml:bond cml:atomRefs2="a1 a2" cml:order="1"/>
  <cml:bond cml:atomRefs2="a2 a3" cml:order="1"/>
  <cml:bond cml:atomRefs2="a1 a4" cml:order="1"/>
  <cml:bond cml:atomRefs2="a1 a5" cml:order="1"/>
  <cml:bond cml:atomRefs2="a1 a6" cml:order="1"/>
  <cml:bond cml:atomRefs2="a2 a7" cml:order="1"/>
  <cml:bond cml:atomRefs2="a2 a8" cml:order="1"/>
  <cml:bond cml:atomRefs2="a3 a9" cml:order="1"/>
</cml:bondArray>
"""
        actual_element = cmlname.namespace_bondarray(self.bond_array_element)
        self.assertMultiLineEqual(
            expected_element, etree.tounicode(actual_element, pretty_print=True)
        )

    def test_namespace_bond(self):
        """Test to see if bond is namespaced.
        """
        expected_element = """<cml:bond xmlns:ssip="http://www-hunter.ch.cam.ac.uk/SSIP" xmlns:cml="http://www.xml-cml.org/schema" cml:atomRefs2="a1 a2" cml:order="1"/>
"""
        actual_element = cmlname.namespace_bond(self.bond_element)
        self.assertMultiLineEqual(
            expected_element, etree.tounicode(actual_element, pretty_print=True)
        )

    def test_namespace_atomarray(self):
        """Test to see if atom array is namespaced.
        """
        expected_element = """<cml:atomArray xmlns:ssip="http://www-hunter.ch.cam.ac.uk/SSIP" xmlns:cml="http://www.xml-cml.org/schema">
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
"""
        actual_element = cmlname.namespace_atomarray(self.atom_array_element)
        self.assertMultiLineEqual(
            expected_element, etree.tounicode(actual_element, pretty_print=True)
        )

    def test_namespace_atom(self):
        """Test to see if expected element is returned.
        """
        expected_element = """<cml:atom xmlns:ssip="http://www-hunter.ch.cam.ac.uk/SSIP" xmlns:cml="http://www.xml-cml.org/schema" cml:id="a1" cml:elementType="C" cml:hydrogenCount="3" cml:x3="-0.888300" cml:y3="0.167000" cml:z3="-0.027300"/>
"""
        actual_element = cmlname.namespace_atom(self.atom_element)
        self.assertMultiLineEqual(
            expected_element, etree.tounicode(actual_element, pretty_print=True)
        )

    def test_read_molecule_string(self):
        """Test to see if molecule is read in.
        """
        self.assertMultiLineEqual(
            self.input_cml, etree.tounicode(self.molecule_element, pretty_print=True)
        )

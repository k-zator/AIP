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
Script for tests of structuretocml methods.

@author: mark
"""
import logging
import unittest
import pathlib
import cmlgenerator.structuretocml as structcml

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.WARN)


class StructureToCMLTestCase(unittest.TestCase):
    """Test case for structure to cml script methods.
    """

    def setUp(self):
        """Set up before tests.
        """
        parent_directory = pathlib.Path(__file__).resolve().parents[0]
        self.mol_filename = (
            (parent_directory / "test_files/ethanol.mol").absolute().as_posix()
        )
        self.expected_mol = """
        RDKit          3D

  9  8  0  0  0  0  0  0  0  0999 V2000
   -0.8883    0.1670   -0.0273 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4658   -0.5116   -0.0368 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4311    0.3229    0.5867 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8487    1.1175   -0.5695 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6471   -0.4704   -0.4896 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1964    0.3978    0.9977 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.7920   -0.7224   -1.0597 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.4246   -1.4559    0.5138 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.4671    1.1550    0.0848 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  1  4  1  0
  1  5  1  0
  1  6  1  0
  2  7  1  0
  2  8  1  0
  3  9  1  0
M  END
"""
        self.expected_cml = """<molecule>
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
        self.expected_mol2 = """@<TRIPOS>MOLECULE
*****
 9 8 0 0 0
SMALL
GASTEIGER

@<TRIPOS>ATOM
      1 C          -0.8883    0.1670   -0.0273 C.3     1  UNL1       -0.0418
      2 C           0.4658   -0.5116   -0.0368 C.3     1  UNL1        0.0414
      3 O           1.4311    0.3229    0.5867 O.3     1  UNL1       -0.3953
      4 H          -0.8487    1.1175   -0.5695 H       1  UNL1        0.0252
      5 H          -1.6471   -0.4704   -0.4896 H       1  UNL1        0.0252
      6 H          -1.1964    0.3978    0.9977 H       1  UNL1        0.0252
      7 H           0.7920   -0.7224   -1.0597 H       1  UNL1        0.0554
      8 H           0.4246   -1.4559    0.5138 H       1  UNL1        0.0554
      9 H           1.4671    1.1550    0.0848 H       1  UNL1        0.2094
@<TRIPOS>BOND
     1     1     2    1
     2     2     3    1
     3     1     4    1
     4     1     5    1
     5     1     6    1
     6     2     7    1
     7     2     8    1
     8     3     9    1
"""

    def tearDown(self):
        """Clean up after tests.
        """
        del self.expected_cml
        del self.expected_mol
        del self.expected_mol2

    def test_generate_inchikey_for_cml(self):
        """Test to see if expected inchikey is returned.
        """
        inchikey = structcml.generate_inchikey_for_cml(self.expected_cml)
        self.assertEqual("LFQSCWFLJHTTHZ-UHFFFAOYSA-N", inchikey)

    def test_create_inchikey_conversion(self):
        """Test to see if conversion is created
        """
        conversion = structcml.create_inchikey_conversion("cml")
        self.assertTrue(conversion is not None)

    def test_perform_conversion_from_string(self):
        """Test to see if expected converion is carried out on string
        """
        actual_cml = structcml.perform_conversion_from_string(self.expected_mol, "mol")
        self.assertMultiLineEqual(self.expected_cml, actual_cml)

    def test_perform_conversion_from_file(self):
        """Test to see if conversion of file is carried out as expected.
        """
        actual_cml = structcml.perform_conversion_from_file(self.mol_filename, "mol")
        self.assertMultiLineEqual(self.expected_cml, actual_cml)

    def test_create_cml_conversion(self):
        """Test to see if conversion is created.
        """
        conversion = structcml.create_cml_conversion("mol", aromatic=False)
        self.assertTrue(conversion is not None)
    
    def test_create_mol2_conversion(self):
        """Test to see if conversion is created.
        """
        conversion = structcml.create_mol2_conversion("mol")
        self.assertTrue(conversion is not None)
    
    def test_perform_mol2_conversion_from_file(self):
        """Test to see if conversion of file is carried out as expected.
        """
        actual_mol2 = structcml.perform_mol2_conversion_from_file(self.mol_filename, "mol")
        self.assertMultiLineEqual(self.expected_mol2, actual_mol2)

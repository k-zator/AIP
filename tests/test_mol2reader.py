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
import filereader.mol2reader as reader
import os

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.WARN)

FIXTURE_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'test_files', '')

mol2_nitromethane = f"{FIXTURE_DIR}/nitromethane.mol2"


class TestMol2Reader(unittest.TestCase):
    """Test case for structure to cml script methods.
    """

    def setUp(self):
        """Set up before tests.
        """
        parent_directory = pathlib.Path(__file__).resolve().parents[0]
        self.mol2_filename = (
            (parent_directory / "test_files/ethanol.mol2").absolute().as_posix()
        )
        self.string_mol2 = """@<TRIPOS>MOLECULE
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
        self.initial_string = "@<TRIPOS>ATOM"
        self.line_0 ="      1 C          -0.8883    0.1670   -0.0273 C.3     1  UNL1       -0.0418"
        self.line_0_split =["1","C","-0.8883" ,"0.1670","-0.0273","C.3","1","UNL1","-0.0418"]
        self.mol2_from_string = reader.Mol2Reader(self.string_mol2, from_file=False)
        self.mol2_from_file = reader.Mol2Reader(self.mol2_filename, from_file=True)
        self.mol2_nitromethane = reader.Mol2Reader(mol2_nitromethane, from_file=True)

    def tearDown(self):
        """Clean up after tests.
        """
        del self.mol2_filename
        del self.string_mol2
        del self.initial_string
        del self.line_0
        del self.line_0_split
        del self.mol2_from_string
        del self.mol2_from_file


    def test_file_as_string(self):
        self.assertEqual(self.mol2_from_file._list_lines_, self.mol2_from_string._list_lines_)
    
    def test_main_lines_index(self):
        self.assertEqual(self.line_0, self.mol2_from_file._main_lines_[0])

    def test_mol2atom(self):
        line_0 = self.mol2_from_file._main_lines_[0]
        atom = reader.Mol2Atom(line_0)
        self.assertEqual(atom.element, "C")
        self.assertEqual(atom.sybyl, "C.3")

    def test_list_atoms(self):
        list_atoms = self.mol2_from_string.list_atoms
        atom = list_atoms[2]
        self.assertEqual(atom.element, "O")
        self.assertEqual(atom.sybyl, "O.3")

    def test_mol2_nitromethane(self):
        self.assertEqual(self.mol2_nitromethane.list_atoms[0].sybyl, "N.pl3")
        self.assertEqual(self.mol2_nitromethane.list_atoms[-1].sybyl, "H")

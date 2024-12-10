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
Created on Sun Apr 24 11:16:45 2016Script contains test case for Molecule class.

@author: mark
"""

import unittest
import logging
import numpy as np
from nwchemcmlutils.DataClasses.Cartesian3D import Cartesian3D
from nwchemcmlutils.DataClasses.DistanceUnits import DistanceUnits
from nwchemcmlutils.DataClasses.Atom import Atom
from nwchemcmlutils.DataClasses.Molecule import Molecule

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.WARN)


class MoleculeTestCase(unittest.TestCase):
    """Test case for the Atom class.
    """

    def setUp(self):
        """Set up for tests.
        """
        cartesian_3d_111_au = Cartesian3D([1, 1, 1], unit=DistanceUnits.atomic_units)
        self.atom_o_1 = Atom(8, cartesian_3d_111_au, id_num=0)
        cartesian_3d_111_ang = Cartesian3D([1, 1, 1], unit=DistanceUnits.angstroms)
        self.atom_o_2 = Atom(8, cartesian_3d_111_ang)
        cartesian_3d_100_au = Cartesian3D([1, 0, 0], unit=DistanceUnits.atomic_units)
        self.atom_h_1 = Atom(1, cartesian_3d_100_au, id_num=1)
        cartesian_3d_010_au = Cartesian3D([0, 1, 0], unit=DistanceUnits.atomic_units)
        self.atom_h_2 = Atom(1, cartesian_3d_010_au, id_num=2)
        self.mol_water_1 = Molecule([self.atom_o_1, self.atom_h_1, self.atom_h_2])
        atom_o_1_ang = Atom(
            8,
            Cartesian3D(
                [0.529177249, 0.529177249, 0.529177249], unit=DistanceUnits.angstroms
            ),
        )
        atom_h_1_ang = Atom(
            1, Cartesian3D([0.529177249, 0.0, 0.0], unit=DistanceUnits.angstroms)
        )
        atom_h_2_ang = Atom(
            1, Cartesian3D([0.0, 0.529177249, 0.0], unit=DistanceUnits.angstroms)
        )
        self.mol_water_2_ang = Molecule(
            [atom_o_1_ang, atom_h_1_ang, atom_h_2_ang], unit=DistanceUnits.angstroms
        )
        self.rot_matrix_45_100 = np.array(
            [[1, 0, 0], [0.0, 0.70710678, 0.70710678], [0.0, -0.70710678, 0.70710678]]
        )

    def tearDown(self):
        """Tear down of tests.
        """
        del self.atom_o_1
        del self.atom_o_2
        del self.atom_h_1
        del self.atom_h_2
        del self.mol_water_1
        del self.rot_matrix_45_100

    def test___str__(self):
        """Test to see if expected
        """
        actual_string = str(self.mol_water_1)
        self.assertEqual(382, len(actual_string))

    def test___repr__(self):
        """Test to see if expected repr is returned.
        """
        actual_repr = repr(self.mol_water_1)
        self.assertEqual(415, len(actual_repr))

    def test___eq__(self):
        """Test to see if == returns true.
        """
        expected_molecule = Molecule(
            set(
                [
                    Atom(
                        8,
                        Cartesian3D([1.0, 1.0, 1.0], unit=DistanceUnits.atomic_units),
                        mass=None,
                    ),
                    Atom(
                        1,
                        Cartesian3D([1.0, 0.0, 0.0], unit=DistanceUnits.atomic_units),
                        mass=None,
                    ),
                    Atom(
                        1,
                        Cartesian3D([0.0, 1.0, 0.0], unit=DistanceUnits.atomic_units),
                        mass=None,
                    ),
                ]
            ),
            unit=DistanceUnits.atomic_units,
        )
        actual_molecule = self.mol_water_1
        self.assertEqual(actual_molecule, expected_molecule)

    def test_add_atom_error(self):
        """Test to see expected error is raised when atom with wrong unit is
        added.
        """
        with self.assertRaises(ValueError):
            self.mol_water_1.addAtom(self.atom_o_2)

    def test_translate_atom(self):
        """Test to see if molecule is correctly translated.
        """
        expected_molecule = Molecule(
            set(
                [
                    Atom(
                        1,
                        Cartesian3D([0.0, 1.0, 1.0], unit=DistanceUnits.atomic_units),
                        mass=None,
                    ),
                    Atom(
                        8,
                        Cartesian3D([1.0, 1.0, 2.0], unit=DistanceUnits.atomic_units),
                        mass=None,
                    ),
                    Atom(
                        1,
                        Cartesian3D([1.0, 0.0, 1.0], unit=DistanceUnits.atomic_units),
                        mass=None,
                    ),
                ]
            ),
            unit=DistanceUnits.atomic_units,
        )
        translation_cartesian_3d = Cartesian3D(
            [0, 0, 1], unit=DistanceUnits.atomic_units
        )
        actual_molecule = self.mol_water_1.translateMolecule(translation_cartesian_3d)
        self.assertEqual(actual_molecule, expected_molecule)

    def test_translate_atom_error(self):
        """Test to see that correct error is raised if Cartesian3D with
        incorrect units is given as an argument. Also checks that AttributeError
        is raised if a Cartesian3D object is not used as the arguement.
        """
        with self.assertRaises(ValueError):
            self.mol_water_1.translateMolecule(
                Cartesian3D([1, 0, 0], unit=DistanceUnits.angstroms)
            )
        with self.assertRaises(AttributeError):
            self.mol_water_1.translateMolecule([1, 0, 0])

    def test_calculate_centroid(self):
        """Test to see if expected Cartesian3D is returned.
        """
        expected_cartesian_3d = Cartesian3D(
            [2.0 / 3.0, 2.0 / 3.0, 1 / 3.0], unit=DistanceUnits.atomic_units
        )
        actual_cartesian_3d = self.mol_water_1.calculateCentroidVector()
        self.assertEqual(actual_cartesian_3d, expected_cartesian_3d)

    def test_rotate_molecule(self):
        """Test to see if molecule is correctly rotated.
        """
        expected_molecule = Molecule(
            [
                Atom(
                    8,
                    Cartesian3D(
                        [1.0, 1.41421356, 0.0], unit=DistanceUnits.atomic_units
                    ),
                    mass=None,
                ),
                Atom(
                    1,
                    Cartesian3D([1.0, 0.0, 0.0], unit=DistanceUnits.atomic_units),
                    mass=None,
                ),
                Atom(
                    1,
                    Cartesian3D(
                        [0.0, 0.70710678, -0.70710678], unit=DistanceUnits.atomic_units
                    ),
                    mass=None,
                ),
            ],
            unit=DistanceUnits.atomic_units,
        )
        actual_molecule = self.mol_water_1.rotateMolecule(self.rot_matrix_45_100)
        self.assertEqual(actual_molecule, expected_molecule)

    def test_rotate_molecule_about_centroid(self):
        """Test to see if molecule is correctly rotated about the centroid.
        """
        expected_molecule = Molecule(
            [
                Atom(
                    8,
                    Cartesian3D(
                        [1.0, 1.37377345, 0.56903559], unit=DistanceUnits.atomic_units
                    ),
                    mass=None,
                ),
                Atom(
                    1,
                    Cartesian3D(
                        [0.0, 0.66666667, -0.13807119], unit=DistanceUnits.atomic_units
                    ),
                    mass=None,
                ),
                Atom(
                    1,
                    Cartesian3D(
                        [1.0, -0.04044011, 0.56903559], unit=DistanceUnits.atomic_units
                    ),
                    mass=None,
                ),
            ],
            unit=DistanceUnits.atomic_units,
        )
        actual_molecule = self.mol_water_1.rotateMoleculeAboutCentroid(
            self.rot_matrix_45_100
        )
        self.assertEqual(actual_molecule, expected_molecule)

    def test_convert_to_au_from_ang(self):
        """Test conversion to atomic units from angstroms.
        """
        expected_molecule = self.mol_water_1
        actual_molecule = self.mol_water_2_ang.convertMoleculeToAU()
        self.assertEqual(actual_molecule, expected_molecule)

    def test_convert_to_au_from_au(self):
        """Test conversion to atomic units from atomic units.
        """
        expected_molecule = self.mol_water_1
        actual_molecule = self.mol_water_1.convertMoleculeToAU()
        self.assertEqual(actual_molecule, expected_molecule)

    def test_convert_to_ang_from_au(self):
        """Test conversion to angstroms to atomic units.
        """
        expected_molecule = self.mol_water_2_ang
        actual_molecule = self.mol_water_1.convertMoleculeToAng()
        self.assertEqual(actual_molecule, expected_molecule)

    def test_convert_to_ang_from_ang(self):
        """Test conversion to angsttroms from angstroms.
        """
        expected_molecule = self.mol_water_2_ang
        actual_molecule = self.mol_water_2_ang.convertMoleculeToAng()
        self.assertEqual(actual_molecule, expected_molecule)

    def test_write_cube_file_lines(self):
        """Test to see if correct cube file lines are written.
        """
        expected_lines = [
            "    8    0.000000    1.000000    1.000000    1.000000",
            "    1    0.000000    1.000000    0.000000    0.000000",
            "    1    0.000000    0.000000    1.000000    0.000000",
        ]
        actual_lines = self.mol_water_1.writeCubeFileLines()
        self.assertListEqual(actual_lines, expected_lines)

    def test_write_cube_file_lines_error(self):
        """Test to see if expected error is raised if atom id's are not all unique.
        """
        with self.assertRaises(ValueError):
            self.mol_water_2_ang.writeCubeFileLines()

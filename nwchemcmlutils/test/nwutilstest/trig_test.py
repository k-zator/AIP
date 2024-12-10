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
Script containing tests for the Trig module.
"""

import unittest
import logging
import math
import numpy as np
import nwchemcmlutils.nwUtils.Trig as Trig

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.WARN)


class TrigTestCase(unittest.TestCase):
    """Test case for the functions in the Trig script.
    """

    def setUp(self):
        """set up for tests.
        """
        self.atoms = np.array([[1, 2, 3], [4, 5, 6]])
        self.water_atoms = np.array(
            [
                [-10.030764, -1.070314, -0.026721],
                [-8.202256, -1.001565, 0.035731],
                [-10.552420, 0.228376, 1.153020],
            ]
        )
        self.rot_matrix = Trig.generateRotationMatrix(np.array([1, 0, 0]), math.pi / 4)
        self.rot_trans_matrix = Trig.rotateMoleculeRelativeToTranslationVector(
            np.array([[1, 0, 0]]), 180, np.array([0, 1, 0]), np.array([1, 1, 1])
        )

    def tearDown(self):
        """tear down for tests.
        """
        del self.atoms
        del self.water_atoms
        del self.rot_matrix
        del self.rot_trans_matrix

    def test_calc_centroid_vector(self):
        """Test to check calculated translation vector matches the expected
        vector.
        """
        expected_array = np.array([2.5, 3.5, 4.5])
        actual_array = Trig.calculateCentroidVector(self.atoms)
        np.testing.assert_array_equal(actual_array, expected_array)

    def test_gen_rotation_matrix(self):
        """Test to see if expected rotation matrix is generated.
        """
        expected_array = np.array(
            [[1, 0, 0], [0.0, 0.70710678, 0.70710678], [0.0, -0.70710678, 0.70710678]]
        )
        actual_array = self.rot_matrix
        np.testing.assert_array_almost_equal(actual_array, expected_array)

    def test_gen_rot_matrix_error(self):
        """Test to see if expected ValueError is raised when we try to rotate
        about [0, 0, 0].
        """
        with self.assertRaises(ValueError) as err:
            Trig.generateRotationMatrix(np.array([0, 0, 0]), 0)
        expected_args = "Can't rotate about [0, 0, 0]"
        actual_args = err.exception.args[0]
        self.assertEqual(actual_args, expected_args)

    def test_translate_atom(self):
        """test to see if atom coordinates are correctly translated.
        """
        expected_array = np.array([2, 3, 4])
        actual_array = Trig.translateAtom(self.atoms[0], np.array([1, 1, 1]))
        np.testing.assert_array_almost_equal(actual_array, expected_array)

    def test_translate_molecule(self):
        """Test to see if molecule is correctly translated.
        """
        expected_array = np.array([[2, 3, 4], [5, 6, 7]])
        actual_array = Trig.translateMolecule(self.atoms, np.array([1, 1, 1]))
        np.testing.assert_array_equal(actual_array, expected_array)

    def test_rotate_atom_by_matrix(self):
        """Test to see if atom is correctly rotated by rotation matrix.
        """
        expected_array = np.array([1.0, 1.41421356, 0.0])
        actual_array = Trig.rotateAtomByMatrix(np.array([1, 1, 1]), self.rot_matrix)
        np.testing.assert_array_almost_equal(actual_array, expected_array)

    def test_rotate_molecule_by_matrix(self):
        """Test to see if water molecule is corrected rotated by rotation matrix.
        """
        expected_array = np.array(
            [
                [-10.030764, -0.77572089, 0.73793169],
                [-8.202256, -0.68294777, 0.73347904],
                [-10.55242, 0.97679448, 0.65382204],
            ]
        )
        actual_array = Trig.rotateMoleculeByMatrix(self.water_atoms, self.rot_matrix)
        np.testing.assert_array_almost_equal(actual_array, expected_array)

    def test_rotate_molecule_360_x(self):
        """Test to see if molecule is correctly rotated about x axis with
        angle of 360 degrees.
        """
        expected_array = np.array([[1, 2, 3], [4, 5, 6]])
        actual_array = Trig.rotateMolecule(self.atoms, 360, np.array([1, 0, 0]))
        np.testing.assert_array_almost_equal(actual_array, expected_array)

    def test_rotate_molecule_0_x(self):
        """Test to see if molecule is correctly rotated about x axis with
        angle of 0 degrees.
        """
        expected_array = np.array([[1, 2, 3], [4, 5, 6]])
        actual_array = Trig.rotateMolecule(self.atoms, 0, np.array([1, 0, 0]))
        np.testing.assert_array_almost_equal(actual_array, expected_array)

    def test_rotate_molecule_180_y(self):
        """Test to see if molecule is correctly rotated about y axis with
        angle of 180 degrees.
        """
        expected_array = np.array([[-1, 2, -3], [-4, 5, -6]])
        actual_array = Trig.rotateMolecule(self.atoms, 180, np.array([0, 1, 0]))
        np.testing.assert_array_almost_equal(actual_array, expected_array)

    def test_rotate_mol_error(self):
        """Test to see if the expected error is raised if we try to rotate
        a molecule about [0, 0, 0].
        """
        with self.assertRaises(ValueError) as err:
            Trig.rotateMolecule(self.atoms, 0, np.array([0, 0, 0]))
        expected_args = "Can't rotate about [0, 0, 0]"
        actual_args = err.exception.args[0]
        self.assertEqual(actual_args, expected_args)

    def test_rot_and_trans_w_trans_v(self):
        """Test to see if the rotation and translation carried out
        returns the expected final position.
        """
        expected_array = np.array([[1, 0, 2]])
        actual_array = self.rot_trans_matrix
        np.testing.assert_array_almost_equal(actual_array, expected_array)

    def test_rotate_mol_w_trans_vec_error(self):
        """Test to see if the expected error is raised if we try to rotate
        a molecule about [0, 0, 0], whilst also translatingthe molecule.
        """
        with self.assertRaises(ValueError) as err:
            Trig.rotateMoleculeRelativeToTranslationVector(
                self.atoms, 0, np.array([0, 0, 0]), np.array([1, 0, 0])
            )
        expected_args = "Can't rotate about [0, 0, 0]"
        actual_args = err.exception.args[0]
        self.assertEqual(actual_args, expected_args)

    def test_rotate_mol_about_centroid(self):
        """Test to see if molecule is correctly rotated about the centroid by
        45 degrees around the x axis (1,0,0).
        """
        expected_array = np.array(
            [
                [-10.030764, -1.22959716, 0.4168641],
                [-8.202256, -1.13682404, 0.41241145],
                [-10.55242, 0.52291821, 0.33275445],
            ]
        )
        actual_array = Trig.rotateMoleculeAboutCentroid(
            self.water_atoms, 45, np.array([1, 0, 0])
        )
        np.testing.assert_array_almost_equal(actual_array, expected_array)

    def test_rotate_mol_centroid_error(self):
        """Test to see if the expected error is raised when we try to rotate
        molecule about centroid by [0, 0, 0]
        """
        with self.assertRaises(ValueError) as err:
            Trig.rotateMoleculeAboutCentroid(self.atoms, 0, np.array([0, 0, 0]))
        expected_args = "Can't rotate about [0, 0, 0]"
        actual_args = err.exception.args[0]
        self.assertEqual(actual_args, expected_args)

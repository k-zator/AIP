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
Tests for GaussianCube class. This contains the tests for the rotate method.

@author: mdd31
"""

import unittest
import math
import os
import numpy as np
import pathlib
from nwchemcmlutils.DataClasses.DistanceUnits import DistanceUnits
from nwchemcmlutils.DataClasses.Cartesian3D import Cartesian3D
from nwchemcmlutils.DataClasses.Atom import Atom
from nwchemcmlutils.DataClasses.Molecule import Molecule
from nwchemcmlutils.DataClasses.PropertyValue import PropertyValue
from nwchemcmlutils.DataClasses.PropertyValueList import PropertyValueList
from nwchemcmlutils.DataClasses.GaussianCube import GaussianCube


class GaussianCubeRotateAddPlotTestCase(unittest.TestCase):
    """Test case for the functions in the GaussianCube object related to
    file reading and writing.
    """

    def setUp(self):
        """set up for tests.
        """
        self.gaussian_cube = GaussianCube()
        self.parent_directory = pathlib.Path(__file__).parents[1]
        self.expected_filename = (
            (self.parent_directory / "resources/expected_cube-0_000.cube")
            .absolute()
            .as_posix()
        )
        self.expected_rot_cube = (
            (self.parent_directory / "resources/expected_rotated_cube.cube")
            .absolute()
            .as_posix()
        )
        self.gaussian_cube.read(self.expected_filename)
        self.rot_matrix_45_100 = np.array(
            [[1, 0, 0], [0.0, 0.70710678, 0.70710678], [0.0, -0.70710678, 0.70710678]]
        )

    def tearDown(self):
        """tear down for tests.
        """
        del self.gaussian_cube
        del self.rot_matrix_45_100

    def test_rotate_by_matrix(self):
        """test to see if expected cube is returned
        """
        expected_cube = GaussianCube()
        expected_cube.read(self.expected_rot_cube)
        actual_cube = self.gaussian_cube.rotateByMatrix(self.rot_matrix_45_100)
        self.assertEqual(actual_cube.atom_coordinates, expected_cube.atom_coordinates)
        self.assertEqual(actual_cube.property_values, expected_cube.property_values)

    def test_rotate_about_axis(self):
        """Test to see if expected cube is returned.
        """
        expected_cube = GaussianCube()
        expected_cube.read(self.expected_rot_cube)
        actual_cube = self.gaussian_cube.rotateAboutAxis(45, np.array([1, 0, 0]))
        self.assertEqual(actual_cube.atom_coordinates, expected_cube.atom_coordinates)
        self.assertEqual(actual_cube.property_values, expected_cube.property_values)

    def test_add_cube_file(self):
        """Test to see if cube files are correctly added.
        """
        expected_cube = GaussianCube()
        expected_cube.read(
            (self.parent_directory / "resources/full-0_000.cube").absolute().as_posix()
        )
        cube_to_add = GaussianCube()
        cube_to_add.read(
            (self.parent_directory / "resources/part_2-0_000.cube")
            .absolute()
            .as_posix()
        )
        actual_cube = self.gaussian_cube.addGausianCube(cube_to_add)
        self.assertEqual(actual_cube, expected_cube)

    def test_add_cube_file_error(self):
        """TEst to see if expeccted errors are raised if we try to add bad
        cube files.
        """
        with self.assertRaises(ValueError) as val_err:
            bad_cube_1 = GaussianCube(dist_units=DistanceUnits.angstroms)
            self.gaussian_cube.addGausianCube(bad_cube_1)
        actual_err_string_1 = val_err.exception.args[0]
        expected_err_string_1 = "Units do not match"
        self.assertEqual(actual_err_string_1, expected_err_string_1)
        with self.assertRaises(ValueError) as val_err:
            bad_cube_2 = GaussianCube(atom_coordinates=Molecule([]))
            self.gaussian_cube.addGausianCube(bad_cube_2)
        actual_err_string_2 = val_err.exception.args[0]
        expected_err_string_2 = "Molecules do not match"
        self.assertEqual(actual_err_string_2, expected_err_string_2)

    def test_calculate_voxel_volume(self):
        """test to see if voxel volume is calculated as expected.
        """
        expected_value = 0.00446655
        actual_value, actual_units = self.gaussian_cube.calculateVoxelVolume()
        self.assertAlmostEqual(expected_value, actual_value)
        self.assertEqual(DistanceUnits.atomic_units, actual_units)

    def test_min_value(self):
        """Test to see if expected property value is returned.
        """
        expected_prop_value = PropertyValue(
            -0.055775,
            Cartesian3D(
                [-12.845512, -2.061264, -0.508571], unit=DistanceUnits.atomic_units
            ),
        )
        actual_prop_value = self.gaussian_cube.returnMinValue()
        self.assertEqual(actual_prop_value, expected_prop_value)

    def test_max_value(self):
        """Test to see if expected property value is returned.
        """
        expected_prop_value = PropertyValue(
            0.014346,
            Cartesian3D(
                [-12.845512, -0.25694, 0.97535], unit=DistanceUnits.atomic_units
            ),
        )
        actual_prop_value = self.gaussian_cube.returnMaxValue()
        self.assertEqual(actual_prop_value, expected_prop_value)

    def test_return_plot_values(self):
        """Test to see if expected plot values are returned.
        """
        expected_plot_values = np.array(
            [
                [
                    -6.884947e00,
                    -6.884947e00,
                    -6.797553e00,
                    -6.797553e00,
                    -6.797553e00,
                    -6.797553e00,
                    -6.797553e00,
                    -6.797553e00,
                    -6.797553e00,
                    -6.797553e00,
                ],
                [
                    -5.69969936e-01,
                    -5.69969936e-01,
                    -1.09077401e00,
                    -1.00397360e00,
                    -9.17172654e-01,
                    -7.43571295e-01,
                    -3.09568161e-01,
                    -1.35966802e-01,
                    -1.35966802e-01,
                    3.76345568e-02,
                ],
                [
                    -9.46227131e-02,
                    -7.37196826e-03,
                    -2.69124203e-01,
                    3.41631540e-01,
                    4.28882285e-01,
                    5.16133030e-01,
                    -5.30876437e-01,
                    -3.56374948e-01,
                    5.16133030e-01,
                    7.98787766e-02,
                ],
                [
                    -3.41250000e-02,
                    -3.09810000e-02,
                    -5.57750000e-02,
                    -3.90070000e-02,
                    -3.28700000e-02,
                    -2.15760000e-02,
                    -4.22330000e-02,
                    -2.99920000e-02,
                    1.43460000e-02,
                    5.40000000e-05,
                ],
            ]
        )
        actual_plot_values = self.gaussian_cube.returnPlotValues(
            units=DistanceUnits.angstroms
        )
        np.testing.assert_array_almost_equal(actual_plot_values, expected_plot_values)

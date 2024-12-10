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
contains test case for PropertyValue class

@author: mdd31
"""


import unittest
import logging
import numpy as np
from nwchemcmlutils.DataClasses.Vector3D import Vector3D
from nwchemcmlutils.DataClasses.Cartesian3D import Cartesian3D
from nwchemcmlutils.DataClasses.DistanceUnits import DistanceUnits
from nwchemcmlutils.DataClasses.PropertyValue import PropertyValue

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.WARN)


class PropertyValueTestCase(unittest.TestCase):
    """Test Case for the PropertyValue class.
    """

    def setUp(self):
        """set up for tests.
        """
        cartesian_3d_111_au = Cartesian3D([1, 1, 1], unit=DistanceUnits.atomic_units)
        self.property_value_au = PropertyValue(1.0, cartesian_3d_111_au)
        cartesian_3d_ang = Cartesian3D(
            [0.529177249, 0.529177249, 0.529177249], unit=DistanceUnits.angstroms
        )
        self.property_value_ang = PropertyValue(1.0, cartesian_3d_ang)
        cartesian_3d_123_au = Cartesian3D([1, 2, 3], unit=DistanceUnits.atomic_units)
        self.property_value_au_123 = PropertyValue(2.0, cartesian_3d_123_au)
        self.property_value_au_111 = PropertyValue(2.0, cartesian_3d_111_au)
        self.rot_matrix_45_100 = np.array(
            [[1, 0, 0], [0.0, 0.70710678, 0.70710678], [0.0, -0.70710678, 0.70710678]]
        )

    def tearDown(self):
        """tear down for tests.
        """
        del self.property_value_au
        del self.property_value_ang
        del self.property_value_au_123
        del self.property_value_au_111
        del self.rot_matrix_45_100

    def test___repr__(self):
        """Test to see if expected string is returned for the represenation of
        the PropertyValue
        """
        expected_repr = """PropertyValue(1.0, Cartesian3D(Vector3D([[1.],
          [1.],
          [1.]]), unit=DistanceUnits.atomic_units))"""
        actual_repr = repr(self.property_value_au)
        self.assertMultiLineEqual(actual_repr, expected_repr)

    def test___str__(self):
        """Test to see if expected string is returned.
        """
        expected_string = """Value = 1.0
Coordinate = Coordinates are:
x=1.00000 y=1.00000 z=1.00000
in units of: DistanceUnits.atomic_units
"""
        actual_string = str(self.property_value_au)
        self.assertEqual(actual_string, expected_string)

    def test___lt__(self):
        """Test that less than operator returns expected boolean.
        """
        self.assertLess(self.property_value_au, self.property_value_au_123)
        self.assertFalse(self.property_value_au < self.property_value_au)

    def test___lt__errors(self):
        """Tests that expected errors are raised during less than comparison.
        """
        with self.assertRaises(TypeError):
            less_than = self.property_value_ang < self.property_value_au_123
        with self.assertRaises(ValueError):
            less_than = self.property_value_au < self.property_value_au_111

    def test___gt__(self):
        """Test that greater than operator returns expected boolean.
        """
        self.assertGreater(self.property_value_au_123, self.property_value_au)
        self.assertFalse(self.property_value_au > self.property_value_au)

    def test___gt__errors(self):
        """Tests that expected errors are raised during greater than comparison.
        """
        with self.assertRaises(TypeError):
            less_than = self.property_value_ang > self.property_value_au_123
        with self.assertRaises(ValueError):
            less_than = self.property_value_au > self.property_value_au_111

    def test___hash__(self):
        """Tests that expected hash is returned.
        """
        actual_hash = hash(self.property_value_au)
        self.assertTrue(isinstance(actual_hash, int))

    def test__md5_hash(self):
        """Test to see if expected md5 sum is returned.
        """
        expected_md5 = "7dcc654b7869c967650c89d1900b05af"
        actual_md5 = self.property_value_au._md5_hash().hexdigest()
        self.assertEqual(actual_md5, expected_md5)

    def test_translate_property_value(self):
        """Test to see if property value is translated.
        """
        expected_property_value = PropertyValue(
            1.0, Cartesian3D([2, 2, 2], unit=DistanceUnits.atomic_units)
        )
        cartesian_3d = Cartesian3D([1, 1, 1], unit=DistanceUnits.atomic_units)
        actual_property_value = self.property_value_au.translatePropertyValue(
            cartesian_3d
        )
        self.assertEqual(actual_property_value, expected_property_value)
        # now we test to see if we get the expected errors
        with self.assertRaises(ValueError) as err_1:
            self.property_value_au.translatePropertyValue(
                Cartesian3D([1, 1, 1], unit=DistanceUnits.angstroms)
            )
        expected_args_1 = "Incorrect Units"
        actual_args_1 = err_1.exception.args[0]
        self.assertEqual(actual_args_1, expected_args_1)
        with self.assertRaises(AttributeError) as err_2:
            self.property_value_au.translatePropertyValue(Vector3D([1, 1, 1]))
        expected_args_2 = "'Vector3D' object has no attribute 'unit'"
        actual_args_2 = err_2.exception.args[0]
        self.assertEqual(actual_args_2, expected_args_2)

    def test_rotate_property_value(self):
        """Test to see if property value is rotated.
        """
        expected_property_value = PropertyValue(
            1.0, Cartesian3D([1.0, 1.41421356, 0.0], unit=DistanceUnits.atomic_units)
        )
        actual_property_value = self.property_value_au.rotatePropertyValue(
            self.rot_matrix_45_100
        )
        self.assertEqual(actual_property_value, expected_property_value)

    def test_convert_to_au_from_ang(self):
        """Test to see if conversion from an PropertyValue with a Cartesian3D position
        in angstroms is converted to one in atomic units.
        """
        expected_property_value = self.property_value_au
        actual_property_value = self.property_value_ang.convertCoordinateUnitsToAU()
        self.assertEqual(actual_property_value, expected_property_value)

    def test_convert_to_au_from_au(self):
        """Test to see if conversion from an atom with a Cartesian3D position
        in atomic units is converted to one in atomic units.
        """
        expected_property_value = self.property_value_au
        actual_property_value = self.property_value_au.convertCoordinateUnitsToAU()
        self.assertEqual(actual_property_value, expected_property_value)

    def test_convert_to_ang_from_au(self):
        """Test to see if conversion from an atom with a Cartesian3D position
        in atomic units is converted to one in angstroms.
        """
        expected_property_value = self.property_value_ang
        actual_property_value = self.property_value_au.convertCoordinateUnitsToAng()
        self.assertEqual(actual_property_value, expected_property_value)

    def test_convert_to_ang_from_ang(self):
        """Test to see if conversion from an atom with a Cartesian3D position
        in angstroms is converted to one in angstroms.
        """
        expected_property_value = self.property_value_ang
        actual_property_value = self.property_value_ang.convertCoordinateUnitsToAng()
        self.assertEqual(actual_property_value, expected_property_value)

    def test_write_cube_file_line(self):
        """Test to see if expected cube file line is written.
        """
        expected_line = "     1.000000     2.000000     3.000000     2.000000"
        actual_line = self.property_value_au_123.writeCubeFileLine()
        self.assertEqual(actual_line, expected_line)

    def test_return_plot_values(self):
        """Test to see if expected tuple is returned.
        """
        expected_array = np.array([1.0, 2.0, 3.0, 2.0])
        actual_array = self.property_value_au_123.returnValuesForPlotting()
        np.testing.assert_array_almost_equal(actual_array, expected_array)

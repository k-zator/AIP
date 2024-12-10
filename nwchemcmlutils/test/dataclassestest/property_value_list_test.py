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
Script contains test case for the PropertyValueList object.

@author: mdd31
"""

import unittest
import logging
import numpy as np
from nwchemcmlutils.DataClasses.Cartesian3D import Cartesian3D
from nwchemcmlutils.DataClasses.DistanceUnits import DistanceUnits
from nwchemcmlutils.DataClasses.PropertyValue import PropertyValue
from nwchemcmlutils.DataClasses.PropertyValueList import PropertyValueList

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.WARN)


class PropertyValueListTestCase(unittest.TestCase):
    """Test case for the Atom class.
    """

    def setUp(self):
        """Set up for tests.
        """
        cartesian_3d_111_au = Cartesian3D([1, 1, 1], unit=DistanceUnits.atomic_units)
        self.prop_val_1_au = PropertyValue(8.0, cartesian_3d_111_au)
        cartesian_3d_111_ang = Cartesian3D([1, 1, 1], unit=DistanceUnits.angstroms)
        self.prop_val_1_ang = PropertyValue(8.0, cartesian_3d_111_ang)
        cartesian_3d_100_au = Cartesian3D([1, 0, 0], unit=DistanceUnits.atomic_units)
        self.prop_val_2_au = PropertyValue(1.0, cartesian_3d_100_au)
        cartesian_3d_010_au = Cartesian3D([0, 1, 0], unit=DistanceUnits.atomic_units)
        self.prop_val_3_au = PropertyValue(1.0, cartesian_3d_010_au)
        self.prop_val_list_1_au = PropertyValueList(
            [self.prop_val_1_au, self.prop_val_2_au, self.prop_val_3_au]
        )
        prop_val_1_ang = PropertyValue(
            8.0,
            Cartesian3D(
                [0.529177249, 0.529177249, 0.529177249], unit=DistanceUnits.angstroms
            ),
        )
        prop_val_2_ang = PropertyValue(
            1.0, Cartesian3D([0.529177249, 0.0, 0.0], unit=DistanceUnits.angstroms)
        )
        prop_val_3_ang = PropertyValue(
            1.0, Cartesian3D([0.0, 0.529177249, 0.0], unit=DistanceUnits.angstroms)
        )
        self.prop_val_list_1_ang = PropertyValueList(
            [prop_val_1_ang, prop_val_2_ang, prop_val_3_ang],
            unit=DistanceUnits.angstroms,
        )
        self.rot_matrix_45_100 = np.array(
            [[1, 0, 0], [0.0, 0.70710678, 0.70710678], [0.0, -0.70710678, 0.70710678]]
        )

    def tearDown(self):
        """Tear down for tests.
        """
        del self.prop_val_1_ang
        del self.prop_val_1_au
        del self.prop_val_2_au
        del self.prop_val_3_au
        del self.prop_val_list_1_ang
        del self.prop_val_list_1_au
        del self.rot_matrix_45_100

    def test___repr__(self):
        """Test to see if correct representation is returned.
        """
        expected_repr = """PropertyValueList([PropertyValue(8.0, Cartesian3D(Vector3D([[1.],
          [1.],
          [1.]]), unit=DistanceUnits.atomic_units)), PropertyValue(1.0, Cartesian3D(Vector3D([[1.],
          [0.],
          [0.]]), unit=DistanceUnits.atomic_units)), PropertyValue(1.0, Cartesian3D(Vector3D([[0.],
          [1.],
          [0.]]), unit=DistanceUnits.atomic_units))], unit=DistanceUnits.atomic_units)"""
        actual_repr = repr(self.prop_val_list_1_au)
        self.assertMultiLineEqual(actual_repr, expected_repr)

    def test___str__(self):
        """Test to see if expected string is returned.
        """
        expected_string = """PropertyValueList: [PropertyValue(8.0, Cartesian3D(Vector3D([[1.],
          [1.],
          [1.]]), unit=DistanceUnits.atomic_units)), PropertyValue(1.0, Cartesian3D(Vector3D([[1.],
          [0.],
          [0.]]), unit=DistanceUnits.atomic_units)), PropertyValue(1.0, Cartesian3D(Vector3D([[0.],
          [1.],
          [0.]]), unit=DistanceUnits.atomic_units))]"""
        actual_string = str(self.prop_val_list_1_au)
        self.assertMultiLineEqual(actual_string, expected_string)

    def test___eq__(self):
        """Test to see if molecules are equal.
        """
        prop_val_1 = PropertyValue(
            8.0, Cartesian3D([1, 1, 1], unit=DistanceUnits.atomic_units)
        )
        prop_val_2 = PropertyValue(
            1.0, Cartesian3D([1, 0, 0], unit=DistanceUnits.atomic_units)
        )
        prop_val_3 = PropertyValue(
            1.0, Cartesian3D([0, 1, 0], unit=DistanceUnits.atomic_units)
        )
        expected_prop_val_list = PropertyValueList(
            [prop_val_1, prop_val_2, prop_val_3], unit=DistanceUnits.atomic_units
        )
        actual_prop_val_list = self.prop_val_list_1_au
        self.assertEqual(actual_prop_val_list, expected_prop_val_list)

    def test_add_property_value_error(self):
        """Test to see if expected error is raised.
        """
        with self.assertRaises(ValueError):
            self.prop_val_list_1_au.addPropertyValue(self.prop_val_1_ang)

    def test_append_property_value_list(self):
        """Test to see if PropertyValueList is correctly appended.
        """
        prop_val_1 = PropertyValue(
            8.0, Cartesian3D([1, 1, 1], unit=DistanceUnits.atomic_units)
        )
        prop_val_2 = PropertyValue(
            1.0, Cartesian3D([1, 0, 0], unit=DistanceUnits.atomic_units)
        )
        prop_val_3 = PropertyValue(
            1.0, Cartesian3D([0, 1, 0], unit=DistanceUnits.atomic_units)
        )
        expected_prop_val_list = PropertyValueList(
            [prop_val_1, prop_val_2, prop_val_3, prop_val_1, prop_val_2, prop_val_3],
            unit=DistanceUnits.atomic_units,
        )
        actual_prop_val_list = self.prop_val_list_1_au.appendPropertyValueList(
            self.prop_val_list_1_au
        )
        self.assertEqual(actual_prop_val_list, expected_prop_val_list)

    def test_append_prop_val_list_error(self):
        """Test to see if expected error is raised.
        """
        with self.assertRaises(ValueError):
            self.prop_val_list_1_au.appendPropertyValueList(self.prop_val_list_1_ang)

    def test_translate_property_value_list(self):
        """Test to see if PropertyValueList is correctly translated.
        """
        prop_val_1 = PropertyValue(
            8.0, Cartesian3D([1, 1, 2], unit=DistanceUnits.atomic_units)
        )
        prop_val_2 = PropertyValue(
            1.0, Cartesian3D([1, 0, 1], unit=DistanceUnits.atomic_units)
        )
        prop_val_3 = PropertyValue(
            1.0, Cartesian3D([0, 1, 1], unit=DistanceUnits.atomic_units)
        )
        expected_prop_val_list = PropertyValueList(
            [prop_val_1, prop_val_2, prop_val_3], unit=DistanceUnits.atomic_units
        )
        translation_cartesian_3d = Cartesian3D(
            [0, 0, 1], unit=DistanceUnits.atomic_units
        )
        actual_prop_val_list = self.prop_val_list_1_au.translatePropertyValueList(
            translation_cartesian_3d
        )
        self.assertEqual(actual_prop_val_list, expected_prop_val_list)

    def test_trans_prop_val_list_error(self):
        """Test to see if correct errors are raised when tring to translate the
        PropertyValueList
        """
        translation_cartesian_3d = Cartesian3D([0, 0, 1], unit=DistanceUnits.angstroms)
        with self.assertRaises(ValueError):
            self.prop_val_list_1_au.translatePropertyValueList(translation_cartesian_3d)
        with self.assertRaises(AttributeError):
            self.prop_val_list_1_au.translatePropertyValueList([1, 0, 0])

    def test_rotate_prop_val_list(self):
        """Test to see if Property_value_list is correctly rotated.
        """
        prop_val_1 = PropertyValue(
            8.0, Cartesian3D([1.0, 1.41421356, 0.0], unit=DistanceUnits.atomic_units)
        )
        prop_val_2 = PropertyValue(
            1.0, Cartesian3D([1.0, 0.0, 0.0], unit=DistanceUnits.atomic_units)
        )
        prop_val_3 = PropertyValue(
            1.0,
            Cartesian3D(
                [0.0, 0.70710678, -0.70710678], unit=DistanceUnits.atomic_units
            ),
        )
        expected_prop_val_list = PropertyValueList(
            [prop_val_1, prop_val_2, prop_val_3], unit=DistanceUnits.atomic_units
        )
        actual_prop_val_list = self.prop_val_list_1_au.rotatePropertyValueList(
            self.rot_matrix_45_100
        )
        self.assertEqual(actual_prop_val_list, expected_prop_val_list)

    def test_rotate_prop_val_list_about_vec(self):
        """Test to see if PropertyValueList can be rotated about a new origin
        when the translation vector between the current frame of reference and
        the new frame of reference is decided. 
        """
        prop_val_1 = PropertyValue(
            8.0,
            Cartesian3D([1.0, 1.37377345, 0.56903559], unit=DistanceUnits.atomic_units),
        )
        prop_val_2 = PropertyValue(
            1.0,
            Cartesian3D(
                [1.0, -0.04044011, 0.56903559], unit=DistanceUnits.atomic_units
            ),
        )
        prop_val_3 = PropertyValue(
            1.0,
            Cartesian3D(
                [0.0, 0.66666667, -0.13807119], unit=DistanceUnits.atomic_units
            ),
        )
        expected_prop_val_list = PropertyValueList(
            [prop_val_1, prop_val_2, prop_val_3], unit=DistanceUnits.atomic_units
        )
        prop_val_list = self.prop_val_list_1_au
        trans_vector = Cartesian3D(
            [2.0 / 3.0, 2.0 / 3.0, 1 / 3.0], unit=DistanceUnits.atomic_units
        )
        actual_prop_val_list = prop_val_list.rotatePropertyValueListAboutVector(
            self.rot_matrix_45_100, trans_vector
        )
        self.assertEqual(actual_prop_val_list, expected_prop_val_list)

    def test_convert_to_au_from_ang(self):
        """Test conversion to atomic units from angstroms.
        """
        expected_prop_val_list = self.prop_val_list_1_au
        actual_prop_val_list = self.prop_val_list_1_ang.convertPropertyValueListToAU()
        self.assertEqual(actual_prop_val_list, expected_prop_val_list)

    def test_convert_to_au_from_au(self):
        """Test conversion to atomic units from atomic units.
        """
        expected_prop_val_list = self.prop_val_list_1_au
        actual_prop_val_list = self.prop_val_list_1_au.convertPropertyValueListToAU()
        self.assertEqual(actual_prop_val_list, expected_prop_val_list)

    def test_convert_to_ang_from_au(self):
        """Test conversion to angstroms to atomic units.
        """
        expected_prop_val_list = self.prop_val_list_1_ang
        actual_prop_val_list = self.prop_val_list_1_au.convertPropertyValueListToAng()
        self.assertEqual(actual_prop_val_list, expected_prop_val_list)

    def test_convert_to_ang_from_ang(self):
        """Test conversion to angsttroms from angstroms.
        """
        expected_prop_val_list = self.prop_val_list_1_ang
        actual_prop_val_list = self.prop_val_list_1_ang.convertPropertyValueListToAng()
        self.assertEqual(actual_prop_val_list, expected_prop_val_list)

    def test_min_value(self):
        """Test to see if expected property value is returned.
        """
        expected_prop_value = PropertyValue(
            1.0, Cartesian3D([1.0, 0.0, 0.0], unit=DistanceUnits.atomic_units)
        )
        actual_prop_value = self.prop_val_list_1_au.returnMinValue()
        self.assertEqual(actual_prop_value, expected_prop_value)

    def test_max_value(self):
        """Test to see if expected property value is returned.
        """
        expected_prop_value = PropertyValue(
            8.0, Cartesian3D([1.0, 1.0, 1.0], unit=DistanceUnits.atomic_units)
        )
        actual_prop_value = self.prop_val_list_1_au.returnMaxValue()
        self.assertEqual(actual_prop_value, expected_prop_value)

    def test_write_cube_file_lines(self):
        """Test to see if correct lines are written
        """
        expected_lines = [
            "     1.000000     1.000000     1.000000     8.000000",
            "     1.000000     0.000000     0.000000     1.000000",
            "     0.000000     1.000000     0.000000     1.000000",
        ]
        actual_lines = self.prop_val_list_1_au.writecubeFileLines()
        self.assertListEqual(actual_lines, expected_lines)

    def test_return_plot_values(self):
        """Test to see if expected plot values are returned.
        """
        expected_array = np.array(
            [[1.0, 1.0, 0.0], [1.0, 0.0, 1.0], [1.0, 0.0, 0.0], [8.0, 1.0, 1.0]]
        )
        actual_array = self.prop_val_list_1_au.returnPlotValues()
        np.testing.assert_array_almost_equal(actual_array, expected_array)

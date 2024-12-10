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
This contains the tests for the Vector3D object.

@author: mark
"""

import unittest
import logging
import math
import decimal
import numpy as np
from nwchemcmlutils.DataClasses.Vector3D import Vector3D

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.WARN)


class Vector3DTestCase(unittest.TestCase):
    """Test case for the Vector3D class.
    """

    def setUp(self):
        """set up for tests.
        """
        self.vec_3d_100 = Vector3D([1, 0, 0])
        self.vec_3d_010 = Vector3D([0, 1, 0])
        self.vec_3d_111 = Vector3D([1, 1, 1])
        self.vec_3d_123 = Vector3D([1, 2, 3])
        self.vec_3d_000 = Vector3D([0.000000, 0.0, 0.0])
        self.vec_3d_000neg = Vector3D([-0.000000, -0.0, -0.0])
        self.vec_3d_000negnear = Vector3D([-0.000010, -0.0, -0.0])
        self.vec_3d_rounding1 = Vector3D([1.000001, 2.222226, 3.333333])
        self.vec_3d_rounding2 = Vector3D([1.0, 2.22222, 3.33333])
        self.rot_matrix_45_100 = np.array(
            [[1, 0, 0], [0.0, 0.70710678, 0.70710678], [0.0, -0.70710678, 0.70710678]]
        )

    def tearDown(self):
        """tear down for the tests.
        """
        del self.vec_3d_010
        del self.vec_3d_100
        del self.vec_3d_111
        del self.vec_3d_000
        del self.vec_3d_000neg
        del self.vec_3d_000negnear
        del self.vec_3d_rounding1
        del self.vec_3d_rounding2
        del self.rot_matrix_45_100

    def test___init__(self):
        """Test to see if blank Vector3D is initialised.
        """
        expected_vec_3d = self.vec_3d_000
        actual_vec_3d = Vector3D([])
        self.assertEqual(actual_vec_3d, expected_vec_3d)

    def test___mul__(self):
        """test to see if mull has been correctly overloaded and returns
        expected result.
        """
        expected_product = 6.0
        actual_product = self.vec_3d_111 * self.vec_3d_123
        LOGGER.info("actual_product: %s", actual_product)
        self.assertAlmostEqual(actual_product[0][0], expected_product)

    def test___hash__(self):
        """Test to see if expected hash is produced.
        """
        actual_hash = hash(self.vec_3d_123)
        self.assertTrue(isinstance(actual_hash, int))

    def test__bin_ith_value(self):
        """Test to see if correct string for rthe binary represntation is returned.
        """
        expected_result = "00111111100000000000000000000000"
        actual_result = self.vec_3d_111._binaryIthValue(0)
        self.assertEqual(actual_result, expected_result)
        expected_result_0 = self.vec_3d_000._binaryIthValue(0)
        actual_result_neg0 = self.vec_3d_000neg._binaryIthValue(0)
        LOGGER.debug(expected_result_0)
        LOGGER.debug(actual_result_neg0)
        self.assertEqual(actual_result_neg0, expected_result_0)

    def test_binary_values(self):
        """Test to see if the correct tuple of binary values is returned.
        """
        expected_tuple = (
            "00111111100000000000000000000000",
            "01000000000000000000000000000000",
            "01000000010000000000000000000000",
        )
        actual_tuple = self.vec_3d_123.binaryTuple()
        self.assertEqual(actual_tuple, expected_tuple)
        expected_tuple_000 = self.vec_3d_000.binaryTuple()
        actual_tuple_000neg = self.vec_3d_000neg.binaryTuple()
        self.assertEqual(actual_tuple_000neg, expected_tuple_000)
        actual_tuple_near000neground = (
            self.vec_3d_000negnear.roundValuesVector3D().binaryTuple()
        )
        self.assertEqual(actual_tuple_near000neground, expected_tuple_000)
        actual_tuple_near000neg = self.vec_3d_000negnear.binaryTuple()
        self.assertNotEqual(actual_tuple_near000neg, expected_tuple_000)

    def test__md5_hash(self):
        """Test to see if expected md5 sum is produced.
        """
        expected_md5 = "fb738fd57f61e87a19014ceea24b724a"
        actual_md5 = self.vec_3d_123._md5_hash().hexdigest()
        self.assertEqual(actual_md5, expected_md5)
        expected_md5_000 = self.vec_3d_000._md5_hash().hexdigest()
        actual_md5_000neg = self.vec_3d_000neg._md5_hash().hexdigest()
        self.assertEqual(actual_md5_000neg, expected_md5_000)
        actual_md5_000nearneg = self.vec_3d_000negnear._md5_hash().hexdigest()
        self.assertEqual(actual_md5_000nearneg, expected_md5_000)

    def test__round_ith_value(self):
        """Test to see if expected value is returned.
        """
        expected_value = 1.00000
        actual_value = self.vec_3d_rounding1._roundIthValue(
            0, decimal_val=decimal.Decimal("0.00001"), rounding=decimal.ROUND_DOWN
        )
        self.assertEqual(actual_value, expected_value)

    def test_round_values_tuple(self):
        """Test to see that expected tuple is returned of rounded numbers.
        """
        expected_value = (1.00000, 2.22222, 3.33333)
        actual_value = self.vec_3d_rounding1.roundValuesTuple()
        self.assertEqual(actual_value, expected_value)

    def test_round_values_vector_3d(self):
        """Test to see that expected Vector3D is returned for the rounded
        numbers.
        """
        expected_vector_3d = Vector3D([1.00, 2.22, 3.33])
        actual_vector_3d = self.vec_3d_rounding1.roundValuesVector3D()
        self.assertEqual(actual_vector_3d, expected_vector_3d)

    def test_length(self):
        """Test to see expected length is returned.
        """
        expected_length = math.sqrt(3)
        actual_length = self.vec_3d_111.length()
        self.assertAlmostEqual(actual_length, expected_length)

    def test_normalise(self):
        """Test to see if expected normalised Vector3D is returned.
        """
        expected_vector_3d = Vector3D(
            [1.0 / math.sqrt(14.0), 2.0 / math.sqrt(14.0), 3.0 / math.sqrt(14.0)]
        )
        actual_vector_3d = self.vec_3d_123.normalise()
        self.assertAlmostEqual(actual_vector_3d, expected_vector_3d)

    def test_cross(self):
        """Test to see if expected cross product is returned.
        """
        expected_vector_3d = Vector3D([0, 0, 1])
        actual_vector_3d = self.vec_3d_100.cross(self.vec_3d_010)
        self.assertAlmostEqual(actual_vector_3d, expected_vector_3d)

    def test_rotate(self):
        """Test to see if expected Vector3D is returned after rotation.
        """
        expected_vector_3d = Vector3D([1.0, 1.41421356, 0.0])
        actual_vector_3d = self.vec_3d_111.rotate(self.rot_matrix_45_100)
        self.assertAlmostEqual(actual_vector_3d, expected_vector_3d)

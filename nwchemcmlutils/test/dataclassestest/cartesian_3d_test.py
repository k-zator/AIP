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
Script containing tests for Cartesian3D object.

@author: mark
"""

import unittest
import logging
import math
import numpy as np
import nwchemcmlutils.DataClasses.Cartesian3D as Cartesian3D
from nwchemcmlutils.DataClasses.DistanceUnits import DistanceUnits


logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.WARN)


class Cartesian3DTestCase(unittest.TestCase):
    """Test case for the Cartesian3D class.
    """

    def setUp(self):
        """set up for tests.
        """
        self.coord_111_au = Cartesian3D.Cartesian3D([1, 1, 1])
        self.coord_111_ang = Cartesian3D.Cartesian3D(
            [1, 1, 1], unit=DistanceUnits.angstroms
        )
        self.coord_123_au = Cartesian3D.Cartesian3D([1, 2, 3])
        self.rot_matrix_45_100 = np.array(
            [[1, 0, 0], [0.0, 0.70710678, 0.70710678], [0.0, -0.70710678, 0.70710678]]
        )

    def tearDown(self):
        """tear down for tests.
        """
        del self.coord_111_au
        del self.coord_111_ang
        del self.coord_123_au
        del self.rot_matrix_45_100

    def test___init__(self):
        """Test to see if objects have been correctly initialised.
        """
        self.assertIsInstance(self.coord_111_au, Cartesian3D.Cartesian3D)
        self.assertIsInstance(self.coord_123_au, Cartesian3D.Cartesian3D)

    def test___repr__(self):
        """Test to see correct representation is returned.
        """
        expected_repr = """Cartesian3D(Vector3D([[1.],
          [1.],
          [1.]]), unit=DistanceUnits.atomic_units)"""
        actual_repr = repr(self.coord_111_au)
        self.assertMultiLineEqual(actual_repr, expected_repr)

    def test___str__(self):
        """Test to check the correct str is returned.
        """
        expected_string = """Coordinates are:
x=1.00000 y=2.00000 z=3.00000
in units of: DistanceUnits.atomic_units
"""
        actual_string = str(self.coord_123_au)
        self.assertEqual(actual_string, expected_string)

    def test___hash__(self):
        """Test to see if expected hash is produced.
        """
        actual_hash = hash(self.coord_111_au)
        self.assertTrue(isinstance(actual_hash, int))

    def test___hash__equivalence(self):
        """Test to see if hashes of 2 vectors that are equal by '==' also have
        the same hash values.
        """
        other_cartesian_3d = Cartesian3D.Cartesian3D([1.000001, 1.0, 1.0])
        self.assertEqual(self.coord_111_au, other_cartesian_3d)
        hash_1 = hash(self.coord_111_au)
        hash_2 = hash(other_cartesian_3d)
        self.assertEqual(hash_1, hash_2)

    def test__md5_hash(self):
        """Test to see if expected md5 sum is produced.
        """
        expected_md5 = "89ec98926d7dcb3dd7c3f6b695264c6c"
        actual_md5 = self.coord_123_au._md5_hash().hexdigest()
        self.assertEqual(actual_md5, expected_md5)

    def test_length(self):
        """Test to see if correct length is returned, and unit.
        """
        expected_unit = DistanceUnits.atomic_units
        expected_length = math.sqrt(3)
        actual_length, actual_unit = self.coord_111_au.length()
        self.assertEqual(actual_unit, expected_unit)
        self.assertEqual(actual_length, expected_length)

    def test_add(self):
        """test to check that the addition of 2 Cartesian3D returns the expected
        value.
        """
        expected_value = Cartesian3D.Cartesian3D([2, 3, 4])
        actual_value = self.coord_111_au.add(self.coord_123_au)
        self.assertEqual(actual_value, expected_value)

    def test_add_value_error(self):
        """Test to see that add returns the expected error when 2 Cartesian3D
        objects are added with the wrong units.
        """
        with self.assertRaises(ValueError):
            self.coord_111_ang.add(self.coord_111_au)

    def test_convert_to_au_from_ang(self):
        """Test conversion to atomic units from angstroms.
        """
        expected_array = Cartesian3D.Cartesian3D(
            [1.88972599, 1.88972599, 1.88972599], unit=DistanceUnits.atomic_units
        )
        actual_cartesian_3d = self.coord_111_ang.convertToAU()
        self.assertAlmostEqual(actual_cartesian_3d, expected_array)

    def test_convert_to_au_from_au(self):
        """Test conversion to atomic units from atomic units.
        """
        expected_cartesian_3d = Cartesian3D.Cartesian3D(
            [1.0, 1.0, 1.0], unit=DistanceUnits.atomic_units
        )
        actual_cartesian_3d = self.coord_111_au.convertToAU()
        self.assertEqual(actual_cartesian_3d, expected_cartesian_3d)

    def test_convert_to_ang_from_au(self):
        """Test conversion to angstroms from atomic units.
        """
        expected_cartesian_3d = Cartesian3D.Cartesian3D(
            [0.529177249, 0.529177249, 0.529177249], unit=DistanceUnits.angstroms
        )
        actual_cartesian_3d = self.coord_111_au.convertToAng()
        self.assertEqual(actual_cartesian_3d, expected_cartesian_3d)

    def test_convert_to_ang_from_ang(self):
        """Test conversion to angstroms from anstroms.
        """
        expected_cartesian_3d = Cartesian3D.Cartesian3D(
            [1.0, 1.0, 1.0], unit=DistanceUnits.angstroms
        )
        actual_cartesian_3d = self.coord_111_ang.convertToAng()
        self.assertEqual(actual_cartesian_3d, expected_cartesian_3d)

    def test_rotate(self):
        """Test rotation returns expected Cartesian3D.
        """
        expected_cartesian_3d = Cartesian3D.Cartesian3D(
            [1.0, 1.41421356, 0.0], unit=DistanceUnits.atomic_units
        )
        actual_cartesian_3d = self.coord_111_au.rotate(self.rot_matrix_45_100)
        self.assertEqual(actual_cartesian_3d, expected_cartesian_3d)

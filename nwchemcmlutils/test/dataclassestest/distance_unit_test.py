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
Script contains test case for the DistanceUnits class.

@author: mark
"""

import unittest
import logging
import nwchemcmlutils.DataClasses.DistanceUnits as DistanceUnits

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.WARN)


class DistanceUnitTestCase(unittest.TestCase):
    """Test case for the DistanceUnit class.
    """

    def setUp(self):
        """set up for tests.
        """
        self.distance_unit_ang = DistanceUnits.DistanceUnits(2)
        self.distance_unit_au = DistanceUnits.DistanceUnits(1)

    def tearDown(self):
        """tear down for tests.
        """
        del self.distance_unit_ang
        del self.distance_unit_au

    def test_init(self):
        """Test initialisation.
        """
        self.assertIsInstance(self.distance_unit_ang, DistanceUnits.DistanceUnits)
        self.assertIsInstance(self.distance_unit_au, DistanceUnits.DistanceUnits)

    def test___repr__(self):
        """Test to see if expected repr is propduced.
        """
        expected_repr = "DistanceUnits.atomic_units"
        actual_repr = repr(self.distance_unit_au)
        self.assertEqual(actual_repr, expected_repr)

    def test_describe(self):
        """Test to see expected output is returned.
        """
        expected_output = ("atomic_units", 1)
        actual_output = self.distance_unit_au.describe()
        self.assertEqual(actual_output, expected_output)

    def test_con_factor_to_ang_from_au(self):
        """Test to see if correct conversion factor is returned.
        """
        expected_conversion_factor = 0.529177249
        actual_conversion_factor = self.distance_unit_au.conversionFactorToAng()
        self.assertAlmostEqual(actual_conversion_factor, expected_conversion_factor)

    def test_con_factor_to_ang_from_ang(self):
        """Test to see if correct conversion factor is returned.
        """
        expected_conversion_factor = 1.0
        actual_conversion_factor = self.distance_unit_ang.conversionFactorToAng()
        self.assertAlmostEqual(actual_conversion_factor, expected_conversion_factor)

    def test_con_factor_to_au_from_au(self):
        """Test to see if correct conversion factor is returned.
        """
        expected_conversion_factor = 1.0
        actual_conversion_factor = self.distance_unit_au.conversionFactorToAU()
        self.assertAlmostEqual(actual_conversion_factor, expected_conversion_factor)

    def test_con_factor_to_au_from_ang(self):
        """Test to see if correct conversion factor is returned.
        """
        expected_conversion_factor = 1.0 / 0.529177249
        actual_conversion_factor = self.distance_unit_ang.conversionFactorToAU()
        self.assertAlmostEqual(actual_conversion_factor, expected_conversion_factor)

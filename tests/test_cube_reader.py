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
import filereader.cube_reader as reader
import numpy as np

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.WARN)

SCALE = 0.529177

class TestCubeReader(unittest.TestCase):
    """Test case for structure to cml script methods.
    """

    def setUp(self):
        """Set up before tests.
        """
        parent_directory = pathlib.Path(__file__).resolve().parents[0]
        cube_filename = (
            (parent_directory / "test_files/example_merged.cube").absolute().as_posix())
        self.cube = reader.CubeReader(cube_filename)

    def tearDown(self):
        """Clean up after tests.
        """
        del self.cube

    def test_line_index(self):
        coordinates = self.cube.xyz[0]
        expected_coordinates = np.array(
                [
                    -5.083109*SCALE,
                    0.941493*SCALE,
                    2.340498*SCALE
                    ])
        self.assertAlmostEqual(coordinates[0], expected_coordinates[0])
        self.assertAlmostEqual(coordinates[1], expected_coordinates[1])
        self.assertAlmostEqual(coordinates[2], expected_coordinates[2])

    def test_maximum_MEP_value(self):
        maximum = self.cube.values.max()
        self.assertEqual(0.403625, maximum)

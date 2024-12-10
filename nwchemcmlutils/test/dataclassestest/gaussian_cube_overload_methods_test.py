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
Script contains the tests for the GaussianCube object's protected methods.
@author: mark
"""

import unittest
import math
import os
import numpy as np
import pathlib

from nwchemcmlutils.DataClasses.GaussianCube import GaussianCube


class GaussianCubeOverloadMethodsTestCase(unittest.TestCase):
    """Test case for the functions in the GaussianCube object related to
    file reading and writing.
    """

    def setUp(self):
        """set up for tests.
        """
        self.parent_directory = pathlib.Path(__file__).parents[1]
        self.expected_filename = (
            (self.parent_directory / "resources/expected_cube-0_000.cube")
            .absolute()
            .as_posix()
        )
        with open(self.expected_filename, "r") as cube_file_in:
            cube_file_contents = cube_file_in.readlines()
        self.cube_file_contents = cube_file_contents
        self.gaussian_cube = GaussianCube()

    def tearDown(self):
        """tear down for tests.
        """
        del self.cube_file_contents
        del self.gaussian_cube

    def test___init__(self):
        """Test to see if correct GaussianCube file is generated on
        initialisation.
        """
        self.assertIsInstance(self.gaussian_cube, GaussianCube)

    def test___repr__(self):
        """Test to see if expected repr is produced.
        """
        expected_string = """GaussianCube(title_line_one='Title',
                        title_line_two='Description of property stored in cubefile',
                        dist_units=DistanceUnits.atomic_units,
                        number_of_atoms=None,
                        grid_origin=None,
                        number_of_x_grid_points=None,
                        number_of_y_grid_points=None,
                        number_of_z_grid_points=None,
                        x_grid_vector=None,
                        y_grid_vector=None,
                        z_grid_vector=None,
                        atom_coordinates=None,
                        property_values=None)
"""
        actual_string = repr(self.gaussian_cube)
        self.assertMultiLineEqual(actual_string, expected_string)

    def test___str__(self):
        """Test to see if the expected string output is produced for the
        initialised object.
        """
        expected_string = (
            "numberOfAtoms: None\nDistance Units: DistanceU"
            + "nits.atomic_units\ngridOrigin: None\nnumberOf"
            + "XGridPoints: None\nxGridVector: None\nnumberOfY"
            + "GridPoints: None\nyGridVector: None\nnumberOfZG"
            + "ridPoints: None\nzGridVector: None\natomCoordin"
            + "ates: None\ngrid: None"
        )
        actual_string = str(self.gaussian_cube)
        self.assertEqual(actual_string, expected_string)

    def test___eq__(self):
        """Test to see if GaussianCube 
        """
        expected_cube = GaussianCube()
        expected_cube.read(self.expected_filename)
        actual_cube = self.gaussian_cube
        actual_cube.read(self.expected_filename)
        self.assertEqual(actual_cube, expected_cube)

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
Tests for GaussianCube class. This contains the tests for the reading and
writing methods.

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


class GaussianCubeReadWriteTestCase(unittest.TestCase):
    """Test case for the functions in the GaussianCube object related to
    file reading and writing.
    """

    def setUp(self):
        """set up for tests.
        """
        self.parent_directory = pathlib.Path(__file__).parents[1]
        with open(
            (self.parent_directory / "resources/expected_cube-0_000.cube")
            .absolute()
            .as_posix(),
            "r",
        ) as cube_file_in:
            cube_file_contents = cube_file_in.readlines()
        self.cube_file_contents = cube_file_contents
        self.gaussian_cube = GaussianCube()
        self.expected_filename = (
            (self.parent_directory / "resources/expected_cube-0_000.cube")
            .absolute()
            .as_posix()
        )

    def tearDown(self):
        """tear down for tests.
        """
        del self.cube_file_contents
        del self.gaussian_cube
        if os.path.isfile("actual_cube-0_000.cube"):
            os.remove("actual_cube-0_000.cube")

    def test__read_title_lines(self):
        """test to see if correct title lines are read in.
        """
        expected_title_line_one = "Marat"
        expected_title_line_two = "Gaussian Cube file"
        title_lines = self.cube_file_contents[0:2]
        self.gaussian_cube._read_title_lines(title_lines)
        actual_title_line_one = self.gaussian_cube.title_line_one
        actual_title_line_two = self.gaussian_cube.title_line_two
        self.assertEqual(actual_title_line_one, expected_title_line_one)
        self.assertEqual(actual_title_line_two, expected_title_line_two)
        with self.assertRaises(ValueError) as err:
            self.gaussian_cube._read_title_lines(title_lines[0])
        expected_args = "Incorrect number of title lines."
        actual_args = err.exception.args[0]
        self.assertEqual(actual_args, expected_args)

    def test__read_atoms_and_origin(self):
        """Test to see if correct number of atoms and grid origin are read in.
        """
        expected_number_of_atoms = 3
        expected_grid_origin = Cartesian3D(
            [-14.331872, -4.849766, -3.806173], unit=DistanceUnits.atomic_units
        )
        atom_origin_line = self.cube_file_contents[2]
        self.gaussian_cube._read_atoms_and_origin(atom_origin_line)
        actual_number_of_atoms = self.gaussian_cube.number_of_atoms
        actual_grid_origin = self.gaussian_cube.grid_origin
        self.assertEqual(actual_number_of_atoms, expected_number_of_atoms)
        self.assertEqual(actual_grid_origin, expected_grid_origin)
        with self.assertRaises(ValueError) as err:
            self.gaussian_cube._read_atoms_and_origin(self.cube_file_contents[0])
        expected_args = "Incorrect number of values in line."
        actual_args = err.exception.args[0]
        self.assertEqual(actual_args, expected_args)

    def test__read_x_grid_points_vector(self):
        """Test to see if correct number of x grid points and x grid vector are
        returned.
        """
        expected_number_of_points = 61
        expected_x_grid_vector = Cartesian3D(
            [0.165151, 0.000000, 0.000000], unit=DistanceUnits.atomic_units
        )
        x_grid_line = self.cube_file_contents[3]
        self.gaussian_cube._read_x_grid_points_vector(x_grid_line)
        actual_number_of_points = self.gaussian_cube.number_of_x_grid_points
        actual_x_grid_vector = self.gaussian_cube.x_grid_vector
        self.assertEqual(actual_number_of_points, expected_number_of_points)
        self.assertEqual(actual_x_grid_vector, expected_x_grid_vector)
        with self.assertRaises(ValueError) as err:
            self.gaussian_cube._read_x_grid_points_vector(self.cube_file_contents[0])
        expected_args = "Incorrect number of values in line."
        actual_args = err.exception.args[0]
        self.assertEqual(actual_args, expected_args)

    def test__read_y_grid_points_vector(self):
        """Test to see if correct number of x grid points and x grid vector are
        returned.
        """
        expected_number_of_points = 55
        expected_y_grid_vector = Cartesian3D(
            [0.000000, 0.164030, 0.000000], unit=DistanceUnits.atomic_units
        )
        y_grid_line = self.cube_file_contents[4]
        self.gaussian_cube._read_y_grid_points_vector(y_grid_line)
        actual_number_of_points = self.gaussian_cube.number_of_y_grid_points
        actual_y_grid_vector = self.gaussian_cube.y_grid_vector
        self.assertEqual(actual_number_of_points, expected_number_of_points)
        self.assertEqual(actual_y_grid_vector, expected_y_grid_vector)
        with self.assertRaises(ValueError) as err:
            self.gaussian_cube._read_y_grid_points_vector(self.cube_file_contents[0])
        expected_args = "Incorrect number of values in line."
        actual_args = err.exception.args[0]
        self.assertEqual(actual_args, expected_args)

    def test__read_z_grid_points_vector(self):
        """Test to see if correct number of x grid points and x grid vector are
        returned.
        """
        expected_number_of_points = 54
        expected_z_grid_vector = Cartesian3D(
            [0.000000, 0.000000, 0.164880], unit=DistanceUnits.atomic_units
        )
        z_grid_line = self.cube_file_contents[5]
        self.gaussian_cube._read_z_grid_points_vector(z_grid_line)
        actual_number_of_points = self.gaussian_cube.number_of_z_grid_points
        actual_z_grid_vector = self.gaussian_cube.z_grid_vector
        self.assertEqual(actual_number_of_points, expected_number_of_points)
        self.assertEqual(actual_z_grid_vector, expected_z_grid_vector)
        with self.assertRaises(ValueError) as err:
            self.gaussian_cube._read_z_grid_points_vector(self.cube_file_contents[0])
        expected_args = "Incorrect number of values in line."
        actual_args = err.exception.args[0]
        self.assertEqual(actual_args, expected_args)

    def test__read_atom_line(self):
        """Test to see if correct array is returned for given line.
        """
        expected_atom = Atom(
            8,
            Cartesian3D(
                [-10.030764, -1.070314, -0.026721], unit=DistanceUnits.atomic_units
            ),
            id_num=0,
        )
        atom_line = self.cube_file_contents[6]
        actual_atomy = self.gaussian_cube._read_atom_line(atom_line)
        self.assertEqual(actual_atomy, expected_atom)
        with self.assertRaises(ValueError) as err:
            self.gaussian_cube._read_atom_line(self.cube_file_contents[0])
        expected_args = "Incorrect number of values."
        actual_args = err.exception.args[0]
        self.assertEqual(actual_args, expected_args)

    def test__read_atom_lines(self):
        """Test to see if correct array is produced.
        """
        atom_1 = Atom(
            8,
            Cartesian3D(
                [-10.030764, -1.070314, -0.026721], unit=DistanceUnits.atomic_units
            ),
            id_num=0,
        )
        atom_2 = Atom(
            1,
            Cartesian3D(
                [-8.202256, -1.001565, 0.035731], unit=DistanceUnits.atomic_units
            ),
            id_num=1,
        )
        atom_3 = Atom(
            1,
            Cartesian3D(
                [-10.552420, 0.228376, 1.153020], unit=DistanceUnits.atomic_units
            ),
            id_num=2,
        )
        expected_molecule = Molecule(
            [atom_1, atom_2, atom_3], unit=DistanceUnits.atomic_units
        )
        atom_lines = self.cube_file_contents[6:9]
        self.gaussian_cube.number_of_atoms = 3
        self.gaussian_cube._read_atom_lines(atom_lines)
        actual_molecule = self.gaussian_cube.atom_coordinates
        self.assertEqual(actual_molecule, expected_molecule)
        with self.assertRaises(ValueError) as err:
            self.gaussian_cube._read_atom_lines(self.cube_file_contents[0])
        expected_args = "Incorrect number of atom_lines."
        actual_args = err.exception.args[0]
        self.assertEqual(actual_args, expected_args)
        self.gaussian_cube.number_of_atoms = None
        with self.assertRaises(AttributeError) as err:
            self.gaussian_cube._read_atom_lines(atom_lines)
        expected_args = "No number of atoms"
        actual_args = err.exception.args[0]
        self.assertEqual(actual_args, expected_args)

    def test__read_property_value_line(self):
        """Test to see if the value_line is read in correctly.
        """
        expected_prop_val = PropertyValue(
            -0.034125,
            Cartesian3D(
                [-13.010663, -1.077087, -0.178811], unit=DistanceUnits.atomic_units
            ),
        )
        value_line = self.cube_file_contents[9]
        actual_prop_val = self.gaussian_cube._read_property_value_line(value_line)
        self.assertEqual(actual_prop_val, expected_prop_val)
        with self.assertRaises(ValueError) as err:
            self.gaussian_cube._read_property_value_line(self.cube_file_contents[0])
        expected_args = "Incorrect number of values."
        actual_args = err.exception.args[0]
        self.assertEqual(actual_args, expected_args)

    def test__read_property_values_lines(self):
        """Test to see if the property values are read in as expected.
        """
        prop_val_1 = PropertyValue(
            -0.034125,
            Cartesian3D(
                [-13.010663, -1.077087, -0.178811], unit=DistanceUnits.atomic_units
            ),
        )
        prop_val_2 = PropertyValue(
            -0.030981,
            Cartesian3D(
                [-13.010663, -1.077087, -0.013931], unit=DistanceUnits.atomic_units
            ),
        )
        prop_val_3 = PropertyValue(
            -0.055775,
            Cartesian3D(
                [-12.845512, -2.061264, -0.508571], unit=DistanceUnits.atomic_units
            ),
        )
        prop_val_4 = PropertyValue(
            -0.039007,
            Cartesian3D(
                [-12.845512, -1.897235, 0.645590], unit=DistanceUnits.atomic_units
            ),
        )
        prop_val_5 = PropertyValue(
            -0.032870,
            Cartesian3D(
                [-12.845512, -1.733205, 0.810470], unit=DistanceUnits.atomic_units
            ),
        )
        prop_val_6 = PropertyValue(
            -0.021576,
            Cartesian3D(
                [-12.845512, -1.405146, 0.975350], unit=DistanceUnits.atomic_units
            ),
        )
        prop_val_7 = PropertyValue(
            -0.042233,
            Cartesian3D(
                [-12.845512, -0.584999, -1.003211], unit=DistanceUnits.atomic_units
            ),
        )
        prop_val_8 = PropertyValue(
            -0.029992,
            Cartesian3D(
                [-12.845512, -0.256940, -0.673451], unit=DistanceUnits.atomic_units
            ),
        )
        prop_val_9 = PropertyValue(
            0.014346,
            Cartesian3D(
                [-12.845512, -0.256940, 0.975350], unit=DistanceUnits.atomic_units
            ),
        )
        prop_val_10 = PropertyValue(
            0.000054,
            Cartesian3D(
                [-12.845512, 0.071119, 0.150949], unit=DistanceUnits.atomic_units
            ),
        )
        expected_prop_list = PropertyValueList(
            [
                prop_val_1,
                prop_val_2,
                prop_val_3,
                prop_val_4,
                prop_val_5,
                prop_val_6,
                prop_val_7,
                prop_val_8,
                prop_val_9,
                prop_val_10,
            ],
            unit=DistanceUnits.atomic_units,
        )
        property_lines = self.cube_file_contents[9:]
        self.gaussian_cube._read_property_values_lines(property_lines)
        actual_prop_list = self.gaussian_cube.property_values
        self.assertEqual(actual_prop_list, expected_prop_list)

    def test_read(self):
        """Test to see if the read function populates the object as expected.
        """
        expected_title_line_one = "Marat"
        expected_title_line_two = "Gaussian Cube file"
        expected_number_of_atoms = 3
        expected_grid_origin = Cartesian3D(
            [-14.331872, -4.849766, -3.806173], unit=DistanceUnits.atomic_units
        )
        expected_number_of_x_points = 61
        expected_x_grid_vector = Cartesian3D(
            [0.165151, 0.000000, 0.000000], unit=DistanceUnits.atomic_units
        )
        expected_number_of_y_points = 55
        expected_y_grid_vector = Cartesian3D(
            [0.000000, 0.164030, 0.000000], unit=DistanceUnits.atomic_units
        )
        expected_number_of_z_points = 54
        expected_z_grid_vector = Cartesian3D(
            [0.000000, 0.000000, 0.164880], unit=DistanceUnits.atomic_units
        )
        atom_1 = Atom(
            8,
            Cartesian3D(
                [-10.030764, -1.070314, -0.026721], unit=DistanceUnits.atomic_units
            ),
            id_num=0,
        )
        atom_2 = Atom(
            1,
            Cartesian3D(
                [-8.202256, -1.001565, 0.035731], unit=DistanceUnits.atomic_units
            ),
            id_num=1,
        )
        atom_3 = Atom(
            1,
            Cartesian3D(
                [-10.552420, 0.228376, 1.153020], unit=DistanceUnits.atomic_units
            ),
            id_num=2,
        )
        expected_molecule = Molecule(
            [atom_1, atom_2, atom_3], unit=DistanceUnits.atomic_units
        )
        prop_val_1 = PropertyValue(
            -0.034125,
            Cartesian3D(
                [-13.010663, -1.077087, -0.178811], unit=DistanceUnits.atomic_units
            ),
        )
        prop_val_2 = PropertyValue(
            -0.030981,
            Cartesian3D(
                [-13.010663, -1.077087, -0.013931], unit=DistanceUnits.atomic_units
            ),
        )
        prop_val_3 = PropertyValue(
            -0.055775,
            Cartesian3D(
                [-12.845512, -2.061264, -0.508571], unit=DistanceUnits.atomic_units
            ),
        )
        prop_val_4 = PropertyValue(
            -0.039007,
            Cartesian3D(
                [-12.845512, -1.897235, 0.645590], unit=DistanceUnits.atomic_units
            ),
        )
        prop_val_5 = PropertyValue(
            -0.032870,
            Cartesian3D(
                [-12.845512, -1.733205, 0.810470], unit=DistanceUnits.atomic_units
            ),
        )
        prop_val_6 = PropertyValue(
            -0.021576,
            Cartesian3D(
                [-12.845512, -1.405146, 0.975350], unit=DistanceUnits.atomic_units
            ),
        )
        prop_val_7 = PropertyValue(
            -0.042233,
            Cartesian3D(
                [-12.845512, -0.584999, -1.003211], unit=DistanceUnits.atomic_units
            ),
        )
        prop_val_8 = PropertyValue(
            -0.029992,
            Cartesian3D(
                [-12.845512, -0.256940, -0.673451], unit=DistanceUnits.atomic_units
            ),
        )
        prop_val_9 = PropertyValue(
            0.014346,
            Cartesian3D(
                [-12.845512, -0.256940, 0.975350], unit=DistanceUnits.atomic_units
            ),
        )
        prop_val_10 = PropertyValue(
            0.000054,
            Cartesian3D(
                [-12.845512, 0.071119, 0.150949], unit=DistanceUnits.atomic_units
            ),
        )
        expected_prop_list = PropertyValueList(
            [
                prop_val_1,
                prop_val_2,
                prop_val_3,
                prop_val_4,
                prop_val_5,
                prop_val_6,
                prop_val_7,
                prop_val_8,
                prop_val_9,
                prop_val_10,
            ],
            unit=DistanceUnits.atomic_units,
        )
        expected_cube = GaussianCube(
            title_line_one=expected_title_line_one,
            title_line_two=expected_title_line_two,
            number_of_atoms=expected_number_of_atoms,
            grid_origin=expected_grid_origin,
            number_of_x_grid_points=expected_number_of_x_points,
            number_of_y_grid_points=expected_number_of_y_points,
            number_of_z_grid_points=expected_number_of_z_points,
            x_grid_vector=expected_x_grid_vector,
            y_grid_vector=expected_y_grid_vector,
            z_grid_vector=expected_z_grid_vector,
            atom_coordinates=expected_molecule,
            property_values=expected_prop_list,
        )

        self.gaussian_cube.read(self.expected_filename)
        actual_cube = self.gaussian_cube
        self.assertEqual(actual_cube, expected_cube)

    def test_write_title_lines_list(self):
        """Test to see if expected title lines are added to list.
        """
        expected_title_list = [" Marat", " Gaussian Cube file"]
        self.gaussian_cube.read(self.expected_filename)
        actual_title_list = self.gaussian_cube._write_title_lines_to_list()
        self.assertListEqual(actual_title_list, expected_title_list)

    def test_write_atom_origin_line(self):
        """Test to see if expected line is produced for the number of atoms and
        grid origin.
        """
        expected_line = "    3  -14.331872   -4.849766   -3.806173"
        self.gaussian_cube.read(self.expected_filename)
        actual_line = self.gaussian_cube._write_atoms_and_origin_line()
        self.assertEqual(actual_line, expected_line)

    def test_write_x_grid_points_vector(self):
        """Test to see if expected line is produced for the number of x grid
        points and x grid vector.
        """
        expected_x_line = "   61    0.165151    0.000000    0.000000"
        self.gaussian_cube.read(self.expected_filename)
        actual_x_line = self.gaussian_cube._write_x_grid_points_vector()
        self.assertEqual(actual_x_line, expected_x_line)

    def test_write_y_grid_points_vector(self):
        """Test to see if expected line is produced for the number of y grid
        points and y grid vector.
        """
        expected_y_line = "   55    0.000000    0.164030    0.000000"
        self.gaussian_cube.read(self.expected_filename)
        actual_y_line = self.gaussian_cube._write_y_grid_points_vector()
        self.assertEqual(actual_y_line, expected_y_line)

    def test_write_z_grid_points_vector(self):
        """Test to see if expected line is produced for the number of z grid
        points and z grid vector.
        """
        expected_z_line = "   54    0.000000    0.000000    0.164880"
        self.gaussian_cube.read(self.expected_filename)
        actual_z_line = self.gaussian_cube._write_z_grid_points_vector()
        self.assertEqual(actual_z_line, expected_z_line)

    def test_write_atom_line_list(self):
        """Test to see if expected atom lines are written.
        """
        expected_atom_line_list = [
            "    8    0.000000  -10.030764   -1.070314   -0.026721",
            "    1    0.000000   -8.202256   -1.001565    0.035731",
            "    1    0.000000  -10.552420    0.228376    1.153020",
        ]
        self.gaussian_cube.read(self.expected_filename)
        actual_atom_line_list = self.gaussian_cube._write_atom_lines()
        self.assertListEqual(actual_atom_line_list, expected_atom_line_list)

    def test_write_property_line_list(self):
        """Test to see if expected property lines are written.
        """
        expected_property_lines_list = [
            "   -13.010663    -1.077087    -0.178811    -0.034125",
            "   -13.010663    -1.077087    -0.013931    -0.030981",
            "   -12.845512    -2.061264    -0.508571    -0.055775",
            "   -12.845512    -1.897235     0.645590    -0.039007",
            "   -12.845512    -1.733205     0.810470    -0.032870",
            "   -12.845512    -1.405146     0.975350    -0.021576",
            "   -12.845512    -0.584999    -1.003211    -0.042233",
            "   -12.845512    -0.256940    -0.673451    -0.029992",
            "   -12.845512    -0.256940     0.975350     0.014346",
            "   -12.845512     0.071119     0.150949     0.000054",
        ]
        self.gaussian_cube.read(self.expected_filename)
        actual_property_lines_list = self.gaussian_cube._write_property_values_lines()
        self.assertListEqual(actual_property_lines_list, expected_property_lines_list)

    def test_write_file_lines(self):
        """Test to see if correct file lines are produced.
        """
        expected_file_lines = [
            " Marat",
            " Gaussian Cube file",
            "    3  -14.331872   -4.849766   -3.806173",
            "   61    0.165151    0.000000    0.000000",
            "   55    0.000000    0.164030    0.000000",
            "   54    0.000000    0.000000    0.164880",
            "    8    0.000000  -10.030764   -1.070314   -0.026721",
            "    1    0.000000   -8.202256   -1.001565    0.035731",
            "    1    0.000000  -10.552420    0.228376    1.153020",
            "   -13.010663    -1.077087    -0.178811    -0.034125",
            "   -13.010663    -1.077087    -0.013931    -0.030981",
            "   -12.845512    -2.061264    -0.508571    -0.055775",
            "   -12.845512    -1.897235     0.645590    -0.039007",
            "   -12.845512    -1.733205     0.810470    -0.032870",
            "   -12.845512    -1.405146     0.975350    -0.021576",
            "   -12.845512    -0.584999    -1.003211    -0.042233",
            "   -12.845512    -0.256940    -0.673451    -0.029992",
            "   -12.845512    -0.256940     0.975350     0.014346",
            "   -12.845512     0.071119     0.150949     0.000054",
        ]
        self.gaussian_cube.read(self.expected_filename)
        actual_file_lines = self.gaussian_cube._write_file_lines()
        self.assertListEqual(actual_file_lines, expected_file_lines)

    def test_write_file_line_string(self):
        """Test to see correct string is produced.
        """
        expected_file_string = """ Marat
 Gaussian Cube file
    3  -14.331872   -4.849766   -3.806173
   61    0.165151    0.000000    0.000000
   55    0.000000    0.164030    0.000000
   54    0.000000    0.000000    0.164880
    8    0.000000  -10.030764   -1.070314   -0.026721
    1    0.000000   -8.202256   -1.001565    0.035731
    1    0.000000  -10.552420    0.228376    1.153020
   -13.010663    -1.077087    -0.178811    -0.034125
   -13.010663    -1.077087    -0.013931    -0.030981
   -12.845512    -2.061264    -0.508571    -0.055775
   -12.845512    -1.897235     0.645590    -0.039007
   -12.845512    -1.733205     0.810470    -0.032870
   -12.845512    -1.405146     0.975350    -0.021576
   -12.845512    -0.584999    -1.003211    -0.042233
   -12.845512    -0.256940    -0.673451    -0.029992
   -12.845512    -0.256940     0.975350     0.014346
   -12.845512     0.071119     0.150949     0.000054
"""
        self.gaussian_cube.read(self.expected_filename)
        actual_file_string = self.gaussian_cube._write_file_lines_string()
        self.assertEqual(actual_file_string, expected_file_string)

    def test_write(self):
        """Test to see if correct file is written out.
        """
        expected_filename = self.expected_filename
        actual_filename = "actual_cube-0_000.cube"
        self.gaussian_cube.read(self.expected_filename)
        self.gaussian_cube.write(actual_filename)
        with open(expected_filename, "r") as expected_file:
            with open(actual_filename, "r") as actual_file:
                actual_file_read_in = actual_file.read()
                expected_file_read_in = expected_file.read()
                self.assertEqual(actual_file_read_in, expected_file_read_in)

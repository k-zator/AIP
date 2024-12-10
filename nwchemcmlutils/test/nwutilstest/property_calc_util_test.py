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
Script contains tests for the property  calc util scripts.

@author: mark
"""


import unittest
import logging
import unittest.mock as mock
import numpy as np
import nwchemcmlutils.nwUtils.property_calc_util as property_calc_util
import nwchemcmlutils.nwUtils.rtdb_util as rtdb_util

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.INFO)


class PropertyCalcUtilTestCase(unittest.TestCase):
    """Test case for the property_calc_util module.
    """

    def setUp(self):
        """set up for the Tests.
        """
        self.rot_axis_000 = np.array([0, 0, 0])
        self.epsiso_prop_dict_1 = {
            "padding": 2.0,
            "step_size": 0.088,
            "iso_surf": 0.002,
            "tol": 0.00003,
        }
        self.epsiso_prop_dict_2 = {
            "padding": 2.0,
            "step_size": 0.088,
            "iso_surf": 0.001,
            "tol": 0.000015,
        }
        self.n_voxel_values = [100, 101, 99, 102, 98, 100]
        self.get_values_1 = [
            rtdb_util.NWChemError,
            [
                -0.37416575,
                0.000000,
                0.65573492,
                1.43619175,
                0.000000,
                0.3855041,
                -1.06013628,
                0.000000,
                -1.04123902,
            ],
            "angstroms",
            1.8897259885789233,
        ]
        self.get_values_2 = [rtdb_util.NWChemError, "angstroms", 1.8897259885789233]
        self.get_values_4 = (
            self.get_values_1 + self.n_voxel_values[0:2] + self.get_values_2
        )
        get_rot_fragment1 = [
            rtdb_util.NWChemError,
            [
                -0.37416575,
                0.46367461,
                0.46367461,
                1.43619175,
                0.27259256,
                0.27259256,
                -1.06013628,
                -0.73626717,
                -0.73626717,
            ],
            "angstroms",
            1.8897259885789233,
        ]
        get_rot_fragment2 = [rtdb_util.NWChemError, "angstroms", 1.8897259885789233]
        self.get_values_5 = (
            self.get_values_1
            + self.get_values_2
            + self.n_voxel_values[0:2]
            + get_rot_fragment1
            + get_rot_fragment2
        )
        self.get_values_6 = self.get_values_5 + self.get_values_5
        self.get_values_7 = self.n_voxel_values[0:2] + self.get_values_5
        self.set_coords_values_1 = [
            "geometry",
            [
                -0.37416575,
                0.000000,
                0.65573492,
                1.43619175,
                0.000000,
                0.3855041,
                -1.06013628,
                0.000000,
                -1.04123902,
            ],
        ]

        self.angles_in_degrees = [45]
        self.axes_of_rotation = [np.array([0, 0, 0]), np.array([1, 0, 0])]

    def tearDown(self):
        """Tear down for tests.
        """
        del self.rot_axis_000
        del self.epsiso_prop_dict_1
        del self.epsiso_prop_dict_2
        del self.get_values_1
        del self.get_values_2
        del self.get_values_4
        del self.get_values_5
        del self.set_coords_values_1
        del self.angles_in_degrees
        del self.axes_of_rotation

    def test_create_cube_filename_stem(self):
        """Test to see if expected cube file stem is returned.
        """
        expected_cubefilename_stem = "water_0.0020_0_000"
        actual_cubefilename_stem = property_calc_util.createCubeFilenameStem(
            "water", 0.002, 0, self.rot_axis_000
        )
        self.assertEqual(actual_cubefilename_stem, expected_cubefilename_stem)

    def test_calc_property(self):
        """Test to see if the expected methods are called with the expected
        names.
        """
        task_energy_mock = mock.Mock()
        property_calc_util.nwchem.attach_mock(task_energy_mock, "task_energy")
        task_prop_mock = mock.Mock()
        property_calc_util.nwchem.attach_mock(task_prop_mock, "task_property")
        property_calc_util.calculateProperty("test_cube_stem", "dft")
        actual_call_arg_list_1 = [
            str(arg_list)
            for arg_list in property_calc_util.nwchem.task_energy.call_args_list
        ]
        expected_call_arg_list_1 = ["call('dft')"]
        self.assertListEqual(actual_call_arg_list_1, expected_call_arg_list_1)
        actual_call_arg_list_2 = [
            str(arg_list)
            for arg_list in property_calc_util.nwchem.task_property.call_args_list
        ]
        expected_call_arg_list_2 = ["call('dft')"]
        self.assertListEqual(actual_call_arg_list_2, expected_call_arg_list_2)

    def test_set_property_epsiso(self):
        """Test checks that expected call to nwchem input parse is made.
        """
        task_input_parse_mock = mock.Mock()
        property_calc_util.nwchem.attach_mock(task_input_parse_mock, "input_parse")
        property_calc_util.setPropertyEPSISO(2.0, 0.088, 0.002, 0.00003)
        actual_call_arg_list_1 = [
            str(arg_list)
            for arg_list in property_calc_util.nwchem.input_parse.call_args_list
        ]
        expected_call_arg_list_1 = [
            "call('\\n# Property module\\n# http://w"
            + "ww.nwchem-sw.org/index.php/Properties\\"
            + "n# New property keyword in use\\nproper"
            + "ty\\n  grid pad 2.000000 step 0.088000\\n"
            + "  espiso iso 0.002000 tol 0.000030\\n"
            + "end\\n')"
        ]
        self.assertListEqual(actual_call_arg_list_1, expected_call_arg_list_1)

    def test_set_and_calc_property_epsiso(self):
        """Test checks that expected arguements are used to call the input parse
        and task opt/prop methods.
        """
        # mock all the nwchem functions.
        task_input_parse_mock = mock.Mock()
        property_calc_util.nwchem.attach_mock(task_input_parse_mock, "input_parse")
        task_energy_mock = mock.Mock()
        property_calc_util.nwchem.attach_mock(task_energy_mock, "task_energy")
        task_prop_mock = mock.Mock()
        property_calc_util.nwchem.attach_mock(task_prop_mock, "task_property")
        rtdb_get_mock = mock.Mock(side_effect=self.n_voxel_values)
        property_calc_util.rtdb_util.nwchem.attach_mock(rtdb_get_mock, "rtdb_get")
        geom_set_coords_mock = mock.Mock(side_effect=self.set_coords_values_1)
        rtdb_util.nwgeom.attach_mock(geom_set_coords_mock, "geom_set_coords")
        # call the function
        n_voxels = property_calc_util.setAndCalculatePropertyEPSISO(
            "water_0.0020_0_000", "dft", self.epsiso_prop_dict_1
        )
        # start checking the arguments used to call the different functions.
        actual_call_arg_list_1 = [
            str(arg_list)
            for arg_list in property_calc_util.nwchem.input_parse.call_args_list
        ]
        expected_call_arg_list_1 = [
            "call('\\n# Property module\\n# http://w"
            + "ww.nwchem-sw.org/index.php/Properties\\"
            + "n# New property keyword in use\\nproper"
            + "ty\\n  grid pad 2.000000 step 0.088000\\n"
            + "  espiso iso 0.002000 tol 0.000030\\n"
            + "end\\n')"
        ]
        self.assertListEqual(actual_call_arg_list_1, expected_call_arg_list_1)
        actual_call_arg_list_2 = [
            str(arg_list)
            for arg_list in property_calc_util.nwchem.task_energy.call_args_list
        ]
        expected_call_arg_list_2 = ["call('dft')"]
        self.assertListEqual(actual_call_arg_list_2, expected_call_arg_list_2)
        actual_call_arg_list_3 = [
            str(arg_list)
            for arg_list in property_calc_util.nwchem.task_property.call_args_list
        ]
        expected_call_arg_list_3 = ["call('dft')"]
        self.assertListEqual(actual_call_arg_list_3, expected_call_arg_list_3)
        self.assertEqual(100, n_voxels)

    def test_set_and_calc_prop_epsiso_list(self):
        """Test checks that the expected arguments are used to call the input parse
        and task opt/prop methods.
        """
        # mock all the nwchem functions.
        task_input_parse_mock = mock.Mock()
        property_calc_util.nwchem.attach_mock(task_input_parse_mock, "input_parse")
        task_energy_mock = mock.Mock()
        property_calc_util.nwchem.attach_mock(task_energy_mock, "task_energy")
        task_prop_mock = mock.Mock()
        property_calc_util.nwchem.attach_mock(task_prop_mock, "task_property")
        rtdb_get_mock = mock.Mock(side_effect=self.n_voxel_values)
        property_calc_util.rtdb_util.nwchem.attach_mock(rtdb_get_mock, "rtdb_get")
        # call the function
        cube_file_stems = property_calc_util.setAndCalcPropEPSISOList(
            "water", "dft", [self.epsiso_prop_dict_1, self.epsiso_prop_dict_2]
        )
        # start checking the arguments used to call the different functions.
        actual_call_arg_list_1 = [
            str(arg_list)
            for arg_list in property_calc_util.nwchem.input_parse.call_args_list
        ]
        expected_call_arg_list_1 = [
            "call('\\n# Property module\\n# http://w"
            + "ww.nwchem-sw.org/index.php/Properties\\"
            + "n# New property keyword in use\\nproper"
            + "ty\\n  grid pad 2.000000 step 0.088000\\n"
            + "  espiso iso 0.002000 tol 0.000030\\n"
            + "end\\n')",
            "call('\\n# Property module\\n# http://w"
            + "ww.nwchem-sw.org/index.php/Properties\\"
            + "n# New property keyword in use\\nproper"
            + "ty\\n  grid pad 2.000000 step 0.088000\\n"
            + "  espiso iso 0.001000 tol 0.000015\\n"
            + "end\\n')",
        ]
        self.assertListEqual(actual_call_arg_list_1, expected_call_arg_list_1)
        actual_call_arg_list_2 = [
            str(arg_list)
            for arg_list in property_calc_util.nwchem.task_energy.call_args_list
        ]
        expected_call_arg_list_2 = ["call('dft')", "call('dft')"]
        self.assertListEqual(actual_call_arg_list_2, expected_call_arg_list_2)
        actual_call_arg_list_3 = [
            str(arg_list)
            for arg_list in property_calc_util.nwchem.task_property.call_args_list
        ]
        expected_call_arg_list_3 = ["call('dft')", "call('dft')"]
        self.assertListEqual(actual_call_arg_list_3, expected_call_arg_list_3)
        actual_cube_file_stems = cube_file_stems
        expected_cube_file_stems = {
            0.002: {"water_0.0020_0_000": 100},
            0.001: {"water_0.0010_0_000": 101},
        }
        self.assertDictEqual(expected_cube_file_stems, actual_cube_file_stems)

    def test_rot_set_calc_prop_epsiso_list_1(self):
        """Test to see if expected arguements are used to call the input_parse
        and task opt/prop methods.
        
        This is done for both the case where the axis is np.array([0, 0, 0])
        to check for correct error handling.
        """
        # mock all the nwchem functions.
        task_input_parse_mock = mock.Mock()
        property_calc_util.nwchem.attach_mock(task_input_parse_mock, "input_parse")
        task_energy_mock = mock.Mock()
        property_calc_util.nwchem.attach_mock(task_energy_mock, "task_energy")
        task_prop_mock = mock.Mock()
        property_calc_util.nwchem.attach_mock(task_prop_mock, "task_property")
        rtdb_get_mock = mock.Mock(side_effect=self.get_values_4)
        property_calc_util.rtdb_util.nwchem.attach_mock(rtdb_get_mock, "rtdb_get")
        geom_set_coords_mock = mock.Mock(side_effect=self.set_coords_values_1)
        rtdb_util.nwgeom.attach_mock(geom_set_coords_mock, "geom_set_coords")
        # call the function with non zero angle, and [0, 0, 0] for axis
        cube_file_stems = property_calc_util.rotateSetCalcPropEPSISOList(
            "water",
            "dft",
            [self.epsiso_prop_dict_1, self.epsiso_prop_dict_2],
            45,
            np.array([0, 0, 0]),
        )
        # start checking the arguments used to call the different functions.
        actual_call_arg_list_1 = [
            str(arg_list)
            for arg_list in property_calc_util.nwchem.input_parse.call_args_list
        ]
        expected_call_arg_list_1 = [
            "call('\\n# Property module\\n# http://w"
            + "ww.nwchem-sw.org/index.php/Properties\\"
            + "n# New property keyword in use\\nproper"
            + "ty\\n  grid pad 2.000000 step 0.088000\\n"
            + "  espiso iso 0.002000 tol 0.000030\\n"
            + "end\\n')",
            "call('\\n# Property module\\n# http://w"
            + "ww.nwchem-sw.org/index.php/Properties\\"
            + "n# New property keyword in use\\nproper"
            + "ty\\n  grid pad 2.000000 step 0.088000\\n"
            + "  espiso iso 0.001000 tol 0.000015\\n"
            + "end\\n')",
        ]
        self.assertListEqual(actual_call_arg_list_1, expected_call_arg_list_1)
        actual_call_arg_list_2 = [
            str(arg_list)
            for arg_list in property_calc_util.nwchem.task_energy.call_args_list
        ]
        expected_call_arg_list_2 = ["call('dft')", "call('dft')"]
        self.assertListEqual(actual_call_arg_list_2, expected_call_arg_list_2)
        actual_call_arg_list_3 = [
            str(arg_list)
            for arg_list in property_calc_util.nwchem.task_property.call_args_list
        ]
        expected_call_arg_list_3 = ["call('dft')", "call('dft')"]
        self.assertListEqual(actual_call_arg_list_3, expected_call_arg_list_3)
        actual_cube_file_stems = cube_file_stems
        expected_cube_file_stems = {
            0.001: {"water_0.0010_0_000": 101},
            0.002: {"water_0.0020_0_000": 100},
        }
        self.assertDictEqual(expected_cube_file_stems, actual_cube_file_stems)

    def test_rot_set_calc_prop_epsiso_list_2(self):
        """Test to see if expected arguements are used to call the input_parse
        and task opt/prop methods.
        
        This is done for both the case where the axis is np.array([1, 0, 0])
        to check for correct error handling.
        """
        # mock all the nwchem functions.
        task_input_parse_mock = mock.Mock()
        property_calc_util.nwchem.attach_mock(task_input_parse_mock, "input_parse")
        task_energy_mock = mock.Mock()
        property_calc_util.nwchem.attach_mock(task_energy_mock, "task_energye")
        task_prop_mock = mock.Mock()
        property_calc_util.nwchem.attach_mock(task_prop_mock, "task_property")
        rtdb_get_mock = mock.Mock(side_effect=self.get_values_5)
        property_calc_util.rtdb_util.nwchem.attach_mock(rtdb_get_mock, "rtdb_get")
        geom_set_coords_mock = mock.Mock()
        property_calc_util.rtdb_util.nwgeom.attach_mock(
            geom_set_coords_mock, "geom_set_coords"
        )
        # call the function
        cube_file_stems = property_calc_util.rotateSetCalcPropEPSISOList(
            "water",
            "dft",
            [self.epsiso_prop_dict_1, self.epsiso_prop_dict_2],
            45,
            np.array([1, 0, 0]),
        )
        # start checking the arguments used to call the different functions.
        actual_call_arg_list_1 = [
            str(arg_list)
            for arg_list in property_calc_util.nwchem.input_parse.call_args_list
        ]
        expected_call_arg_list_1 = [
            "call('\\n# Property module\\n# http://w"
            + "ww.nwchem-sw.org/index.php/Properties\\"
            + "n# New property keyword in use\\nproper"
            + "ty\\n  grid pad 2.000000 step 0.088000\\n"
            + "  espiso iso 0.002000 tol 0.000030\\n"
            + "end\\n')",
            "call('\\n# Property module\\n# http://w"
            + "ww.nwchem-sw.org/index.php/Properties\\"
            + "n# New property keyword in use\\nproper"
            + "ty\\n  grid pad 2.000000 step 0.088000\\n"
            + "  espiso iso 0.001000 tol 0.000015\\n"
            + "end\\n')",
        ]
        self.assertListEqual(actual_call_arg_list_1, expected_call_arg_list_1)
        actual_call_arg_list_2 = [
            str(arg_list)
            for arg_list in property_calc_util.nwchem.task_energy.call_args_list
        ]
        expected_call_arg_list_2 = [
            "call('dft')",
            "call('dft')",
            "call('dft')",
            "call('dft')",
        ]
        self.assertListEqual(actual_call_arg_list_2, expected_call_arg_list_2)
        actual_call_arg_list_3 = [
            str(arg_list)
            for arg_list in property_calc_util.nwchem.task_property.call_args_list
        ]
        expected_call_arg_list_3 = ["call('dft')", "call('dft')"]
        self.assertListEqual(actual_call_arg_list_3, expected_call_arg_list_3)
        actual_call_arg_list_4 = [
            arg_list for arg_list in rtdb_util.nwgeom.geom_set_coords.call_args_list
        ]
        expected_call_arg_list_4 = [
            (
                "geometry",
                [
                    -0.19800000225502176,
                    0.24536605380629931,
                    0.24536605380629928,
                    0.75999999930149575,
                    0.14424978273712952,
                    0.14424978273712949,
                    -0.56100000021549368,
                    -0.3896158365434288,
                    -0.38961583654342874,
                ],
            ),
            (
                "geometry",
                [
                    -0.19800000225502176,
                    -2.7755575615628914e-17,
                    0.34700000209192716,
                    0.75999999930149575,
                    -2.7755575615628914e-17,
                    0.20399999665766425,
                    -0.56100000021549368,
                    5.5511151231257827e-17,
                    -0.55099999874959149,
                ],
            ),
        ]
        for i, act_call_arg_list in enumerate(actual_call_arg_list_4):
            LOGGER.debug("act call list : %s", act_call_arg_list)
            for j, call_arg in enumerate(act_call_arg_list[0]):
                if isinstance(call_arg, list):
                    np.testing.assert_array_almost_equal(
                        np.array(expected_call_arg_list_4[i][j]), np.array(call_arg)
                    )
                else:
                    self.assertEqual(expected_call_arg_list_4[i][j], call_arg)
        actual_cube_file_stems = cube_file_stems
        expected_cube_file_stems = {
            0.001: {"water_0.0010_45_100": 101},
            0.002: {"water_0.0020_45_100": 100},
        }
        self.assertDictEqual(actual_cube_file_stems, expected_cube_file_stems)

    def test_rot_set_calc_prop_epsiso_list_3(self):
        """Test to see if expected arguements are used to call the input_parse
        and task opt/prop methods.
        
        This is done for both the case where the axis is np.array([1, 0, 0])
        and with angle = 0 to check for correct error handling.
        """
        # mock all the nwchem functions.
        task_input_parse_mock = mock.Mock()
        property_calc_util.nwchem.attach_mock(task_input_parse_mock, "input_parse")
        task_energy_mock = mock.Mock()
        property_calc_util.nwchem.attach_mock(task_energy_mock, "task_energy")
        task_prop_mock = mock.Mock()
        property_calc_util.nwchem.attach_mock(task_prop_mock, "task_property")
        rtdb_get_mock = mock.Mock(side_effect=self.get_values_6)
        property_calc_util.rtdb_util.nwchem.attach_mock(rtdb_get_mock, "rtdb_get")
        geom_set_coords_mock = mock.Mock(side_effect=self.set_coords_values_1)
        property_calc_util.rtdb_util.nwgeom.attach_mock(
            geom_set_coords_mock, "geom_set_coords"
        )
        # call the function
        cube_file_stems = property_calc_util.rotateSetCalcPropEPSISOList(
            "water",
            "dft",
            [self.epsiso_prop_dict_1, self.epsiso_prop_dict_2],
            0,
            np.array([1, 0, 0]),
        )
        # start checking the arguments used to call the different functions.
        actual_call_arg_list_1 = [
            str(arg_list)
            for arg_list in property_calc_util.nwchem.input_parse.call_args_list
        ]
        expected_call_arg_list_1 = [
            "call('\\n# Property module\\n# http://w"
            + "ww.nwchem-sw.org/index.php/Properties\\"
            + "n# New property keyword in use\\nproper"
            + "ty\\n  grid pad 2.000000 step 0.088000\\n"
            + "  espiso iso 0.002000 tol 0.000030\\n"
            + "end\\n')",
            "call('\\n# Property module\\n# http://w"
            + "ww.nwchem-sw.org/index.php/Properties\\"
            + "n# New property keyword in use\\nproper"
            + "ty\\n  grid pad 2.000000 step 0.088000\\n"
            + "  espiso iso 0.001000 tol 0.000015\\n"
            + "end\\n')",
        ]
        self.assertListEqual(actual_call_arg_list_1, expected_call_arg_list_1)
        actual_call_arg_list_2 = [
            str(arg_list)
            for arg_list in property_calc_util.nwchem.task_energy.call_args_list
        ]
        expected_call_arg_list_2 = ["call('dft')", "call('dft')"]
        self.assertListEqual(actual_call_arg_list_2, expected_call_arg_list_2)
        actual_call_arg_list_3 = [
            str(arg_list)
            for arg_list in property_calc_util.nwchem.task_property.call_args_list
        ]
        expected_call_arg_list_3 = ["call('dft')", "call('dft')"]
        self.assertListEqual(actual_call_arg_list_3, expected_call_arg_list_3)
        actual_cube_file_stems = cube_file_stems
        expected_cube_file_stems = {
            0.001: {"water_0.0010_0_000": 101},
            0.002: {"water_0.0020_0_000": 100},
        }
        self.assertDictEqual(actual_cube_file_stems, expected_cube_file_stems)

    def test_multi_orientation_calc_epsiso(self):
        """Test to see if expected arguments are used to call functions, and
        expected file stems are returned.
        """
        # mock all the nwchem functions.
        task_input_parse_mock = mock.Mock()
        property_calc_util.nwchem.attach_mock(task_input_parse_mock, "input_parse")
        task_energy_mock = mock.Mock()
        property_calc_util.nwchem.attach_mock(task_energy_mock, "task_energy")
        task_prop_mock = mock.Mock()
        property_calc_util.nwchem.attach_mock(task_prop_mock, "task_property")
        rtdb_get_mock = mock.Mock(side_effect=self.get_values_7)
        property_calc_util.rtdb_util.nwchem.attach_mock(rtdb_get_mock, "rtdb_get")
        geom_set_coords_mock = mock.Mock()
        property_calc_util.rtdb_util.nwgeom.attach_mock(
            geom_set_coords_mock, "geom_set_coords"
        )
        # now we do the call the function.
        cube_file_stems = property_calc_util.multiOrientationCalcPropEPSISOList(
            "water",
            "dft",
            [self.epsiso_prop_dict_1, self.epsiso_prop_dict_2],
            self.axes_of_rotation,
            self.angles_in_degrees,
        )
        actual_call_arg_list_1 = [
            str(arg_list)
            for arg_list in property_calc_util.nwchem.input_parse.call_args_list
        ]
        expected_call_arg_list_1 = [
            "call('\\n# Property module\\n# http://w"
            + "ww.nwchem-sw.org/index.php/Properties\\"
            + "n# New property keyword in use\\nproper"
            + "ty\\n  grid pad 2.000000 step 0.088000\\n"
            + "  espiso iso 0.002000 tol 0.000030\\n"
            + "end\\n')",
            "call('\\n# Property module\\n# http://w"
            + "ww.nwchem-sw.org/index.php/Properties\\"
            + "n# New property keyword in use\\nproper"
            + "ty\\n  grid pad 2.000000 step 0.088000\\n"
            + "  espiso iso 0.001000 tol 0.000015\\n"
            + "end\\n')",
            "call('\\n# Property module\\n# http://w"
            + "ww.nwchem-sw.org/index.php/Properties\\"
            + "n# New property keyword in use\\nproper"
            + "ty\\n  grid pad 2.000000 step 0.088000\\n"
            + "  espiso iso 0.002000 tol 0.000030\\n"
            + "end\\n')",
            "call('\\n# Property module\\n# http://w"
            + "ww.nwchem-sw.org/index.php/Properties\\"
            + "n# New property keyword in use\\nproper"
            + "ty\\n  grid pad 2.000000 step 0.088000\\n"
            + "  espiso iso 0.001000 tol 0.000015\\n"
            + "end\\n')",
        ]
        self.assertListEqual(actual_call_arg_list_1, expected_call_arg_list_1)
        actual_call_arg_list_2 = [
            str(arg_list)
            for arg_list in property_calc_util.nwchem.task_energy.call_args_list
        ]
        expected_call_arg_list_2 = [
            "call('dft')",
            "call('dft')",
            "call('dft')",
            "call('dft')",
        ]
        self.assertListEqual(actual_call_arg_list_2, expected_call_arg_list_2)
        actual_call_arg_list_3 = [
            str(arg_list)
            for arg_list in property_calc_util.nwchem.task_property.call_args_list
        ]
        expected_call_arg_list_3 = [
            "call('dft')",
            "call('dft')",
            "call('dft')",
            "call('dft')",
        ]
        self.assertListEqual(actual_call_arg_list_3, expected_call_arg_list_3)
        actual_call_arg_list_4 = [
            arg_list for arg_list in rtdb_util.nwgeom.geom_set_coords.call_args_list
        ]
        expected_call_arg_list_4 = [
            (
                "geometry",
                [
                    -0.19800000225502176,
                    0.24536605380629931,
                    0.24536605380629928,
                    0.75999999930149575,
                    0.14424978273712952,
                    0.14424978273712949,
                    -0.56100000021549368,
                    -0.3896158365434288,
                    -0.38961583654342874,
                ],
            ),
            (
                "geometry",
                [
                    -0.19800000225502176,
                    -2.7755575615628914e-17,
                    0.34700000209192716,
                    0.75999999930149575,
                    -2.7755575615628914e-17,
                    0.20399999665766425,
                    -0.56100000021549368,
                    5.5511151231257827e-17,
                    -0.55099999874959149,
                ],
            ),
        ]
        for i, act_call_arg_list in enumerate(actual_call_arg_list_4):
            LOGGER.debug("act call list : %s", act_call_arg_list)
            for j, call_arg in enumerate(act_call_arg_list[0]):
                if isinstance(call_arg, list):
                    np.testing.assert_array_almost_equal(
                        np.array(expected_call_arg_list_4[i][j]), np.array(call_arg)
                    )
                else:
                    self.assertEqual(expected_call_arg_list_4[i][j], call_arg)
        actual_cube_file_stems = cube_file_stems
        expected_cube_file_stems = {
            0.001: {"water_0.0010_0_000.cube": 101, "water_0.0010_45_100.cube": 101},
            0.002: {"water_0.0020_0_000.cube": 100, "water_0.0020_45_100.cube": 100},
        }
        self.assertDictEqual(actual_cube_file_stems, expected_cube_file_stems)

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
Script contains tests for the rtdb_util functions.

@author: mdd31
"""

import unittest
import logging
import unittest.mock as mock
import numpy as np
import nwchemcmlutils.nwUtils.rtdb_util as rtdb_util

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.DEBUG)


class RTDBUtilTestCase(unittest.TestCase):
    """Test case for the rtdb_util module.
    """

    def setUp(self):
        """set up for the Tests.
        """
        self.maxDiff = None
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
        self.get_values_3 = [rtdb_util.NWChemError, ["O", "H", "H"]]
        self.get_values_4 = self.get_values_1 + self.get_values_2
        self.get_values_5 = self.get_values_3 + self.get_values_1
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
        self.input_array = np.array(
            [
                [-0.198000, 0.000000, 0.347000],
                [0.760000, 0.000000, 0.204000],
                [-0.561000, 0.000000, -0.551000],
            ]
        )

    def tearDown(self):
        """tear down for the tests.
        """
        del self.get_values_1
        del self.get_values_2
        del self.get_values_3
        del self.get_values_4
        del self.set_coords_values_1
        del self.input_array

    def test_geom_get_coords(self):
        """Test to see if get coords returns the expected values.
        This requires mocking- run tests to see if correct functions are called.
        """
        rtdb_get_mock = mock.Mock(side_effect=self.get_values_1)
        rtdb_util.nwchem.attach_mock(rtdb_get_mock, "rtdb_get")
        actual_coords = rtdb_util.geomGetCoordsInAngToNumpyArray("geometry")
        rtdb_util.nwchem.rtdb_get.assert_called_with("geometry:geometry:angstrom_to_au")
        actual_call_arg_list = [
            str(arg_list) for arg_list in rtdb_util.nwchem.rtdb_get.call_args_list
        ]
        LOGGER.info(rtdb_util.nwchem.rtdb_get.call_args_list)
        expected_call_arg_list = [
            "call('geometry')",
            "call('geometry:geometry:coords')",
            "call('geometry:geometry:user units')",
            "call('geometry:geometry:angstrom_to_au')",
        ]
        self.assertListEqual(actual_call_arg_list, expected_call_arg_list)
        expected_coords = np.array(
            [
                [-0.198000, 0.000000, 0.347000],
                [0.760000, 0.000000, 0.204000],
                [-0.561000, 0.000000, -0.551000],
            ]
        )
        np.testing.assert_array_almost_equal(actual_coords, expected_coords)

    def test_geom_set_coords(self):
        """Test to see if set coords calls the expected methods with the
        expected arguments.
        """
        # set up the rtdb_set and rtdb_get mocks for this test.
        rtdb_get_mock = mock.Mock(side_effect=self.get_values_2)
        rtdb_util.nwchem.attach_mock(rtdb_get_mock, "rtdb_get")
        geom_set_coords_mock = mock.Mock(side_effect=self.set_coords_values_1)
        rtdb_util.nwgeom.attach_mock(geom_set_coords_mock, "geom_set_coords")
        # perform the setting of coordinates.
        rtdb_util.geomSetCoordsFromNumpyArray("geometry", self.input_array)
        # test to see if expected arguments for rtdb_get are used.
        actual_call_arg_list_1 = [
            str(arg_list) for arg_list in rtdb_util.nwchem.rtdb_get.call_args_list
        ]
        LOGGER.info(rtdb_util.nwchem.rtdb_get.call_args_list)
        expected_call_arg_list_1 = [
            "call('geometry')",
            "call('geometry:geometry:user units')",
            "call('geometry:geometry:angstrom_to_au')",
        ]
        self.assertListEqual(actual_call_arg_list_1, expected_call_arg_list_1)
        # test to see if expected arguements for rtdb_set are used.
        actual_call_arg_list_2 = [
            arg_list for arg_list in rtdb_util.nwgeom.geom_set_coords.call_args_list
        ]
        expected_call_arg_list_2 = [
            (
                "geometry",
                [
                    -0.19800000000000001,
                    0.0,
                    0.34699999999999998,
                    0.76000000000000001,
                    0.0,
                    0.20399999999999999,
                    -0.56100000000000005,
                    0.0,
                    -0.55100000000000005,
                ],
            )
        ]
        for i, act_call_arg_list in enumerate(actual_call_arg_list_2):
            LOGGER.debug("act call list : %s", act_call_arg_list)
            for j, call_arg in enumerate(act_call_arg_list[0]):
                if isinstance(call_arg, list):
                    np.testing.assert_array_almost_equal(
                        np.array(expected_call_arg_list_2[i][j]), np.array(call_arg)
                    )
                else:
                    self.assertEqual(expected_call_arg_list_2[i][j], call_arg)

    def test_get_atom_tags(self):
        """Test to see if the expected atom tags are returned and the expected
        calls to rtdb_util.nwchem.rtdb_get are made.
        """
        # set up the rtdb_get mock for this test.
        rtdb_get_mock = mock.Mock(side_effect=self.get_values_3)
        rtdb_util.nwchem.attach_mock(rtdb_get_mock, "rtdb_get")
        actual_atom_tags = rtdb_util.geomGetAtomTags("driverinitial")
        # test to see if expected arguments for rtdb_get are used.
        actual_call_arg_list_1 = [
            str(arg_list) for arg_list in rtdb_util.nwchem.rtdb_get.call_args_list
        ]
        LOGGER.info(rtdb_util.nwchem.rtdb_get.call_args_list)
        expected_call_arg_list_1 = [
            "call('driverinitial')",
            "call('geometry:driverinitial:tags')",
        ]
        self.assertListEqual(actual_call_arg_list_1, expected_call_arg_list_1)
        expected_atom_tags = ["O", "H", "H"]
        self.assertListEqual(actual_atom_tags, expected_atom_tags)

    def test_merge_tags_and_coords(self):
        """Test to see if the merging of the tags and coordinates produces
        expected dictionary.
        """
        expected_atom_dict = {
            0: {"coordinates": np.array([-0.198, 0.0, 0.347]), "element": "O"},
            1: {"coordinates": np.array([0.76, 0.0, 0.204]), "element": "H"},
            2: {"coordinates": np.array([-0.561, 0.0, -0.551]), "element": "H"},
        }
        atom_tags = ["O", "H", "H"]
        atom_coords = np.array(
            [
                [-0.198000, 0.000000, 0.347000],
                [0.760000, 0.000000, 0.204000],
                [-0.561000, 0.000000, -0.551000],
            ]
        )
        actual_atom_dict = rtdb_util.mergeTagAndCoords(atom_tags, atom_coords)
        self.assertListEqual(
            list(actual_atom_dict.keys()), list(expected_atom_dict.keys())
        )
        for key in actual_atom_dict.keys():
            self.assertEqual(
                actual_atom_dict[key]["element"], expected_atom_dict[key]["element"]
            )
            np.testing.assert_array_almost_equal(
                actual_atom_dict[key]["coordinates"],
                expected_atom_dict[key]["coordinates"],
            )

    def test_get_merge_tags_coords(self):
        """Test to see if the tags and coords are properly merged.
        """
        # Set up mock
        rtdb_get_mock = mock.Mock(side_effect=self.get_values_5)
        rtdb_util.nwchem.attach_mock(rtdb_get_mock, "rtdb_get")
        # now we do the reading, seeing if expected arguements are used.
        actual_atom_dict = rtdb_util.geomGetAndMergeTagsAndCoords(
            "driverinitial", "geometry"
        )
        # check args used.
        actual_call_arg_list_1 = [
            str(arg_list) for arg_list in rtdb_util.nwchem.rtdb_get.call_args_list
        ]
        expected_call_arg_list_1 = [
            "call('driverinitial')",
            "call('geometry:driverinitial:tags')",
            "call('geometry')",
            "call('geometry:geometry:coords')",
            "call('geometry:geometry:user units')",
            "call('geometry:geometry:angstrom_to_au')",
        ]
        self.assertListEqual(actual_call_arg_list_1, expected_call_arg_list_1)
        expected_atom_dict = {
            0: {"coordinates": np.array([-0.198, 0.0, 0.347]), "element": "O"},
            1: {"coordinates": np.array([0.76, 0.0, 0.204]), "element": "H"},
            2: {"coordinates": np.array([-0.561, 0.0, -0.551]), "element": "H"},
        }
        self.assertListEqual(
            list(actual_atom_dict.keys()), list(expected_atom_dict.keys())
        )
        for key in actual_atom_dict.keys():
            self.assertEqual(
                actual_atom_dict[key]["element"], expected_atom_dict[key]["element"]
            )
            np.testing.assert_array_almost_equal(
                actual_atom_dict[key]["coordinates"],
                expected_atom_dict[key]["coordinates"],
            )

    def test_rotate_molecule_about_centroid(self):
        """Test to see if the molecule is rotated about the centroid.
        """
        rtdb_get_mock = mock.Mock(side_effect=self.get_values_4)
        rtdb_util.nwchem.attach_mock(rtdb_get_mock, "rtdb_get")
        geom_set_coords_mock = mock.Mock(side_effect=self.set_coords_values_1)
        rtdb_util.nwchem.attach_mock(geom_set_coords_mock, "geom_set_coords")
        # now we do rotation.
        rtdb_util.rotateMoleculeAboutCentroid("geometry", 45, np.array([1, 0, 0]))
        actual_call_arg_list_1 = [
            str(arg_list) for arg_list in rtdb_util.nwchem.rtdb_get.call_args_list
        ]
        LOGGER.info(rtdb_util.nwchem.rtdb_get.call_args_list)
        expected_call_arg_list_1 = [
            "call('geometry')",
            "call('geometry:geometry:coords')",
            "call('geometry:geometry:user units')",
            "call('geometry:geometry:angstrom_to_au')",
            "call('geometry')",
            "call('geometry:geometry:user units')",
            "call('geometry:geometry:angstrom_to_au')",
        ]
        self.assertListEqual(actual_call_arg_list_1, expected_call_arg_list_1)
        # test to see if expected arguements for rtdb_set are used.
        actual_call_arg_list_2 = [
            arg_list for arg_list in rtdb_util.nwgeom.geom_set_coords.call_args_list
        ]
        expected_call_arg_list_2 = [
            (
                "geometry",
                [
                    -0.19800000000000001,
                    0.0,
                    0.34699999999999998,
                    0.76000000000000001,
                    0.0,
                    0.20399999999999999,
                    -0.56100000000000005,
                    0.0,
                    -0.55100000000000005,
                ],
            ),
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
        ]
        for i, act_call_arg_list in enumerate(actual_call_arg_list_2):
            LOGGER.debug("act call list 0: %s", act_call_arg_list[0])
            for j, call_arg in enumerate(act_call_arg_list[0]):
                if isinstance(call_arg, list):
                    np.testing.assert_array_almost_equal(
                        np.array(expected_call_arg_list_2[i][j]), np.array(call_arg)
                    )
                else:
                    self.assertEqual(expected_call_arg_list_2[i][j], call_arg)

    def test_rotate_mol_about_centroid_error(self):
        """Test to see if expected error is raised if we try to rotate about
        [0, 0, 0]
        """
        rtdb_get_mock = mock.Mock(side_effect=self.get_values_4)
        rtdb_util.nwchem.attach_mock(rtdb_get_mock, "rtdb_get")
        geom_set_coords_mock = mock.Mock(side_effect=self.set_coords_values_1)
        rtdb_util.nwchem.attach_mock(geom_set_coords_mock, "geom_set_coords")
        with self.assertRaises(ValueError) as err:
            rtdb_util.rotateMoleculeAboutCentroid("geometry", 0, np.array([0, 0, 0]))
        expected_args = "Can't rotate about [0, 0, 0]"
        actual_args = err.exception.args[0]
        self.assertEqual(actual_args, expected_args)

    def test_get_number_of_voxels_in_surface(self):
        """Test to see if expected number is returned
        """
        rtdb_get_mock = mock.Mock(side_effect=[100, rtdb_util.NWChemError])
        rtdb_util.nwchem.attach_mock(rtdb_get_mock, "rtdb_get")
        expected_number = 100
        self.assertEqual(expected_number, rtdb_util.get_number_of_voxels_in_surface())
        self.assertEqual(None, rtdb_util.get_number_of_voxels_in_surface())

    def test_set_cube_file_name(self):
        """Test to see if expected call arguements for writing the cube file
        names are called.
        """
        rtdb_put_mock = mock.Mock()
        rtdb_util.nwchem.attach_mock(rtdb_put_mock, "rtdb_put")
        rtdb_util.setCubeFileName("cubename.cube")
        actual_call_arg_list_1 = [
            str(arg_list) for arg_list in rtdb_util.nwchem.rtdb_put.call_args_list
        ]
        expected_call_arg_list_1 = ["call('prop:grid:output', 'cubename.cube')"]
        self.assertListEqual(actual_call_arg_list_1, expected_call_arg_list_1)

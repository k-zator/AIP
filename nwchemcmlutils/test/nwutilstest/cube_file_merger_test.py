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
Script contains tests for the cube_file_merger test case.

@author: mark
"""

import unittest
import logging
import os
import pathlib
import numpy as np
from nwchemcmlutils.DataClasses.GaussianCube import GaussianCube
from nwchemcmlutils.DataClasses.DistanceUnits import DistanceUnits
import nwchemcmlutils.nwUtils.cube_file_merger as cube_file_merger

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.WARN)


class CubeFileMergerTestCase(unittest.TestCase):
    """Test case for the property_calc_util module.
    """

    def setUp(self):
        """set up for the Tests.
        """
        self.maxDiff = None
        cube_file = GaussianCube()
        self.parent_directory = pathlib.Path(__file__).parents[1]
        cube_file.read(
            (self.parent_directory / "resources/full-0_000.cube").absolute().as_posix()
        )
        self.water_002_0_000_cube = cube_file
        cube_file_2 = cube_file_merger.readCubeFile(
            (
                self.parent_directory
                / "resources/cubefiles/notrot/pose_0.0020_-45_001.cube"
            )
            .absolute()
            .as_posix()
        )
        self.water_002_45_001_cube = cube_file_2
        cube_file_3 = cube_file_merger.readCubeFile(
            (
                self.parent_directory
                / "resources/cubefiles/rot/pose_0.0020_-45_001rot.cube"
            )
            .absolute()
            .as_posix()
        )
        self.water_002_r45_001_cube = cube_file_3
        self.axes_of_rotation = [
            np.array([0, 0, 0]),
            np.array([1, 0, 0]),
            np.array([0, 1, 0]),
            np.array([0, 0, 1]),
        ]
        self.angles_in_degrees = [0, -45, -120]
        cube_stem = (
            (self.parent_directory / "resources/cubefiles/notrot/pose")
            .absolute()
            .as_posix()
        )
        self.cube_filename_dict = cube_file_merger.createCubeFileNameDict(
            cube_stem, [0.002], self.angles_in_degrees, self.axes_of_rotation
        )
        cube_file_list = self.cube_filename_dict[0.002]
        self.cube_file_dict = cube_file_merger.readCubeFileList(cube_file_list)
        self.cube_stem = (
            (self.parent_directory / "resources/cubefiles/rot/pose")
            .absolute()
            .as_posix()
        )
        self.cube_stem_not_rot = (
            (self.parent_directory / "resources/cubefiles/notrot/pose")
            .absolute()
            .as_posix()
        )
        self.number_of_voxels = {
            self.cube_stem_not_rot + "_0.0020_-120_010.cube": 1000,
            self.cube_stem_not_rot + "_0.0020_-45_100.cube": 1001,
            self.cube_stem_not_rot + "_0.0020_-120_001.cube": 999,
            self.cube_stem_not_rot + "_0.0020_-45_001.cube": 1005,
            self.cube_stem_not_rot + "_0.0020_0_000.cube": 995,
            self.cube_stem_not_rot + "_0.0020_-45_010.cube": 1000,
            self.cube_stem_not_rot + "_0.0020_-120_100.cube": 1000,
        }

    def tearDown(self):
        """Tear down for tests.
        """
        del self.water_002_0_000_cube
        del self.axes_of_rotation
        del self.angles_in_degrees
        del self.cube_file_dict
        del self.water_002_45_001_cube
        del self.water_002_r45_001_cube
        del self.cube_filename_dict
        prod_file_list = [
            "resources/cubefiles/notrot/pose_0.0020_0_000rot.cube",
            "resources/cubefiles/notrot/pose_0.0020_-45_100rot.cube",
            "resources/cubefiles/notrot/pose_0.0020_-120_001rot.cube",
            "resources/cubefiles/notrot/pose_0.0020_-120_010rot.cube",
            "resources/cubefiles/notrot/pose_0.0020_-120_100rot.cube",
            "resources/cubefiles/notrot/pose_0.0020_-45_010rot.cube",
            "resources/cubefiles/notrot/pose_0.0020_-45_001rot.cube",
            "resources/cubefiles/notrot/pose_0.0020_merged.cube",
        ]
        prod_file_list = [
            (self.parent_directory / prod_file).absolute().as_posix()
            for prod_file in prod_file_list
        ]
        for prod_file in prod_file_list:
            if os.path.isfile(prod_file):
                os.remove(prod_file)

    def test_read_cube_file(self):
        """Test to see if the expected cube file 
        """
        expected_cube = self.water_002_0_000_cube
        actual_cube = cube_file_merger.readCubeFile(
            (self.parent_directory / "resources/full-0_000.cube").absolute().as_posix()
        )
        self.assertEqual(actual_cube, expected_cube)

    def test_rotate_cube_file(self):
        """Test to see if the expected rotated cube is returned.
        """
        expected_cube = self.water_002_r45_001_cube
        actual_cube = cube_file_merger.rotateCubeFile(
            self.water_002_45_001_cube, np.array([0, 0, 1]), 45
        )
        self.assertEqual(actual_cube, expected_cube)

    def test_read_cube_file_list(self):
        """Test to see if list of cube files is read correctly.
        """
        expected_dict = {
            (
                self.parent_directory
                / "resources/cubefiles/notrot/pose_0.0020_-45_001.cube"
            )
            .absolute()
            .as_posix(): self.water_002_45_001_cube
        }
        cube_filename_list = [
            (
                self.parent_directory
                / "resources/cubefiles/notrot/pose_0.0020_-45_001.cube"
            )
            .absolute()
            .as_posix()
        ]
        actual_dict = cube_file_merger.readCubeFileList(cube_filename_list)
        self.assertDictEqual(actual_dict, expected_dict)

    def test_create_cube_filename_dict(self):
        """Test to see if expected dictionary is returned.
        """
        expected_dict = {
            0.002: {
                "pose_0.0020_0_000.cube",
                "pose_0.0020_-45_100.cube",
                "pose_0.0020_-120_100.cube",
                "pose_0.0020_-45_010.cube",
                "pose_0.0020_-120_010.cube",
                "pose_0.0020_-45_001.cube",
                "pose_0.0020_-120_001.cube",
            }
        }
        acutal_dict = cube_file_merger.createCubeFileNameDict(
            "pose", [0.002], self.angles_in_degrees, self.axes_of_rotation
        )
        self.assertDictEqual(acutal_dict, expected_dict)

    def test_rotate_write_cube_file(self):
        """Test to see if cube file is rotated as expected.
        """
        expected_cube = self.water_002_r45_001_cube
        cube_file_merger.rotateWriteCubeFile(
            self.water_002_45_001_cube,
            (
                self.parent_directory
                / "resources/cubefiles/notrot/pose_0.0020_-45_001rot.cube"
            )
            .absolute()
            .as_posix(),
            np.array([0, 0, 1]),
            45,
        )
        actual_cube = cube_file_merger.readCubeFile(
            (
                self.parent_directory
                / "resources/cubefiles/notrot/pose_0.0020_-45_001rot.cube"
            )
            .absolute()
            .as_posix()
        )
        self.assertEqual(actual_cube, expected_cube)

    def test_rotate_back_write_dict(self):
        """Test to see if the cube files are correctly rotated.
        """
        exp_cube_stem = (
            (self.parent_directory / "resources/cubefiles/rot/pose")
            .absolute()
            .as_posix()
        )
        expected_cube_names = cube_file_merger.createCubeFileNameDict(
            exp_cube_stem, [0.002], self.angles_in_degrees, self.axes_of_rotation
        )
        expected_cube_names = expected_cube_names[0.002]
        expected_cube_names = [
            exp_cube.replace(".cube", "rot.cube") for exp_cube in expected_cube_names
        ]
        expected_cube_dict_a = cube_file_merger.readCubeFileList(expected_cube_names)
        expected_cube_dict = {}
        for key, cube_file in expected_cube_dict_a.items():
            new_key = key.replace(
                (self.parent_directory / "resources/cubefiles/rot/")
                .absolute()
                .as_posix(),
                "",
            )
            expected_cube_dict[new_key] = cube_file
        rotated_files = cube_file_merger.rotateBackandWriteCubeFileDict(
            self.cube_file_dict, 0.002, self.angles_in_degrees, self.axes_of_rotation
        )
        actual_cube_dict_a = cube_file_merger.readCubeFileList(rotated_files)
        actual_cube_dict = {}
        for key, cube_file in actual_cube_dict_a.items():
            new_key = key.replace(
                (self.parent_directory / "resources/cubefiles/notrot/")
                .absolute()
                .as_posix(),
                "",
            )
            actual_cube_dict[new_key] = cube_file
        self.assertDictEqual(actual_cube_dict, expected_cube_dict)
        exp_rot_file_list = sorted(
            [
                (
                    self.parent_directory
                    / "resources/cubefiles/notrot/pose_0.0020_0_000rot.cube"
                )
                .absolute()
                .as_posix(),
                (
                    self.parent_directory
                    / "resources/cubefiles/notrot/pose_0.0020_-45_100rot.cube"
                )
                .absolute()
                .as_posix(),
                (
                    self.parent_directory
                    / "resources/cubefiles/notrot/pose_0.0020_-120_001rot.cube"
                )
                .absolute()
                .as_posix(),
                (
                    self.parent_directory
                    / "resources/cubefiles/notrot/pose_0.0020_-120_010rot.cube"
                )
                .absolute()
                .as_posix(),
                (
                    self.parent_directory
                    / "resources/cubefiles/notrot/pose_0.0020_-120_100rot.cube"
                )
                .absolute()
                .as_posix(),
                (
                    self.parent_directory
                    / "resources/cubefiles/notrot/pose_0.0020_-45_010rot.cube"
                )
                .absolute()
                .as_posix(),
                (
                    self.parent_directory
                    / "resources/cubefiles/notrot/pose_0.0020_-45_001rot.cube"
                )
                .absolute()
                .as_posix(),
            ]
        )
        self.assertListEqual(sorted(rotated_files), exp_rot_file_list)

    def test_merge_cube_file_dict(self):
        """Test to see if expected cube file is created from merger.
        """
        cube_stem = (
            (self.parent_directory / "resources/cubefiles/rot/pose")
            .absolute()
            .as_posix()
        )
        cube_names = cube_file_merger.createCubeFileNameDict(
            cube_stem, [0.002], self.angles_in_degrees, self.axes_of_rotation
        )
        cube_names = cube_names[0.002]
        cube_names = [cube.replace(".cube", "rot.cube") for cube in cube_names]
        actual_cube, voxel_volumes = cube_file_merger.mergeCubeFileList(cube_names)
        expected_cube = cube_file_merger.readCubeFile(
            (self.parent_directory / "resources/cubefiles/rot/pose_0.0020_merged.cube")
            .absolute()
            .as_posix()
        )
        self.assertEqual(actual_cube, expected_cube)
        LOGGER.info("voxel volumes: %s", voxel_volumes)
        expected_volumes = {
            self.cube_stem
            + "_0.0020_-120_010": (0.0044550107345088007, DistanceUnits.atomic_units),
            self.cube_stem
            + "_0.0020_-45_100": (0.0044550107345088007, DistanceUnits.atomic_units),
            self.cube_stem
            + "_0.0020_-120_001": (0.0044550107345088007, DistanceUnits.atomic_units),
            self.cube_stem
            + "_0.0020_-45_001": (0.0044550107345088007, DistanceUnits.atomic_units),
            self.cube_stem
            + "_0.0020_0_000": (0.0044550107345088007, DistanceUnits.atomic_units),
            self.cube_stem
            + "_0.0020_-45_010": (0.0044550107345088007, DistanceUnits.atomic_units),
            self.cube_stem
            + "_0.0020_-120_100": (0.0044550107345088007, DistanceUnits.atomic_units),
        }
        self.assertListEqual(
            sorted(expected_volumes.keys()), sorted(voxel_volumes.keys())
        )
        for filename, voxel_volume in expected_volumes.items():
            self.assertAlmostEqual(voxel_volume[0], voxel_volumes[filename][0])
            self.assertEqual(voxel_volume[1], voxel_volumes[filename][1])

    def test_rotate_and_merge_dict(self):
        """Test to see if expected cube is produced.
        """
        actual_cube, voxel_volumes = cube_file_merger.rotateandMergeCubeFileDict(
            self.cube_file_dict, 0.002, self.angles_in_degrees, self.axes_of_rotation
        )
        expected_cube = cube_file_merger.readCubeFile(
            (self.parent_directory / "resources/cubefiles/rot/pose_0.0020_merged.cube")
            .absolute()
            .as_posix()
        )
        self.assertEqual(actual_cube, expected_cube)
        expected_volumes = {
            self.cube_stem_not_rot
            + "_0.0020_-120_010": (0.0044550107345088007, DistanceUnits.atomic_units),
            self.cube_stem_not_rot
            + "_0.0020_-45_100": (0.0044550107345088007, DistanceUnits.atomic_units),
            self.cube_stem_not_rot
            + "_0.0020_-120_001": (0.0044550107345088007, DistanceUnits.atomic_units),
            self.cube_stem_not_rot
            + "_0.0020_-45_001": (0.0044550107345088007, DistanceUnits.atomic_units),
            self.cube_stem_not_rot
            + "_0.0020_0_000": (0.0044550107345088007, DistanceUnits.atomic_units),
            self.cube_stem_not_rot
            + "_0.0020_-45_010": (0.0044550107345088007, DistanceUnits.atomic_units),
            self.cube_stem_not_rot
            + "_0.0020_-120_100": (0.0044550107345088007, DistanceUnits.atomic_units),
        }
        self.assertListEqual(
            sorted(expected_volumes.keys()), sorted(voxel_volumes.keys())
        )
        for filename, voxel_volume in expected_volumes.items():
            self.assertAlmostEqual(voxel_volume[0], voxel_volumes[filename][0])
            self.assertEqual(voxel_volume[1], voxel_volumes[filename][1])

    def test_calculate_mean_volume(self):
        """Test to see if expected value is returned.
        """
        expected_mean_volume = 4.455010734508801
        expected_volumes = {
            self.cube_stem_not_rot
            + "_0.0020_-120_010": (0.0044550107345088007, DistanceUnits.atomic_units),
            self.cube_stem_not_rot
            + "_0.0020_-45_100": (0.0044550107345088007, DistanceUnits.atomic_units),
            self.cube_stem_not_rot
            + "_0.0020_-120_001": (0.0044550107345088007, DistanceUnits.atomic_units),
            self.cube_stem_not_rot
            + "_0.0020_-45_001": (0.0044550107345088007, DistanceUnits.atomic_units),
            self.cube_stem_not_rot
            + "_0.0020_0_000": (0.0044550107345088007, DistanceUnits.atomic_units),
            self.cube_stem_not_rot
            + "_0.0020_-45_010": (0.0044550107345088007, DistanceUnits.atomic_units),
            self.cube_stem_not_rot
            + "_0.0020_-120_100": (0.0044550107345088007, DistanceUnits.atomic_units),
        }

        actual_mean_volume, unit = cube_file_merger.calculate_mean_volume(
            expected_volumes, self.number_of_voxels
        )
        self.assertEqual(DistanceUnits.atomic_units, unit)
        self.assertAlmostEqual(expected_mean_volume, actual_mean_volume)

    def test_read_rot_merge_write(self):
        """Test to see if files are correctly read, rotated, merged and the
        resulting file is written out.
        """
        expected_cube = cube_file_merger.readCubeFile(
            (
                self.parent_directory
                / "resources/cubefiles/rot/pose_0.0020_merged_withvol.cube"
            )
            .absolute()
            .as_posix()
        )
        cube_file_merger.readRotateMergeWriteCubeFileList(
            (self.parent_directory / "resources/cubefiles/notrot/pose")
            .absolute()
            .as_posix(),
            self.cube_filename_dict[0.002],
            0.002,
            self.angles_in_degrees,
            self.axes_of_rotation,
            self.number_of_voxels,
        )
        actual_cube = cube_file_merger.readCubeFile(
            (
                self.parent_directory
                / "resources/cubefiles/notrot/pose_0.0020_merged.cube"
            )
            .absolute()
            .as_posix()
        )
        self.assertEqual(actual_cube, expected_cube)

    def test_read_rot_merge_write_dict(self):
        """Test to see if files are correctly read, rotated, merged and
        resulting file are written out when looping over a list of iso surfaces.
        """
        expected_cube = cube_file_merger.readCubeFile(
            (
                self.parent_directory
                / "resources/cubefiles/rot/pose_0.0020_merged_withvol.cube"
            )
            .absolute()
            .as_posix()
        )
        cube_file_merger.readRotateMergeWriteCubeFileDict(
            (self.parent_directory / "resources/cubefiles/notrot/pose")
            .absolute()
            .as_posix(),
            {0.002: self.number_of_voxels},
            [0.002],
            self.angles_in_degrees,
            self.axes_of_rotation,
            test=True,
        )
        actual_cube = cube_file_merger.readCubeFile(
            (
                self.parent_directory
                / "resources/cubefiles/notrot/pose_0.0020_merged.cube"
            )
            .absolute()
            .as_posix()
        )
        self.assertEqual(actual_cube, expected_cube)

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
Script for tests of OptCalcCreator functions.

@author: mark
"""

import unittest
import logging
import os
import numpy as np
from lxml import etree
import pathlib
import nwchemcmlutils.calcCreator.CMLFileReader as CMLFileReader
import nwchemcmlutils.calcCreator.SubmitScriptCreator as SubmitScriptCreator

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.WARN)


class SubmitScriptCreatorTestCase(unittest.TestCase):
    """Test case for the OptCalcCreator module.
    """

    def setUp(self):
        """Set up for tests.
        """
        self.maxDiff = None
        self.parent_directory = pathlib.Path(__file__).parents[1]
        self.example_cml_filename = (
            (self.parent_directory / "resources/examplecml.cml").absolute().as_posix()
        )
        self.example_cml_file = CMLFileReader.readCMLFile(self.example_cml_filename)
        self.molecule_cml_list = CMLFileReader.extractMoleculeList(
            self.example_cml_file
        )
        self.molecule_cml = self.molecule_cml_list[0]
        self.molecule_cml_string = CMLFileReader.writeMoleculeCMLForNwinFile(
            self.molecule_cml
        )
        self.epsiso_prop_dict_list = [
            {"padding": 2.0, "step_size": 0.088, "iso_surf": 0.002, "tol": 0.00003},
            {"padding": 2.0, "step_size": 0.088, "iso_surf": 0.001, "tol": 0.000015},
        ]
        self.axes_of_rotation = [
            np.array([0, 0, 0]),
            np.array([1, 0, 0]),
            np.array([0, 1, 0]),
            np.array([0, 0, 1]),
        ]
        self.angles_in_degrees = [0, 45, 120]

    def tearDown(self):
        """Tear down afterr test.
        """
        del self.example_cml_file
        del self.molecule_cml_list
        del self.molecule_cml
        del self.molecule_cml_string
        del self.epsiso_prop_dict_list
        del self.axes_of_rotation
        del self.angles_in_degrees
        if os.path.isfile(
            "CALCS/XLYOFNOQVPJJNP-UHFFFAOYSA-N/XLYOFNOQVPJJNP-UHFFFAOYSA-N.nwin"
        ):
            os.remove(
                "CALCS/XLYOFNOQVPJJNP-UHFFFAOYSA-N/XLYOFNOQVPJJNP-UHFFFAOYSA-N.nwin"
            )
        if os.path.isfile(
            "CALCS/XLYOFNOQVPJJNP-UHFFFAOYSA-N/XLYOFNOQVPJJNP-UHFFFAOYSA-N_submit.sh"
        ):
            os.remove(
                "CALCS/XLYOFNOQVPJJNP-UHFFFAOYSA-N/XLYOFNOQVPJJNP-UHFFFAOYSA-N_submit.sh"
            )
        if os.path.isdir("CALCS/XLYOFNOQVPJJNP-UHFFFAOYSA-N/"):
            os.removedirs("CALCS/XLYOFNOQVPJJNP-UHFFFAOYSA-N/")
        if os.path.isfile("CALCS/job_submit.sh"):
            os.remove("CALCS/job_submit.sh")
        if os.path.isdir("CALCS/"):
            os.removedirs("CALCS/")

    def test_read_cml_file_to_mol_list(self):
        """Test to see if file is read in as expected.
        """
        expected_mol_list = self.molecule_cml_list
        actual_mol_list = SubmitScriptCreator.readCMLFileToMolList(
            self.example_cml_filename
        )
        self.assertEqual(len(actual_mol_list), len(expected_mol_list))
        for index, cml_mol in enumerate(expected_mol_list):
            exp_mol_string = etree.tounicode(cml_mol)
            act_mol_cml_string = etree.tounicode(actual_mol_list[index])
            self.assertMultiLineEqual(act_mol_cml_string, exp_mol_string)

    def test_read_cml_file_list_to_mol_list(self):
        """Test to see if a list of filenames is read correctly.
        """
        expected_mol_list = self.molecule_cml_list
        actual_mol_list = SubmitScriptCreator.readCMLFileListToMolList(
            [self.example_cml_filename]
        )
        self.assertEqual(len(actual_mol_list), len(expected_mol_list))
        for index, cml_mol in enumerate(expected_mol_list):
            exp_mol_string = etree.tounicode(cml_mol)
            act_mol_cml_string = etree.tounicode(actual_mol_list[index])
            self.assertMultiLineEqual(act_mol_cml_string, exp_mol_string)

    def test_create_reactor_dir(self):
        """Test to see if expected reactor directory is produced.
        """
        mk_dir = SubmitScriptCreator.createReactorDirectory("CALCS")
        dir_exists = os.path.isdir("CALCS")
        self.assertTrue(dir_exists)

    def test_create_react_dir_and_files(self):
        """Test to see if expected files are created.
        """
        expected_nwin_filename = (
            (self.parent_directory / "resources/expected_input3.nwin")
            .absolute()
            .as_posix()
        )
        actual_nwin_filename = (
            "CALCS/XLYOFNOQVPJJNP-UHFFFAOYSA-N/XLYOFNOQVPJJNP-UHFFFAOYSA-N.nwin"
        )
        expected_submit_filename = (
            (self.parent_directory / "resources/expected_input_submit.sh")
            .absolute()
            .as_posix()
        )
        actual_submit_filename = (
            "CALCS/XLYOFNOQVPJJNP-UHFFFAOYSA-N/XLYOFNOQVPJJNP-UHFFFAOYSA-N_submit.sh"
        )
        fail_list = SubmitScriptCreator.createReactorAndFiles(
            "CALCS", self.molecule_cml_list, 8, queue_system="TORQUE"
        )
        self.assertListEqual(fail_list, [])
        with open(expected_nwin_filename, "r") as expected_nwin_file:
            with open(actual_nwin_filename, "r") as actual_nwin_file:
                actual_nwin_file_read_in = actual_nwin_file.read()
                expected_nwin_file_read_in = expected_nwin_file.read()
                self.assertMultiLineEqual(
                    actual_nwin_file_read_in, expected_nwin_file_read_in
                )
        with open(expected_submit_filename, "r") as expected_submit_file:
            with open(actual_submit_filename, "r") as actual_submit_file:
                actual_submit_file_read_in = actual_submit_file.read()
                expected_submit_file_read_in = expected_submit_file.read()
                self.assertMultiLineEqual(
                    actual_submit_file_read_in, expected_submit_file_read_in
                )

    def test_create_react_dir_and_files_epsiso(self):
        """Test to see if expected files are created.
        """
        self.maxDiff = None
        expected_nwin_filename = (
            (self.parent_directory / "resources/expected_input2.nwin")
            .absolute()
            .as_posix()
        )
        actual_nwin_filename = (
            "CALCS/XLYOFNOQVPJJNP-UHFFFAOYSA-N/XLYOFNOQVPJJNP-UHFFFAOYSA-N.nwin"
        )
        expected_submit_filename = (
            (self.parent_directory / "resources/expected_input_submit_epsiso.sh")
            .absolute()
            .as_posix()
        )
        actual_submit_filename = (
            "CALCS/XLYOFNOQVPJJNP-UHFFFAOYSA-N/XLYOFNOQVPJJNP-UHFFFAOYSA-N_submit.sh"
        )
        epsiso_data = {
            "epsiso_param_dict_list": self.epsiso_prop_dict_list,
            "axes_of_rotation": self.axes_of_rotation,
            "angles_in_degrees": self.angles_in_degrees,
        }
        fail_list = SubmitScriptCreator.createReactorAndFiles(
            "CALCS",
            self.molecule_cml_list,
            8,
            queue_system="TORQUE",
            epsiso_calc=True,
            epsiso_data=epsiso_data,
        )
        self.assertListEqual(fail_list, [])
        with open(expected_nwin_filename, "r") as expected_nwin_file:
            with open(actual_nwin_filename, "r") as actual_nwin_file:
                actual_nwin_file_read_in = actual_nwin_file.read()
                expected_nwin_file_read_in = expected_nwin_file.read()
                self.assertMultiLineEqual(
                    actual_nwin_file_read_in, expected_nwin_file_read_in
                )
        with open(expected_submit_filename, "r") as expected_submit_file:
            with open(actual_submit_filename, "r") as actual_submit_file:
                actual_submit_file_read_in = actual_submit_file.read()
                expected_submit_file_read_in = expected_submit_file.read()
                self.assertMultiLineEqual(
                    actual_submit_file_read_in, expected_submit_file_read_in
                )

    def test_create_reactor_files_submit(self):
        """Test to see if the reactor and files are created along with submit
        script.
        """
        expected_nwin_filename = (
            (self.parent_directory / "resources/expected_input3.nwin")
            .absolute()
            .as_posix()
        )
        actual_nwin_filename = (
            "CALCS/XLYOFNOQVPJJNP-UHFFFAOYSA-N/XLYOFNOQVPJJNP-UHFFFAOYSA-N.nwin"
        )
        expected_submit_filename = (
            (self.parent_directory / "resources/expected_input_submit.sh")
            .absolute()
            .as_posix()
        )
        actual_submit_filename = (
            "CALCS/XLYOFNOQVPJJNP-UHFFFAOYSA-N/XLYOFNOQVPJJNP-UHFFFAOYSA-N_submit.sh"
        )
        expected_job_file = (
            (self.parent_directory / "resources/expected_reactor_submit.sh")
            .absolute()
            .as_posix()
        )
        actual_job_file = "CALCS/job_submit.sh"
        fail_list, sub_file = SubmitScriptCreator.createReactorFilesSubmit(
            "job_submit.sh", "CALCS", self.molecule_cml_list, 8, queue_system="TORQUE"
        )
        self.assertListEqual(fail_list, [])
        self.assertEqual(sub_file, 0)
        with open(expected_nwin_filename, "r") as expected_nwin_file:
            with open(actual_nwin_filename, "r") as actual_nwin_file:
                actual_nwin_file_read_in = actual_nwin_file.read()
                expected_nwin_file_read_in = expected_nwin_file.read()
                self.assertMultiLineEqual(
                    actual_nwin_file_read_in, expected_nwin_file_read_in
                )
        with open(expected_submit_filename, "r") as expected_submit_file:
            with open(actual_submit_filename, "r") as actual_submit_file:
                actual_submit_file_read_in = actual_submit_file.read()
                expected_submit_file_read_in = expected_submit_file.read()
                self.assertMultiLineEqual(
                    actual_submit_file_read_in, expected_submit_file_read_in
                )
        with open(expected_job_file, "r") as expected_job_file:
            with open(actual_job_file, "r") as actual_job_file:
                actual_job_file_read_in = actual_job_file.read()
                expected_job_file_read_in = expected_job_file.read()
                self.assertMultiLineEqual(
                    actual_job_file_read_in, expected_job_file_read_in
                )

    def test_create_jobs_from_cml(self):
        """Test to see if the reactor and files are created from input CML file
        list.
        """
        expected_nwin_filename = (
            (self.parent_directory / "resources/expected_input3.nwin")
            .absolute()
            .as_posix()
        )
        actual_nwin_filename = (
            "CALCS/XLYOFNOQVPJJNP-UHFFFAOYSA-N/XLYOFNOQVPJJNP-UHFFFAOYSA-N.nwin"
        )
        expected_submit_filename = (
            (self.parent_directory / "resources/expected_input_submit.sh")
            .absolute()
            .as_posix()
        )
        actual_submit_filename = (
            "CALCS/XLYOFNOQVPJJNP-UHFFFAOYSA-N/XLYOFNOQVPJJNP-UHFFFAOYSA-N_submit.sh"
        )
        expected_job_file = (
            (self.parent_directory / "resources/expected_reactor_submit.sh")
            .absolute()
            .as_posix()
        )
        actual_job_file = "CALCS/job_submit.sh"
        fail_list, sub_file = SubmitScriptCreator.createJobsFromCML(
            "job_submit.sh",
            "CALCS",
            [self.example_cml_filename],
            8,
            queue_system="TORQUE",
        )
        self.assertListEqual(fail_list, [])
        self.assertEqual(sub_file, 0)
        with open(expected_nwin_filename, "r") as expected_nwin_file:
            with open(actual_nwin_filename, "r") as actual_nwin_file:
                actual_nwin_file_read_in = actual_nwin_file.read()
                expected_nwin_file_read_in = expected_nwin_file.read()
                self.assertMultiLineEqual(
                    actual_nwin_file_read_in, expected_nwin_file_read_in
                )
        with open(expected_submit_filename, "r") as expected_submit_file:
            with open(actual_submit_filename, "r") as actual_submit_file:
                actual_submit_file_read_in = actual_submit_file.read()
                expected_submit_file_read_in = expected_submit_file.read()
                self.assertMultiLineEqual(
                    actual_submit_file_read_in, expected_submit_file_read_in
                )
        with open(expected_job_file, "r") as expected_job_file:
            with open(actual_job_file, "r") as actual_job_file:
                actual_job_file_read_in = actual_job_file.read()
                expected_job_file_read_in = expected_job_file.read()
                self.assertMultiLineEqual(
                    actual_job_file_read_in, expected_job_file_read_in
                )

    def test_create_jobs_from_cml_epsiso(self):
        """Test to see if the reactor and files are created from input CML file
        list.
        """
        expected_nwin_filename = (
            (self.parent_directory / "resources/expected_input2.nwin")
            .absolute()
            .as_posix()
        )
        actual_nwin_filename = (
            "CALCS/XLYOFNOQVPJJNP-UHFFFAOYSA-N/XLYOFNOQVPJJNP-UHFFFAOYSA-N.nwin"
        )
        expected_submit_filename = (
            (self.parent_directory / "resources/expected_input_submit_epsiso.sh")
            .absolute()
            .as_posix()
        )
        actual_submit_filename = (
            "CALCS/XLYOFNOQVPJJNP-UHFFFAOYSA-N/XLYOFNOQVPJJNP-UHFFFAOYSA-N_submit.sh"
        )
        expected_job_file = (
            (self.parent_directory / "resources/expected_reactor_submit.sh")
            .absolute()
            .as_posix()
        )
        actual_job_file = "CALCS/job_submit.sh"
        epsiso_data = {
            "epsiso_param_dict_list": self.epsiso_prop_dict_list,
            "axes_of_rotation": self.axes_of_rotation,
            "angles_in_degrees": self.angles_in_degrees,
        }
        fail_list, sub_file = SubmitScriptCreator.createJobsFromCML(
            "job_submit.sh",
            "CALCS",
            [self.example_cml_filename],
            8,
            queue_system="TORQUE",
            epsiso_calc=True,
            epsiso_data=epsiso_data,
        )
        self.assertListEqual(fail_list, [])
        self.assertEqual(sub_file, 0)
        with open(expected_nwin_filename, "r") as expected_nwin_file:
            with open(actual_nwin_filename, "r") as actual_nwin_file:
                actual_nwin_file_read_in = actual_nwin_file.read()
                expected_nwin_file_read_in = expected_nwin_file.read()
                self.assertMultiLineEqual(
                    actual_nwin_file_read_in, expected_nwin_file_read_in
                )
        with open(expected_submit_filename, "r") as expected_submit_file:
            with open(actual_submit_filename, "r") as actual_submit_file:
                actual_submit_file_read_in = actual_submit_file.read()
                expected_submit_file_read_in = expected_submit_file.read()
                self.assertMultiLineEqual(
                    actual_submit_file_read_in, expected_submit_file_read_in
                )
        with open(expected_job_file, "r") as expected_job_file:
            with open(actual_job_file, "r") as actual_job_file:
                actual_job_file_read_in = actual_job_file.read()
                expected_job_file_read_in = expected_job_file.read()
                self.assertMultiLineEqual(
                    actual_job_file_read_in, expected_job_file_read_in
                )

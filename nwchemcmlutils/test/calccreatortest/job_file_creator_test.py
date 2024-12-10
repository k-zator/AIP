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
import nwchemcmlutils.calcCreator.jobFileCreator as jobFileCreator

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.WARN)


class JobFileCreatorTestCase(unittest.TestCase):
    """Test case for the OptCalcCreator module.
    """

    def setUp(self):
        """Set up for tests.
        """
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
        if os.path.isfile("XLYOFNOQVPJJNP-UHFFFAOYSA-N.nwin"):
            os.remove("XLYOFNOQVPJJNP-UHFFFAOYSA-N.nwin")
        if os.path.isfile("XLYOFNOQVPJJNP-UHFFFAOYSA-N_submit.sh"):
            os.remove("XLYOFNOQVPJJNP-UHFFFAOYSA-N_submit.sh")

    def test_write_nwin_file(self):
        """Test to see if expected file is written.
        """
        self.maxDiff = None
        expected_filename = (
            (self.parent_directory / "resources/expected_input3.nwin")
            .absolute()
            .as_posix()
        )
        actual_filename = "XLYOFNOQVPJJNP-UHFFFAOYSA-N.nwin"
        jobFileCreator.writeNwinFile(
            molecule_cml=self.molecule_cml, memory_limit=8, directory=""
        )
        with open(expected_filename, "r") as expected_file:
            with open(actual_filename, "r") as actual_file:
                actual_file_read_in = actual_file.read()
                expected_file_read_in = expected_file.read()
                self.assertMultiLineEqual(actual_file_read_in, expected_file_read_in)

    def test_write_slurm_submit_file(self):
        """Function tests the writing of submit file.
        """
        self.maxDiff = None
        expected_filename = (
            (self.parent_directory / "resources/expected_slurm3_submit.sh")
            .absolute()
            .as_posix()
        )
        actual_filename = "XLYOFNOQVPJJNP-UHFFFAOYSA-N_submit.sh"
        jobFileCreator.writeSLURMSubmitFile(
            molecule_cml=self.molecule_cml, directory=""
        )
        with open(expected_filename, "r") as expected_file:
            with open(actual_filename, "r") as actual_file:
                actual_file_read_in = actual_file.read()
                expected_file_read_in = expected_file.read()
                self.assertMultiLineEqual(expected_file_read_in, actual_file_read_in)

    def test_write_ziggy_submit_file(self):
        """Function tests the writing of submit file.
        """
        expected_filename = (
            (self.parent_directory / "resources/expected_input_submit.sh")
            .absolute()
            .as_posix()
        )
        actual_filename = "XLYOFNOQVPJJNP-UHFFFAOYSA-N_submit.sh"
        jobFileCreator.writeZiggySubmitFile(
            molecule_cml=self.molecule_cml, directory=""
        )
        with open(expected_filename, "r") as expected_file:
            with open(actual_filename, "r") as actual_file:
                actual_file_read_in = actual_file.read()
                expected_file_read_in = expected_file.read()
                self.assertMultiLineEqual(actual_file_read_in, expected_file_read_in)

    def test_write_ziggy_submit_file_epsiso(self):
        """Test to see if expected submit file is written when the job is an
        epsiso calculation.
        """
        expected_filename_epsiso = (
            (self.parent_directory / "resources/expected_input_submit_epsiso.sh")
            .absolute()
            .as_posix()
        )
        actual_filename_epsiso = "XLYOFNOQVPJJNP-UHFFFAOYSA-N_submit.sh"
        jobFileCreator.writeZiggySubmitFile(
            molecule_cml=self.molecule_cml, directory="", epsiso_calc=True
        )
        with open(expected_filename_epsiso, "r") as expected_file:
            with open(actual_filename_epsiso, "r") as actual_file:
                actual_file_read_in = actual_file.read()
                expected_file_read_in = expected_file.read()
                self.assertMultiLineEqual(actual_file_read_in, expected_file_read_in)

    def test_writenwin_and_slurm_submit_file(self):
        """Test to see if expected file is produced.
        """
        expected_nwin_filename = (
            (self.parent_directory / "resources/expected_input3.nwin")
            .absolute()
            .as_posix()
        )
        actual_nwin_filename = "XLYOFNOQVPJJNP-UHFFFAOYSA-N.nwin"
        expected_submit_filename = (
            (self.parent_directory / "resources/expected_slurm3_submit.sh")
            .absolute()
            .as_posix()
        )
        actual_submit_filename = "XLYOFNOQVPJJNP-UHFFFAOYSA-N_submit.sh"
        jobFileCreator.writeNwinAndSubmitFiles(
            molecule_cml=self.molecule_cml, memory_limit=8, directory=""
        )
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

    def test_writenwin_and_ziggy_submit_file(self):
        """Test to see both files are written as expected.
        """
        expected_nwin_filename = (
            (self.parent_directory / "resources/expected_input3.nwin")
            .absolute()
            .as_posix()
        )
        actual_nwin_filename = "XLYOFNOQVPJJNP-UHFFFAOYSA-N.nwin"
        expected_submit_filename = (
            (self.parent_directory / "resources/expected_input_submit.sh")
            .absolute()
            .as_posix()
        )
        actual_submit_filename = "XLYOFNOQVPJJNP-UHFFFAOYSA-N_submit.sh"
        jobFileCreator.writeNwinAndSubmitFiles(
            molecule_cml=self.molecule_cml,
            memory_limit=8,
            queue_system="TORQUE",
            directory="",
        )
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

    def test_writenwin_and_ziggy_submit_file_epsiso(self):
        """Test to see both files are written as expected.
        """
        self.maxDiff = None
        expected_nwin_filename = (
            (self.parent_directory / "resources/expected_input2.nwin")
            .absolute()
            .as_posix()
        )
        actual_nwin_filename = "XLYOFNOQVPJJNP-UHFFFAOYSA-N.nwin"
        expected_submit_filename = (
            (self.parent_directory / "resources/expected_input_submit_epsiso.sh")
            .absolute()
            .as_posix()
        )
        actual_submit_filename = "XLYOFNOQVPJJNP-UHFFFAOYSA-N_submit.sh"
        epsiso_data = {
            "epsiso_param_dict_list": self.epsiso_prop_dict_list,
            "axes_of_rotation": self.axes_of_rotation,
            "angles_in_degrees": self.angles_in_degrees,
        }
        jobFileCreator.writeNwinAndSubmitFiles(
            molecule_cml=self.molecule_cml,
            memory_limit=8,
            queue_system="TORQUE",
            directory="",
            epsiso_calc=True,
            epsiso_data=epsiso_data,
        )
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

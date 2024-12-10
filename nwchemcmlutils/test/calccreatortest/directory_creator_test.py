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
Script for tests of DirectorycCreator functions.

@author: mark
"""

import unittest
import logging
import os
from lxml import etree
import pathlib
import nwchemcmlutils.calcCreator.CMLFileReader as CMLFileReader
import nwchemcmlutils.calcCreator.DirectoryCreator as DirectoryCreator

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.WARN)


class DirectoryCreatorTestCase(unittest.TestCase):
    """Test case for the DirectorycCreator module.
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

    def tearDown(self):
        """Tear down afterr test.
        """
        del self.example_cml_file
        del self.molecule_cml_list
        del self.molecule_cml
        del self.molecule_cml_string
        if os.path.isfile(
            "XLYOFNOQVPJJNP-UHFFFAOYSA-N/XLYOFNOQVPJJNP-UHFFFAOYSA-N.nwin"
        ):
            os.remove("XLYOFNOQVPJJNP-UHFFFAOYSA-N/XLYOFNOQVPJJNP-UHFFFAOYSA-N.nwin")
        if os.path.isfile(
            "XLYOFNOQVPJJNP-UHFFFAOYSA-N/XLYOFNOQVPJJNP-UHFFFAOYSA-N_submit.sh"
        ):
            os.remove(
                "XLYOFNOQVPJJNP-UHFFFAOYSA-N/XLYOFNOQVPJJNP-UHFFFAOYSA-N_submit.sh"
            )
        if os.path.isdir("XLYOFNOQVPJJNP-UHFFFAOYSA-N"):
            os.removedirs("XLYOFNOQVPJJNP-UHFFFAOYSA-N")
        if os.path.isfile("actual_bash_file.sh"):
            os.remove("actual_bash_file.sh")

    def test_create_directory(self):
        """Test to see if directory is created as expected, also check to see
        if expected error is raised when the directory already exists.
        """
        expected_directory_name = "XLYOFNOQVPJJNP-UHFFFAOYSA-N"
        dir_made = DirectoryCreator.createDirectory("XLYOFNOQVPJJNP-UHFFFAOYSA-N")
        is_dir = os.path.isdir(expected_directory_name)
        self.assertTrue(is_dir)
        with self.assertRaises(ValueError) as err:
            dir_made = DirectoryCreator.createDirectory(
                (self.parent_directory / "resources/cubefiles").absolute().as_posix()
            )
        expected_args = "Directory already exists."
        actual_args = err.exception.args[0]
        self.assertEqual(actual_args, expected_args)

    def test_create_directory_cml(self):
        """Test to see if directory is created.
        """
        expected_directory_name = "XLYOFNOQVPJJNP-UHFFFAOYSA-N"
        dir_made = DirectoryCreator.createDirectoryCML(self.molecule_cml, "")
        is_dir = os.path.isdir(expected_directory_name)
        self.assertTrue(is_dir)

    def test_mkdir_w_nwin_and_ziggy_submit_file(self):
        """Test to see both files are written as expected.
        """
        self.maxDiff = None
        expected_nwin_filename = (
            (self.parent_directory / "resources/expected_input3.nwin")
            .absolute()
            .as_posix()
        )
        actual_nwin_filename = (
            "XLYOFNOQVPJJNP-UHFFFAOYSA-N/XLYOFNOQVPJJNP-UHFFFAOYSA-N.nwin"
        )
        expected_submit_filename = (
            (self.parent_directory / "resources/expected_input_submit.sh")
            .absolute()
            .as_posix()
        )
        actual_submit_filename = (
            "XLYOFNOQVPJJNP-UHFFFAOYSA-N/XLYOFNOQVPJJNP-UHFFFAOYSA-N_submit.sh"
        )
        DirectoryCreator.createDirectoryAndFiles(
            self.molecule_cml, 8, "", queue_system="TORQUE"
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
        with self.assertRaises(ValueError) as err:
            DirectoryCreator.createDirectoryAndFiles(self.molecule_cml, 8, "")
        expected_args = "Directory already exists."
        actual_args = err.exception.args[0]
        self.assertEqual(actual_args, expected_args)

    def test_mkdirs_and_files(self):
        """Test to see if files can be genreated from a list of molecule cml
        elements.
        """
        self.maxDiff = None
        expected_nwin_filename = (
            (self.parent_directory / "resources/expected_input3.nwin")
            .absolute()
            .as_posix()
        )
        actual_nwin_filename = (
            "XLYOFNOQVPJJNP-UHFFFAOYSA-N/XLYOFNOQVPJJNP-UHFFFAOYSA-N.nwin"
        )
        expected_submit_filename = (
            (self.parent_directory / "resources/expected_input_submit.sh")
            .absolute()
            .as_posix()
        )
        actual_submit_filename = (
            "XLYOFNOQVPJJNP-UHFFFAOYSA-N/XLYOFNOQVPJJNP-UHFFFAOYSA-N_submit.sh"
        )
        fail_list = DirectoryCreator.createDirectoriesAndFiles(
            self.molecule_cml_list, 8, "", queue_system="TORQUE"
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
        new_fail_list = DirectoryCreator.createDirectoriesAndFiles(
            self.molecule_cml_list, 8, ""
        )
        self.assertListEqual(new_fail_list, ["XLYOFNOQVPJJNP-UHFFFAOYSA-N"])

    def test_write_bash_script(self):
        """Test to see if expected bash script is produced.
        """
        expected_bash_filename = (
            (self.parent_directory / "resources/expected_fullSubmit.sh")
            .absolute()
            .as_posix()
        )
        actual_bash_filename = "actual_bash_file.sh"
        DirectoryCreator.createDirectoriesAndFiles(self.molecule_cml_list, 8, "")
        DirectoryCreator.writeBashSubmitfile(actual_bash_filename, "")
        with open(expected_bash_filename, "r") as expected_bash_file:
            with open(actual_bash_filename, "r") as actual_bash_file:
                expected_bash_file_read_in = expected_bash_file.read()
                actual_bash_file_read_in = actual_bash_file.read()
                self.assertMultiLineEqual(
                    actual_bash_file_read_in, expected_bash_file_read_in
                )

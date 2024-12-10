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
Script for cleaning up files after calculation of merged cube file.

@author: mark
"""

import logging
import glob
import os

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.WARN)


def clean_up_files(molecule_name, expected_number_merged):
    """This cleans up all the files, if the correct numer of merged files
    have been found.
    """
    merged_file_list = glob.glob("{0}_*_merged.cube".format(molecule_name))
    if expected_number_merged == len(merged_file_list):
        delete_calc_files(molecule_name)
        delete_initial_cube_files(molecule_name)


def delete_initial_cube_files(molecule_name):
    """This dleetes the cube files
    """
    cube_files = glob.glob("{0}_*_*_*.cube".format(molecule_name))
    for filename in cube_files:
        delete_file(filename)


def delete_calc_files(molecule_name):
    """This cleans up the calculation files.
    """
    filenames = [
        "nwchem.py*",
        "{0}*.hess".format(molecule_name),
        "{0}*.movecs".format(molecule_name),
        "{0}*.db".format(molecule_name),
    ]
    for filename in filenames:
        delete_file(filename)


def delete_file(filename):
    """This trys to delete the given file.
    """
    if os.path.isfile(filename):
        os.remove(filename)

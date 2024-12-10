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
Script to create the directories

@author: mark
"""

import logging
import os
import glob
import nwchemcmlutils.calcCreator.CMLFileReader as CMLFileReader
import nwchemcmlutils.calcCreator.jobFileCreator as jobFileCreator

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.WARN)


def createDirectory(directory_name):
    """Creates a directory with name directory_name if it doesn't exist,
    if it exists it raises a ValueError.
    """
    if not os.path.isdir(directory_name):
        directory = os.mkdir(directory_name)
        return directory
    else:
        LOGGER.warning("Directory already exists")
        raise ValueError("Directory already exists.")


def createDirectoryCML(molecule_cml, calc_dir):
    """creates a directory with name equal to the StdInChiKey
    """
    if calc_dir != "":
        dir_prefix = calc_dir + "/"
    else:
        dir_prefix = ""
    stdinchikey = CMLFileReader.extractStdInChiKey(molecule_cml)
    LOGGER.debug("StdInChIKey: %s", stdinchikey)
    directory_name = dir_prefix + stdinchikey
    directory = createDirectory(directory_name)
    return directory_name


def createDirectoryAndFiles(molecule_cml, memory_limit, calc_dir, **kwargs):
    """Creates a directory with name equal to the StdInChiKey, then creates an
    nwin and submit script inside the directory.
    """
    try:
        directory_name = createDirectoryCML(molecule_cml, calc_dir) + "/"
    except ValueError as err:
        LOGGER.warning("Directory already exists.")
        LOGGER.warning("We do not want to overwrite files. So raising error.")
        raise err
    else:
        files_out = jobFileCreator.writeNwinAndSubmitFiles(
            molecule_cml=molecule_cml,
            memory_limit=memory_limit,
            directory=directory_name,
            **kwargs
        )
        return files_out


def createDirectoriesAndFiles(molecule_cml_list, memory_limit, calc_dir, **kwargs):
    """Creates directory and nwin and submit script for all files in the list.
    """
    failed_file_creation_list = []
    for molecule_cml in molecule_cml_list:
        try:
            createDirectoryAndFiles(
                molecule_cml=molecule_cml,
                memory_limit=memory_limit,
                calc_dir=calc_dir,
                **kwargs
            )
        except ValueError:
            directory_name = CMLFileReader.extractStdInChiKey(molecule_cml)
            failed_file_creation_list.append(directory_name)
    return failed_file_creation_list


def writeBashSubmitfile(filename, calc_dir):
    """function creates a bash script to be run, which will submit all submit
    scripts in the lower directories.
    """
    if calc_dir != "":
        dir_prefix = calc_dir + "/"
    else:
        dir_prefix = ""
    submit_file_list = glob.glob(dir_prefix + "*/*_submit.sh")
    submit_file_list = sorted(submit_file_list)
    with open(dir_prefix + filename, "w") as bash_file:
        LOGGER.debug("Write first line.")
        bash_file.write("#!/bin/bash\n")
        for submit_file_w_dir in submit_file_list:
            # we want to first change to the directory contatining the file,
            # then submit the script, then exit back.
            submit_directory, submit_file = submit_file_w_dir.split("/")[-2:]
            submit_line = "cd {0} && sbatch {1} && cd ../\n".format(
                submit_directory, submit_file
            )
            bash_file.write(submit_line)
            bash_file.write("sleep 1\n")
    return 0

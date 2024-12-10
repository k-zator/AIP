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

Script provides interface for job file creation- parameters are fed to it
and then the corresponding jobs are created.

@author: mark
"""

import logging
import copy
from lxml import etree
import nwchemcmlutils.calcCreator.CMLFileReader as CMLFileReader
import nwchemcmlutils.calcCreator.DirectoryCreator as DirectoryCreator

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.DEBUG)

CML_SCHEMA = CMLFileReader.createCMLSchema(
    "http://www-hunter.ch.cam.ac.uk/schema/cmlschema_AIP.xsd"
)


def readCMLFileToMolList(filename, **kwargs):
    """Function reads in a cml file and creates a molecule list.
    """
    LOGGER.info("Reading CML file.")
    LOGGER.info(filename)
    molecule_cml_file = CMLFileReader.readCMLFile(filename)
    LOGGER.info("Creating molecule list")
    cml_schema = kwargs.get("cmlschema", CML_SCHEMA)
    #try:
    CMLFileReader.validateXMLFile(copy.deepcopy(molecule_cml_file), cml_schema)
    #except etree.DocumentInvalid as e:
    #LOGGER.error(e)
    #return None
    molecule_cml_list = CMLFileReader.extractMoleculeList(molecule_cml_file)
    LOGGER.info("Created molecule list.")
    LOGGER.debug("Found %i cml molecules", len(molecule_cml_list))
    return molecule_cml_list


def readCMLFileListToMolList(filename_list, **kwargs):
    """Function reads a list of cml files, and then extracts all the molecules
    in them, returning one catenated list.
    """
    LOGGER.info("Reading in each file, and creating a list")
    merged_molecule_cml_list = []
    for filename in filename_list:
        LOGGER.debug("Reading CML file")
        mol_cml_list = readCMLFileToMolList(filename, **kwargs)
        LOGGER.debug("Adding each molecule cml to merged list")
        for mol_cml in mol_cml_list:
            merged_molecule_cml_list.append(mol_cml)
    LOGGER.info("Created merged molecule list")
    return merged_molecule_cml_list


def createReactorDirectory(reactor_name):
    """Creates directory within which all sub directories are created.
    """
    mk_dir = DirectoryCreator.createDirectory(reactor_name + "/")
    return mk_dir


def createReactorAndFiles(
    reactor_name, molecule_cml_list, memory_limit, *import_args, **kwargs
):
    """Function creates the reactor directory and also all the sub directories for the jobs,
    based on the arguements given.
    """
    reactor_dir = createReactorDirectory(reactor_name)
    failed_files = DirectoryCreator.createDirectoriesAndFiles(
        molecule_cml_list, memory_limit, reactor_name, *import_args, **kwargs
    )
    return failed_files


def createReactorFilesSubmit(
    submit_filename,
    reactor_name,
    molecule_cml_list,
    memory_limit,
    *import_args,
    **kwargs
):
    """Function creates the reactor directory and the files for submission, as
    well as the script to submit all the job files.
    """
    failed_files = createReactorAndFiles(
        reactor_name, molecule_cml_list, memory_limit, *import_args, **kwargs
    )
    submit_file = DirectoryCreator.writeBashSubmitfile(submit_filename, reactor_name)
    return failed_files, submit_file


def createJobsFromCML(
    submit_filename, reactor_name, cml_file_list, memory_limit, *import_args, **kwargs
):
    """Function reads in the cml from the given cml files, then creates the
    reactor directory, and the directories containing all the calculations,
    and the submit file for all the jobs.
    """
    molecule_cml_list = readCMLFileListToMolList(cml_file_list, **kwargs)
    failed_files, submit_file = createReactorFilesSubmit(
        submit_filename,
        reactor_name,
        molecule_cml_list,
        memory_limit,
        *import_args,
        **kwargs
    )
    return failed_files, submit_file

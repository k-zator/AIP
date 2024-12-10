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
Script for generating the job files. This contains the default values for
surfaces, so all that need to be supplied is the CML files.

@author: mark
"""

import logging
import glob
import numpy as np
import os
import pathlib
import json
import nwchemcmlutils.calcCreator.SubmitScriptCreator as SubmitCreator

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.WARN)

# Default values for Calculations

epsiso_prop_dict_list_one_surf = [
    {"padding": 4.0, "step_size": 0.088, "iso_surf": 0.002, "tol": 0.00003}

]

epsiso_prop_dict_list_tri_surf = [
        {'padding': 4.0, 'iso_surf': 0.0104, 'step_size': 0.04, 'tol': 3e-05},  
        {'padding': 4.0, 'iso_surf': 0.0300, 'step_size': 0.04, 'tol': 0.0001}, 
        {'padding': 4.0, 'iso_surf': 0.002, 'step_size': 0.088, 'tol': 3e-05}
        ]

axes_of_rotation = [
    np.array([0, 0, 0]),
    np.array([1, 0, 0]),
    np.array([0, 1, 0]),
    np.array([0, 0, 1]),
]
angles_in_degrees = [0, 45, 120]

epsiso_data_one_surf = {
    "epsiso_param_dict_list": epsiso_prop_dict_list_one_surf,
    "axes_of_rotation": axes_of_rotation,
    "angles_in_degrees": angles_in_degrees,
}

epsiso_data_tri_surf = {
    "epsiso_param_dict_list": epsiso_prop_dict_list_tri_surf,
    "axes_of_rotation": axes_of_rotation,
    "angles_in_degrees": angles_in_degrees,
}

# Default memory for running on Ziggy. All other defaults are set in slurmsubmitgenerator.
ZIGGY_MEMORY_LIMIT = 8  # this should be per CPU.

RESOURCE_PATH = pathlib.Path(__file__).parent / "resources"

# Default values for reactor directory and bulk submission script.
SUBMIT_FILENAME = "job_submit.sh"

REACTOR_NAME = "EPSISOCALCS"


def get_parameters_from_json_file(json_filename):
    """Get parameters from a JSON file.

    Parameters
    ----------
    json_filename : str
        JSON filename.

    Returns
    -------
    dict
        keyword arguments for calculation set up.

    """
    parameter_dict = None
    with open(json_filename, "r") as json_stream:
        parameter_dict = json.load(json_stream)
    return parameter_dict

# Default assignments for Cambridge HPC to Hunter account.

CAM_HPC_DEFAULTS = get_parameters_from_json_file((RESOURCE_PATH /"camhpcparameters.json").absolute().as_posix())

# Default Ziggy parameters.

ZIGGY_DEFAULTS = get_parameters_from_json_file((RESOURCE_PATH /"ziggyparameters.json").absolute().as_posix())


def get_cml_filename_list(directory):
    """Function for getting cml files. Generates regex to find cml files based on directory given.
    """
    cml_regex = ""
    if directory == "":
        cml_regex = "*.cml"
    elif directory[-1] == "/":
        cml_regex = directory + "*.cml"
    else:
        cml_regex = directory + "/*.cml"
    LOGGER.info("Regex: %s", cml_regex)
    return glob.glob(cml_regex)


def create_jobs_from_cml_files(cml_filename_list, **kwargs):
    """This generates jobs from the given filename.
    """
    return SubmitCreator.createJobsFromCML(
        kwargs.pop("submit_filename", SUBMIT_FILENAME),
        kwargs.pop("reactor_name", REACTOR_NAME),
        cml_filename_list,
        kwargs.pop("memory_limit", ZIGGY_MEMORY_LIMIT),
        epsiso_calc=True,
        epsiso_data=kwargs.pop("epsiso_data", epsiso_data_one_surf),
        **kwargs
    )


def create_jobs_from_directory(**kwargs):
    """Create jobs from given directory.
    """
    cml_filenames = get_cml_filename_list(directory=kwargs.pop("cml_directory", ""))
    scratch_dir_default = "./scratch"
    LOGGER.info("CML files found: %i", len(cml_filenames))
    LOGGER.info(cml_filenames)
    return create_jobs_from_cml_files(
        cml_filename_list=cml_filenames, scratch_dir=scratch_dir_default, **kwargs
    )

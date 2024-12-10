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
Script for parsing command line input.

Contains main function that is run when called from the command line.

The CLI parser has two main modes of operation:

:Authors:
    Mark Driver <mdd31>
"""

import logging
import argparse
import nwchemcmlutils.epsisocalccreator as epscreator

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.WARN)


def process_args(args):
    """Process CLI arguments.

    Parameters
    ----------
    args : argparse.ArgumentParser
        CLI arguments.

    Returns
    -------
    obj
        result of processing function.

    """
    return args.func(args)


def process_job_creator_args(args):
    """Process CLI arguments for NWChem job creation.

    Parameters
    ----------
    args : argparse.ArgumentParser
        CLI arguments.

    Returns
    -------
    fail_list : list
        List of molecules that failed to generate input files.
    sub_file : int
        Result of submission file creation.

    """
    CML_DIRECTORY = args.cml_directory
    reactor_name = args.reactor_name
    geom_opt = args.geom_opt
    fail_list = None
    sub_file = None
    job_parameters = {}
    if args.json_filename is not None:
        job_parameters = epscreator.get_parameters_from_json_file(args.json_filename)

    if args.trisurf:
        job_parameters["epsiso_data"] = epscreator.epsiso_data_tri_surf
    elif args.monosurf:
        job_parameters["epsiso_data"] = epscreator.epsiso_data_one_surf
    if args.hpc:
        job_parameters.update(**epscreator.CAM_HPC_DEFAULTS)
    elif args.ziggy:
        job_parameters.update(**epscreator.ZIGGY_DEFAULTS)
    
    fail_list, sub_file = epscreator.create_jobs_from_directory(
            cml_directory=CML_DIRECTORY,
            reactor_name=reactor_name,
            geom_opt=geom_opt,
            **job_parameters
    )
    
    return fail_list, sub_file



def create_argparser():
    """Parser for NWChem job creation.

    Returns
    -------
    argparse.ArgumentParser
        CLI option parser.

    """
    job_create_description = """nwchemcmlutils can generate the input scripts for the calculations, and also the functionality
used to perform operations during jobs.

When called as a module it uses a default CLI for setting up input scripts for
NWChem calculations for a couple of basic options. For non standard systems and
more control over job parameters, please import epsiocalccreator in a script.

The CLI allows you to select whether to generate inputs for University of
Cambridge HPC or the Hunter group Ziggy cluster. Options of mono or tri MEPS surface are also possible.
The CLI also features the ability to feed in a json file containing keyowrd:value pairs.
See documentation for more information aobut accepted parameters."""
    parser_job_create = argparse.ArgumentParser(
        description=job_create_description,
    )
    parser_job_create.add_argument(
        "cml_directory", type=str, help="Directory containing CML files."
    )
    parser_job_create.add_argument(
        "-r", "--reactor_name", type=str, help="reactor name"
    )
    parser_job_create.add_argument(
        "--geom_opt",
        dest="geom_opt",
        action="store_true",
        help="Run geometry optimisation to generate 3D structure CML. Default mode, is not required to be specified.",
    )
    parser_job_create.add_argument(
        "--no_geom_opt",
        dest="geom_opt",
        action="store_false",
        help="no geometry optimisation, just single point calculation to generate 3D structure CML.",
    )
    parser_job_create.add_argument("--json-file", "--json",
                                   dest="json_filename",
                                   type=str,
                              help="Method to feed in parameters from a JSON file to modify SLURM and other parameters.")
    ziggy_or_hpc = parser_job_create.add_mutually_exclusive_group()
    ziggy_or_hpc.add_argument(
        "-z", "--ziggy", action="store_true", help="Use SLURM defaults for Ziggy cluster."
    )
    ziggy_or_hpc.add_argument(
        "--hpc", "--camhpc", action="store_true", help="Use SLURM defaults for CSD3 HPC"
    )
    ziggy_or_hpc.add_argument("--other-hpc", action="store_true",
                              help="Use custom SLURM parameters from JSON file.")
    mono_or_tri = parser_job_create.add_mutually_exclusive_group()
    mono_or_tri.add_argument(
        "-m",
        "--monosurf",
        action="store_true",
        help="sets epsiso to calculate single surface",
    )
    mono_or_tri.add_argument(
        "-t",
        "--trisurf",
        action="store_true",
        help="sets epsiso to calculate multiple surfaces",
    )
    mono_or_tri.add_argument("-c","--customsurf", action="store_true",
                             help="Use EPSISO parameters stored in JSON file.")
    parser_job_create.set_defaults(func=process_job_creator_args, geom_opt=True)
    return parser_job_create


def main():
    """Main method: creates parser to process CLI input and execution.

    Returns
    -------
    fail_list, sub_file : list, int
        Returns list of molecules where input was not generated.

    """
    # Create parser
    parser = create_argparser()
    # parse args
    args = parser.parse_args()
    return process_args(args)

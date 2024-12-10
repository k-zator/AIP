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
Script to read in the cml and write the nwchem input files and ziggy submission
scripts.

@author: mark
"""

import logging
from lxml import etree
import os
import glob
import nwchemcmlutils.nwinFileGen.NWChemInputGenerator as NWInGen
import nwchemcmlutils.nwinFileGen.ZiggySubmitGenerator as ZiggyGen
import nwchemcmlutils.nwinFileGen.slurmsubmitgenerator as slurmgen
import nwchemcmlutils.calcCreator.CMLFileReader as CMLFileReader

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.WARN)


def writeNwinFile(
    molecule_cml,
    memory_limit,
    method_type="dft",
    tag_name="driverinitial",
    geom_name="geometry",
    directory="",
    **kwargs
):
    """writes nwin file in the given directory.
    """
    filename_stem = CMLFileReader.extractStdInChiKey(molecule_cml)
    molecule_cml_string = CMLFileReader.writeMoleculeCMLForNwinFile(molecule_cml)
    file_out = NWInGen.writeNwinFile(
        filename_stem=filename_stem,
        memory_limit=memory_limit,
        molecule_cml_string=molecule_cml_string,
        method_type=method_type,
        tag_name=tag_name,
        geom_name=geom_name,
        directory=directory,
        molecule_name=filename_stem,
        **kwargs
    )
    return file_out


def writeSLURMSubmitFile(molecule_cml, directory="", **kwargs):
    """Function writes ziggy submit file.
    """
    filename_stem = CMLFileReader.extractStdInChiKey(molecule_cml)
    # test to see if we are doing epsiso calculation- if true we name job
    if kwargs.get("epsiso_calc", False):
        job_name = "epsiso" + filename_stem
        kwargs["cube_clean"] = kwargs.get("cube_clean", True)
    else:
        job_name = "opt" + filename_stem
    file_out = slurmgen.writeSubmitFile(
        filename_stem=filename_stem,
        molecule_name=filename_stem,
        slurm_type=kwargs.pop("slurm_type", "ziggy"),
        directory=directory,
        job_name=job_name,
        **kwargs
    )
    return file_out


def writeZiggySubmitFile(molecule_cml, directory="", **kwargs):
    """Function writes ziggy submit file.
    """
    filename_stem = CMLFileReader.extractStdInChiKey(molecule_cml)
    # test to see if we are doing epsiso calculation- if true we name job
    if kwargs.get("epsiso_calc", False):
        job_name = "epsiso" + filename_stem
        kwargs["cube_clean"] = kwargs.get("cube_clean", True)
    else:
        job_name = "opt" + filename_stem
    file_out = ZiggyGen.writeSubmitFile(
        filename_stem=filename_stem,
        molecule_name=filename_stem,
        directory=directory,
        job_name=job_name,
        **kwargs
    )
    return file_out


def writeNwinAndSubmitFiles(
    molecule_cml,
    memory_limit,
    queue_system="SLURM",
    method_type="dft",
    tag_name="driverinitial",
    geom_name="geometry",
    directory="",
    **kwargs
):
    """Function writes the submit file and nwin file for a molecule.
    """
    nwin_file = writeNwinFile(
        molecule_cml,
        memory_limit,
        method_type,
        tag_name,
        geom_name,
        directory=directory,
        **kwargs
    )
    if queue_system == "SLURM":
        submit_file = writeSLURMSubmitFile(molecule_cml, directory=directory, **kwargs)
    elif queue_system == "TORQUE":
        submit_file = writeZiggySubmitFile(molecule_cml, directory=directory, **kwargs)
    else:
        raise ValueError("don't recognise this queueing system")
    return nwin_file, submit_file

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
script for creating NWChem input files.

@author: mark
"""

import logging
import nwchemcmlutils.nwinFileGen.NWChemCalcSetUp as NWChemCalcSetUp
import nwchemcmlutils.nwinFileGen.NWChemPythonSec as NWChemPythonSec

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.WARN)


def writeNwinFileContents(
    memory_limit,
    molecule_cml_string,
    method_type,
    tag_name,
    geom_name,
    filename_stem,
    epsiso_calc=False,
    epsiso_data=None,
    *import_args,
    **kwargs
):
    """function writes the text for an NWChem nwin file.
    """
    file_lines = ""
    file_lines += NWChemCalcSetUp.writeNwinSetUpLines(memory_limit, **kwargs)
    file_lines += NWChemCalcSetUp.writeDFTAndDriverLines(**kwargs)
    file_lines += NWChemPythonSec.setPythonBlock(
        molecule_cml_string,
        method_type,
        tag_name,
        geom_name,
        filename_stem,
        epsiso_calc,
        epsiso_data,
        *import_args,
        **kwargs
    )
    file_lines += NWChemCalcSetUp.setTask(task=kwargs.get("task", "python"))
    return file_lines


def writeNwinFileContentsToFile(filename, nwin_contents):
    """function writes the given contents to file.
    """
    with open(filename, "w") as nwin_file:
        nwin_file.write(nwin_contents)
    return 0


def writeNwinFile(
    filename_stem,
    memory_limit,
    molecule_cml_string,
    method_type,
    tag_name,
    geom_name,
    epsiso_calc=False,
    epsiso_data=None,
    directory="",
    *import_args,
    **kwargs
):
    """Function writes nwchem nwin file.
    """
    file_contents = writeNwinFileContents(
        memory_limit,
        molecule_cml_string,
        method_type,
        tag_name,
        geom_name,
        filename_stem,
        epsiso_calc,
        epsiso_data,
        *import_args,
        **kwargs
    )
    filename = directory + filename_stem + ".nwin"
    LOGGER.info("Writing to file")
    file_out = writeNwinFileContentsToFile(
        filename=filename, nwin_contents=file_contents
    )
    return file_out

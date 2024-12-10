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
Script for use during calculation for the reading in of molecule cml,
setting the geometry and basis sets, and then running a geometry optimisation,
and writing out the optimised structure to cml.

@author: mark
"""

import logging
import copy
import nwchemcmlutils.nwUtils.rtdb_util as rtdb_util
import nwchemcmlutils.nwUtils.cml_reading_util as cml_reading_util
import nwchemcmlutils.nwUtils.cml_writing_util as cml_writing_util
import nwchemcmlutils.nwinFileGen.NWChemCalcSetUp as NWChemCalcSetUp
from nwchemcmlutils.nwUtils.nwchemerror import NWChemError


logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.WARN)

try:
    import nwchem
    import nwgeom
except ImportError as err:
    LOGGER.warning(err)
    import unittest.mock as mock

    nwchem = mock.Mock()


def readCMLSetGeomAndBasis(molecule_cml_string, **basis_functions):
    """Function sets the geometry from cml, ready for calculation.
    """
    LOGGER.info("Reading in molecule cml")
    molecule_cml = cml_reading_util.createElementTreeFromString(molecule_cml_string)
    atom_array = cml_reading_util.extractAtomArray(molecule_cml)
    # now to read in the coordinates
    LOGGER.info("Setting coordinates.")
    atom_lines = cml_reading_util.createAtomLineListFromArray(atom_array)
    geometry = NWChemCalcSetUp.setGeometry(atom_lines)
    nwchem.input_parse(geometry)
    # now to set basis
    LOGGER.info("Setting basis functions.")
    mol_element_set = cml_reading_util.createElementSetFromArray(atom_array)
    basis = NWChemCalcSetUp.setBasis(mol_element_set, **basis_functions)
    nwchem.input_parse(basis)
    return molecule_cml


def geomOptCalc(method_type):
    """Function runs nwchem.task_optimise with the given method type.
    """
    nwchem.task_optimize(method_type)


def extractCoordWriteToCML(init_molecule_cml, tag_name, geom_name):
    """Function extracts the tag and coordinate information from the rtdb, and
    then creates a new cml molecule, and writes the new coordinates to file.
    """
    new_geom_dict = rtdb_util.geomGetAndMergeTagsAndCoords(tag_name, geom_name)
    new_molecule_cml = copy.deepcopy(init_molecule_cml)
    new_mol_atom_array = cml_reading_util.extractAtomArray(new_molecule_cml)
    if nwchem.ga_nodeid() == 0:
        cml_writing_util.setAtomCoordinatesForAtomArray(
            new_geom_dict, new_mol_atom_array
        )
    return new_molecule_cml


def writeCMLToFile(filename, molecule_cml):
    """Function writes cml to file.
    """
    if nwchem.ga_nodeid() == 0:
        return cml_writing_util.writeMoleculeCMLToFile(molecule_cml, filename)


def readCMLOptAndWriteOutputToFile(
    molecule_cml_string,
    method_type,
    tag_name,
    geom_name,
    filename,
    geom_opt=True,
    **basis_functions
):
    """Function acts as a wrapper for the other functions in this module- it
    reads in the coordinates and basis for the molecule, runs the optimisation
    and then outputs the cml.
    """
    init_molecule_cml = readCMLSetGeomAndBasis(molecule_cml_string, **basis_functions)
    if geom_opt:
        geomOptCalc(method_type)
        LOGGER.debug("rtdb contents:")
        LOGGER.debug(nwchem.rtdb_print(0))
        LOGGER.debug("tag name: %s", tag_name)
        new_molecule_cml = extractCoordWriteToCML(
            init_molecule_cml, tag_name, geom_name
        )
        out_file = writeCMLToFile(filename, new_molecule_cml)
    else:
        nwchem.task_energy(method_type)
        out_file = writeCMLToFile(filename, init_molecule_cml)
    return out_file

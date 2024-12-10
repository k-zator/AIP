# -*- coding: utf-8 -*-
#    cmlgenerator creates fully namespaced CML for molecules from input structures.
#    Copyright (C) 2019  Mark D. Driver
#
#    cmlgenerator is free software: you can redistribute it and/or modify
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
"""
Script for making CML from SMILES or from 3D structure and writing to file.

:Authors:
    Mark Driver <mdd31>
"""

import logging
from lxml import etree
import cmlgenerator.structuretocml as structcml
import cmlgenerator.smilestocml as smicml
import cmlgenerator.cmlnamespacing as cmlname
from filereader.cml_reader import CmlReader
from filereader.mol2reader import Mol2Reader

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.WARN)


def generate_cml_file_from_smiles(smiles_string, add_sybyl=False, add_aip_atom_types=False, **kwargs):
    """Generate CML file from SMILES string.

    Parameters
    ----------
    smiles_string : str
        SMILES string.
    max_confs : int, optional
        Maximum number of conformers to search for. The default is 10.
    MaxIters : int, optional
        Maximum number of iterations in calculation. The default is 1000.
    directory : str, optional
        directory to put output file, if different to CWD. The default is None.

    Returns
    -------
    None

    """
    directory = kwargs.pop("directory", None)
    inchikey, molecule_element = create_cml_from_smiles(smiles_string, add_sybyl, add_aip_atom_types, **kwargs)
    return cmlname.write_cml_to_file(molecule_element,
                                     create_filename_frominchikey(inchikey),
                                     directory=directory)


def generate_cml_file_from_file(input_file, fileformat, directory=None, add_sybyl=False, add_aip_atom_types=False, aromatic=False):
    """Generate CML file from other 3D input structure.

    Parameters
    ----------
    input_file : str
        input filename.
    file_format : str
        input file format.
    directory : str, optional
        directory to put output file, if different to CWD. The default is None.

    Returns
    -------
    None

    """
    inchikey, molecule_element = create_cml_from_file(input_file, fileformat, add_sybyl, add_aip_atom_types, aromatic)
    LOGGER.info("mol element")
    LOGGER.info(molecule_element)
    return cmlname.write_cml_to_file(molecule_element,
                                     create_filename_frominchikey(inchikey),
                                     directory)


def create_cml_from_smiles(smiles_string, add_sybyl = False, add_aip_atom_types = False, aromatic = False, **kwargs):
    """Create fully namespace qualified CML from SMILES.

    Parameters
    ----------
    smiles_string : str
        SMILES string.
    max_confs : int, optional
        Maximum number of conformers to search for. The default is 10.
    MaxIters : int, optional
        Maximum number of iterations in calculation. The default is 1000.

    Returns
    -------
    inchikey : str
        InChIKey
    cml : etree._Element
        cml:molecule with qualified attributes.

    """
    inchikey, molecule_cml = smicml.generate_cml_and_inchikey_from_smiles(
        smiles_string, **kwargs
    )


    if add_sybyl or add_aip_atom_types:
        if aromatic:
            molecule_cml2 = structcml.perform_conversion_from_string(molecule_cml, in_format="cml", aromatic=aromatic)
            cml = CmlReader(molecule_cml2, from_file=False, ns=None)
        else:
            cml = CmlReader(molecule_cml, from_file=False, ns=None)
        molecule_mol2 = structcml.perform_mol2_conversion_from_string(molecule_cml, in_format="cml")
        mol2 = Mol2Reader(molecule_mol2, from_file = False)
        if add_aip_atom_types:
            molecule_cml = etree.tostring(cml.get_cml_with_sybyl_and_aip_atom_types(mol2, aromatic), encoding="unicode")
        else:
            molecule_cml = etree.tostring(cml.get_cml_with_sybyl_info(mol2), encoding="unicode")
    return inchikey, cmlname.convert_molecule_from_string(inchikey, molecule_cml)


def create_cml_from_file(input_file, file_format, add_sybyl=False, add_aip_atom_types=False, aromatic=False):
    """Create fully namespace qualified CML 

    Parameters
    ----------
    input_file : str
        input filename.
    file_format : str
        input file format.

    Returns
    -------
    inchikey : str
        InChIKey
    cml : etree._Element
        cml:molecule with qualified attributes.

    """
    molecule_cml_orig = structcml.perform_conversion_from_file(input_file, file_format)
    molecule_cml = structcml.perform_conversion_from_string(molecule_cml_orig, "cml", aromatic)
    if add_sybyl:
        cml = CmlReader(molecule_cml, from_file=False, ns=None)
        molecule_mol2 = structcml.perform_mol2_conversion_from_string(molecule_cml_orig, in_format="cml")
        mol2 = Mol2Reader(molecule_mol2, from_file = False)
        if add_aip_atom_types:
            molecule_cml = etree.tostring(cml.get_cml_with_sybyl_and_aip_atom_types(mol2, aromatic), encoding="unicode")
        else:    
            molecule_cml = etree.tostring(cml.get_cml_with_sybyl_info(mol2), encoding="unicode")
    inchikey = structcml.generate_inchikey_for_cml(molecule_cml)
    LOGGER.info("inchikey: %s", inchikey)
    LOGGER.debug("molecule cml:")
    LOGGER.debug(molecule_cml)
    return inchikey, cmlname.convert_molecule_from_string(inchikey, molecule_cml)


def create_filename_frominchikey(inchikey):
    """Create CML filename from InChIKey.

    Parameters
    ----------
    inchikey : str
        InChIKey.

    Returns
    -------
    str
        output CML filename.

    """
    return inchikey + ".cml"

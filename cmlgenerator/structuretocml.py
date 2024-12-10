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
Script for generating CML from other 3D structure files accepted by openbabel.

See Also:
---------
openbabel : check documentation for accepted file formats.


:Authors:
    Mark Driver <mdd31>
"""
import logging
from openbabel import openbabel

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.WARN)


def generate_inchikey_for_cml(input_cml):
    """Generate InChIKey for cml in string.

    Parameters
    ----------
    input_cml : str
        CML of molecule.

    Returns
    -------
    str
        InChIKey.

    """
    conversion = create_inchikey_conversion("cml")
    ob_mol = openbabel.OBMol()
    if conversion.ReadString(ob_mol, input_cml):
        return conversion.WriteString(ob_mol).strip()


def create_inchikey_conversion(in_format):
    """Create converter for conversion to InChIKey.

    Parameters
    ----------
    in_format : str
        input format.

    Returns
    -------
    conversion : openbabel.OBConversion
        Conversion to inchikey.

    """
    conversion = openbabel.OBConversion()
    conversion.SetInFormat(in_format)
    conversion.SetOutFormat("inchi")
    conversion.SetOptions("K", conversion.OUTOPTIONS)
    return conversion


def perform_conversion_from_string(input_string, in_format, aromatic=False):
    """Convert string to CML with unqualified attributes.
    

    Parameters
    ----------
    input_string : str
        input representation.
    in_format : str
        string format.

    Returns
    -------
    str
        CML molecule representation.

    """
    conversion = create_cml_conversion(in_format, aromatic)
    #conversion.SetOptions("A", conversion.OUTOPTIONS)
    ob_mol = openbabel.OBMol()
    LOGGER.debug("Attempting conversion.")
    LOGGER.debug("InputString: %s", input_string)
    if conversion.ReadString(ob_mol, input_string):
        return conversion.WriteString(ob_mol)


def perform_conversion_from_file(infile_name, in_format, aromatic=False):
    """Convert input file to CML with unqualified attributes.

    Parameters
    ----------
    infile_name : str
        filename.
    in_format : str
        file format.

    Returns
    -------
    str
        CML molecule representation.

    """
    conversion = create_cml_conversion(in_format, aromatic)
    ob_mol = openbabel.OBMol()
    LOGGER.debug("Attempting conversion from %s.", in_format)
    LOGGER.debug("filename: %s", infile_name)
    if conversion.ReadFile(ob_mol, infile_name):
        return conversion.WriteString(ob_mol)


def create_cml_conversion(in_format, aromatic):
    """Create converter object which generates CML.

    Parameters
    ----------
    in_format : str
        input structure format.

    Returns
    -------
    conversion : openbabel.OBConversion
        conversion to CML from input format.

    """
    conversion = openbabel.OBConversion()
    conversion.SetInFormat(in_format)
    conversion.SetOutFormat("cml")
    if aromatic:
        conversion.SetOptions("A", conversion.OUTOPTIONS)
    # omit XML and namespace declarations
    conversion.SetOptions('x', openbabel.OBConversion.OUTOPTIONS)
    return conversion



def create_mol2_conversion(in_format):
    """Create converter object which generates CML.

    Parameters
    ----------
    in_format : str
        input structure format.

    Returns
    -------
    conversion : openbabel.OBConversion
        conversion to CML from input format.

    """
    conversion = openbabel.OBConversion()
    conversion.SetInFormat(in_format)
    conversion.SetOutFormat("mol2")
    return conversion

def perform_mol2_conversion_from_string(input_string, in_format):
    """Convert string to CML with unqualified attributes.
    

    Parameters
    ----------
    input_string : str
        input representation.
    in_format : str
        string format.

    Returns
    -------
    str
        CML molecule representation.

    """
    conversion = create_mol2_conversion(in_format)
    ob_mol = openbabel.OBMol()
    LOGGER.debug("Attempting conversion.")
    LOGGER.debug("InputString: %s", input_string)
    if conversion.ReadString(ob_mol, input_string):
        return conversion.WriteString(ob_mol)

def perform_mol2_conversion_from_file(infile_name, in_format):
    """Convert input file to CML with unqualified attributes.

    Parameters
    ----------
    infile_name : str
        filename.
    in_format : str
        file format.

    Returns
    -------
    str
        CML molecule representation.

    """
    conversion = create_mol2_conversion(in_format)
    ob_mol = openbabel.OBMol()
    LOGGER.debug("Attempting conversion from %s.", in_format)
    LOGGER.debug("filename: %s", infile_name)
    if conversion.ReadFile(ob_mol, infile_name):
        return conversion.WriteString(ob_mol)

def create_mol_conversion(in_format):
    """Create converter object which generates CML.

    Parameters
    ----------
    in_format : str
        input structure format.

    Returns
    -------
    conversion : openbabel.OBConversion
        conversion to CML from input format.

    """
    conversion = openbabel.OBConversion()
    conversion.SetInFormat(in_format)
    conversion.SetOutFormat("mol")
    return conversion

def perform_mol_conversion_from_string(input_string, in_format):
    """Convert string to CML with unqualified attributes.
    

    Parameters
    ----------
    input_string : str
        input representation.
    in_format : str
        string format.

    Returns
    -------
    str
        CML molecule representation.

    """
    conversion = create_mol_conversion(in_format)
    ob_mol = openbabel.OBMol()
    LOGGER.debug("Attempting conversion.")
    LOGGER.debug("InputString: %s", input_string)
    if conversion.ReadString(ob_mol, input_string):
        return conversion.WriteString(ob_mol)

def perform_mol_conversion_from_file(infile_name, in_format):
    """Convert input file to CML with unqualified attributes.

    Parameters
    ----------
    infile_name : str
        filename.
    in_format : str
        file format.

    Returns
    -------
    str
        CML molecule representation.

    """
    conversion = create_mol_conversion(in_format)
    ob_mol = openbabel.OBMol()
    LOGGER.debug("Attempting conversion from %s.", in_format)
    LOGGER.debug("filename: %s", infile_name)
    if conversion.ReadFile(ob_mol, infile_name):
        return conversion.WriteString(ob_mol)


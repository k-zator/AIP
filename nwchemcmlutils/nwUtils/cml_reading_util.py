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
Script for reading in CML block and creating 3D coordinate block ready for
input file.

Script also contains functions for converting outputted 3D coordinate block
information to CML

:Authors:
    Mark Driver <mdd31>
"""

import logging
import copy
from lxml import etree

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.WARN)

NAMESPACE_DICT = {"cml": "http://www.xml-cml.org/schema"}

CML_NAMESPACE_STRING = "{http://www.xml-cml.org/schema}"


def createElementTreeFromString(molecule_cml_string):
    """Create ElementTree from CML in a string.

    Parameters
    ----------
    molecule_cml_string : str
        Molecule CML.

    Returns
    -------
    molecule_cml : lxml.etree._ElementTree
        ElementTree representation of molecule XML.

    """
    molecule_cml = etree.XML(molecule_cml_string)
    return molecule_cml


def extractAtomArray(molecule_cml):
    """Extract atom array element from molecule CML.

    Parameters
    ----------
    molecule_cml : lxml.etree._ElementTree
        ElementTree representation of molecule XML.

    Raises
    ------
    ValueError
        Raised if multiple arrays are found.

    Returns
    -------
    atom_array : lxml.etree._ElementTree
        ElementTree representation of atom array XML.

    """
    atom_array = molecule_cml.findall("cml:atomArray", namespaces=NAMESPACE_DICT)
    if len(atom_array) == 1:
        atom_array = atom_array[0]
        LOGGER.debug("atomArray is %s", etree.tounicode(atom_array))
        return atom_array
    else:
        LOGGER.info("Number of atomArray elements is %i", len(atom_array))
        raise ValueError("Incorrect number of atom arrays.")


def extractBondArray(molecule_cml):
    """Extract bond array element from molecule CML.

    Parameters
    ----------
    molecule_cml : lxml.etree._ElementTree
        ElementTree representation of molecule XML.

    Raises
    ------
    ValueError
        Raised if multiple arrays are found.

    Returns
    -------
    bond_array : lxml.etree._ElementTree
        ElementTree representation of bond array XML.

    """
    bond_array = molecule_cml.findall("cml:bondArray", namespaces=NAMESPACE_DICT)
    if len(bond_array) == 1:
        bond_array = bond_array[0]
        LOGGER.debug("bondArray is %s", etree.tounicode(bond_array))
        return bond_array
    else:
        LOGGER.info("Number of bondArray elements is %i", len(bond_array))
        raise ValueError("Incorrect number of bond arrays.")


def createAtomList(atom_array_cml):
    """Get list of atoms present in atom array.

    Parameters
    ----------
    atom_array_cml : etree._Element
        CML atom array element.

    Returns
    -------
    atom_cml_list : list
        List of Atom CML elements.

    """
    atom_cml_list = atom_array_cml.findall("cml:atom", namespaces=NAMESPACE_DICT)
    LOGGER.debug("found %i atom elements", len(atom_cml_list))
    return atom_cml_list


def getAtomElementType(atom_cml):
    """Get Element of atom.

    Parameters
    ----------
    atom_cml : lxml.etree._ElementTree
        ElementTree representation of atom XML.

    Raises
    ------
    KeyError
        Raised if no attribute present.

    Returns
    -------
    atom_element_type : str
        Element name.

    """
    atom_element_type = atom_cml.get(CML_NAMESPACE_STRING + "elementType")
    if atom_element_type:
        return atom_element_type
    else:
        raise KeyError("no elementType")


def getAtomX3(atom_cml):
    """Get atom x3 coordinate.

    Parameters
    ----------
    atom_cml : lxml.etree._ElementTree
        ElementTree representation of atom XML.

    Raises
    ------
    KeyError
        Raised if no attribute present.

    Returns
    -------
    atom_x3 : str
        atom x3 coordinate.

    """
    atom_x3 = atom_cml.get(CML_NAMESPACE_STRING + "x3")
    if atom_x3:
        return atom_x3
    else:
        raise KeyError("no x3")


def getAtomY3(atom_cml):
    """Get atom y3 coordinate.

    Parameters
    ----------
    atom_cml : lxml.etree._ElementTree
        ElementTree representation of atom XML.

    Raises
    ------
    KeyError
        Raised if no attribute present.

    Returns
    -------
    atom_y3 : str
        atom y3 coordinate.

    """
    atom_y3 = atom_cml.get(CML_NAMESPACE_STRING + "y3")
    if atom_y3:
        return atom_y3
    else:
        raise KeyError("no y3")


def getAtomZ3(atom_cml):
    """Get atom z3 coordinate.

    Parameters
    ----------
    atom_cml : lxml.etree._ElementTree
        ElementTree representation of atom XML.

    Raises
    ------
    KeyError
        Raised if no attribute present.

    Returns
    -------
    atom_z3 : str
        atom z3 coordinate.

    """
    atom_z3 = atom_cml.get(CML_NAMESPACE_STRING + "z3")
    if atom_z3:
        return atom_z3
    else:
        raise KeyError("no z3")


def getAtomIsotopeNumber(atom_cml):
    """Get atom isotope number.

    Parameters
    ----------
    atom_cml : lxml.etree._ElementTree
        ElementTree representation of atom XML.

    Raises
    ------
    KeyError
        Raised if no attribute present.

    Returns
    -------
    atom_isotope_number : str
        Atom isotope number.

    """
    atom_isotope_number = atom_cml.get(CML_NAMESPACE_STRING + "isotopeNumber")
    if atom_isotope_number:
        return atom_isotope_number
    else:
        raise KeyError("no isotopeNumber")


def getAtomAttributes(atom_cml):
    """Get atom information- position and element information.

    Parameters
    ----------
    atom_cml : Tlxml.etree._ElementTree
        ElementTree representation of atom XML.

    Raises
    ------
    KeyError
        Raised if no attribute present for a field.

    Returns
    -------
    atom_dict : dict
        contains information about atom position and element type.

    """
    LOGGER.info("reading attributes for atom")
    try:
        atom_element_type = getAtomElementType(atom_cml)
        atom_x3 = getAtomX3(atom_cml)
        atom_y3 = getAtomY3(atom_cml)
        atom_z3 = getAtomZ3(atom_cml)
    except KeyError as err:
        raise err
    else:
        atom_dict = {
            "element": atom_element_type,
            "x3": atom_x3,
            "y3": atom_y3,
            "z3": atom_z3,
        }
    try:
        atom_isotope_number = getAtomIsotopeNumber(atom_cml)
    except KeyError:
        LOGGER.debug("no isotope number found")
        pass
    else:
        LOGGER.debug("Atom has isoptope number attribute")
        atom_dict["isotope_number"] = atom_isotope_number
    return atom_dict


def getAtomListAttributes(atom_cml_list):
    """Get atom information for all atoms in array.

    Parameters
    ----------
    atom_cml_list : list
        List of Atom CML elements.

    Returns
    -------
    atom_list_attribute_dict : dict
        Dictionary of all atom information- key is index in atom array.

    """
    LOGGER.info("Reading in Dictionary")
    atom_list_attribute_dict = {}
    i = 0
    for atom_cml in atom_cml_list:
        atom_attribute_dict = getAtomAttributes(atom_cml)
        atom_list_attribute_dict[i] = copy.deepcopy(atom_attribute_dict)
        i += 1
    return atom_list_attribute_dict


def convertAtomListAttribToAtomLineList(atom_list_attribute_dict):
    """Convert atom information to list of strings for input into NWChem.

    Parameters
    ----------
    atom_list_attribute_dict : dict
        Dictionary of all atom information- key is index in atom array.

    Returns
    -------
    atom_line_list : list
        List of formatted strings containing atom information for input into NWChem.

    """
    atom_line_list = []
    for number in atom_list_attribute_dict.keys():
        atom_dict = atom_list_attribute_dict[number]
        atom_line = "    %s" % (atom_dict["element"])
        atom_line += " " * 7
        if atom_dict["x3"][0] == "-":
            atom_line += atom_dict["x3"]
        else:
            atom_line += " " + atom_dict["x3"]
        atom_line += " " * 7
        if atom_dict["y3"][0] == "-":
            atom_line += atom_dict["y3"]
        else:
            atom_line += " " + atom_dict["y3"]
        atom_line += " " * 7
        if atom_dict["z3"][0] == "-":
            atom_line += atom_dict["z3"]
        else:
            atom_line += " " + atom_dict["z3"]
        if "isotope_number" in atom_dict.keys():
            atom_line += " mass" + " " * 8
            atom_line += atom_dict["isotope_number"]
        atom_line += "\n"
        atom_line_list.append(atom_line)
    return atom_line_list


def getElementSetFromAttributeDict(atom_list_attribute_dict):
    """Get set of unique atom element types.

    Parameters
    ----------
    atom_list_attribute_dict : dict
        Dictionary of all atom information- key is index in atom array.

    Returns
    -------
    element_set : set
        Elements in molecule.

    """
    element_set = set()
    for number in atom_list_attribute_dict.keys():
        element = atom_list_attribute_dict[number]["element"]
        element_set.add(element)
    return element_set


def createAtomLineListFromArray(atom_array_cml):
    """Create atom line entries for atom arrays.

    Parameters
    ----------
    atom_array_cml : etree._Element
        CML atom array element.

    Returns
    -------
    atom_line_list : list
        List of formatted strings containing atom information for input into NWChem.

    """
    atom_cml_list = createAtomList(atom_array_cml)
    atom_list_attribute_dict = getAtomListAttributes(atom_cml_list)
    atom_line_list = convertAtomListAttribToAtomLineList(atom_list_attribute_dict)
    return atom_line_list


def createElementSetFromArray(atom_array_cml):
    """Find unique elements in atom array CML.

    Parameters
    ----------
    atom_array_cml : etree._Element
        CML atom array element.

    Returns
    -------
    element_set : set
        Elements in molecule.

    """
    atom_cml_list = createAtomList(atom_array_cml)
    atom_list_attribute_dict = getAtomListAttributes(atom_cml_list)
    element_set = getElementSetFromAttributeDict(atom_list_attribute_dict)
    return element_set

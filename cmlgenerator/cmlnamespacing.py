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
# -*- coding: utf-8 -*-
"""
Script for converting standard CML to CML with namespace qualified attributes.

Also updates to CML format required for NWChem calculation.

:Authors:
    Mark Driver <mdd31>
"""

import logging
from lxml import etree
import pathlib
#import xmlvalidator.xmlvalidation as xmlvalidation

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.WARN)

CML_NAMESPACE_DICT = {
    "ssip": "http://www-hunter.ch.cam.ac.uk/SSIP",
    "cml": "http://www.xml-cml.org/schema",
}
CML_NAME = "{{{}}}".format(CML_NAMESPACE_DICT["cml"])

SSIP_NAME = "{{{}}}".format(CML_NAMESPACE_DICT["ssip"])


def convert_unnamespaced_file(filename, out_directory):
    ""
    #Assume file is named by inchikey- without using openbabel not able to generate InChIKey.
    inpath = pathlib.Path(filename)
    inchikey = inpath.name.replace(inpath.suffix,"")
    converted_molecule = convert_molecule_from_file(inchikey, filename)
    return write_cml_to_file(converted_molecule, inpath.name,
                             directory=out_directory)


def convert_molecule_from_file(inchikey, filename):
    """Generate molecule from CML file.

    Parameters
    ----------
    inchikey : str
        molecule inchikey.
    filename : str
        input CML molecule.

    Returns
    -------
    converted_molecule : etree._Element
        cml:molecule with qualified attributes.

    """
    molecule_element = etree.parse(filename)
    converted_molecule = convert_molecule(inchikey, molecule_element)
    LOGGER.debug("converted element:")
    LOGGER.debug(converted_molecule)
    return converted_molecule

def write_cml_to_file(molecule_element, filename, directory=None):
    """Write CML to file after asserting that it is valid against schema.

    Parameters
    ----------
    molecule_element : etree.Element
        CML molecule element.
    filename : str
        output filename.
    directory : str, optional
        directory. The default is None.

    Returns
    -------
    None

    """
    element_tree = etree.ElementTree(molecule_element)
    #xmlvalidation.validate_xml(element_tree, xmlvalidation.CML_SCHEMA)
    output_filename = (pathlib.Path(directory if directory != None else "") / filename).absolute().as_posix()
    return element_tree.write(
        output_filename, encoding="UTF-8", xml_declaration=True, pretty_print=True
    )

def convert_molecule_from_string(inchikey, cml_string):
    """Convert cml molecule to fully namespace qualified version with InChIKey.

    Parameters
    ----------
    inchikey : str
        inchikey.
    cml_string : inchikey
        CML string.

    Returns
    -------
    converted_molecule : etree._Element
        cml:molecule with qualified attributes.

    """
    molecule_element = read_molecule_string(cml_string)
    LOGGER.info("read in molecule string to namespace.")
    LOGGER.debug(etree.tounicode(molecule_element))
    converted_molecule = convert_molecule(inchikey, molecule_element)
    LOGGER.debug("converted element:")
    LOGGER.debug(converted_molecule)
    return converted_molecule


def convert_molecule(inchikey, molecule_element):
    """Namespace qualify molecule element and append InChIKey ot attributes.

    Parameters
    ----------
    inchikey : str
        InChIKey.
    molecule_element : etree._Element
        cml:molecule element with unqualified attributes.

    Returns
    -------
    converted_molecule : etree._Element
        cml:molecule with qualified attributes.

    """
    converted_molecule = namespace_molecule(molecule_element)
    append_inchikey(converted_molecule, inchikey)
    return converted_molecule


def append_inchikey(molecule_element, inchikey):
    """Append InChIKey as attribute to cml:molecule element.
    
    This sets both the cml:id and ssip:stdInChIKey attributes.

    Parameters
    ----------
    molecule_element : etree._Element
        cml:molecule element.
    inchikey : str
        InChIKey.

    Returns
    -------
    None.

    """
    molecule_element.set(CML_NAME + "id", inchikey)
    molecule_element.set(SSIP_NAME + "stdInChIKey", inchikey)


def namespace_molecule(molecule_element):
    """Namespace qualify cml:molecule element including children.

    Parameters
    ----------
    molecule_element : etree._Element
        cml:molecule element with unqualified attributes.

    Returns
    -------
    namespaced_molecule : TYPE
        DESCRIPTION.

    """
    namespaced_molecule = etree.Element(CML_NAME + "molecule", nsmap=CML_NAMESPACE_DICT)
    atom_array_element = molecule_element.xpath("atomArray")[0]
    namespaced_molecule.append(namespace_atomarray(atom_array_element))
    bond_array_elements = molecule_element.xpath("bondArray")
    if len(bond_array_elements) == 1:
        namespaced_molecule.append(namespace_bondarray(bond_array_elements[0]))
    return namespaced_molecule


def namespace_bondarray(bondarray_element):
    """Namespace qualify cml:bondArray element including children.

    Parameters
    ----------
    bondarray_element : etree._Element
        cml:bondArray element with unqualified attributes.

    Returns
    -------
    namespaced_bond_array : etree._Element
        Fully namespaced cml:bondarray element.

    """
    namespaced_bond_array = etree.Element(
        CML_NAME + "bondArray", nsmap=CML_NAMESPACE_DICT
    )
    bond_elements = bondarray_element.xpath("bond")
    for bond_element in bond_elements:
        namespaced_bond_array.append(namespace_bond(bond_element))
    return namespaced_bond_array


def namespace_bond(bond_element):
    """Namespace qualify cml:bond element.

    Parameters
    ----------
    bond_element : etree._Element
        cml:bond element with unqualified attributes.

    Returns
    -------
    namespaced_bond : etree._Element
        Fully namespaced cml:bond element.

    """
    namespaced_bond = etree.Element(CML_NAME + "bond", nsmap=CML_NAMESPACE_DICT)
    for attribute, value in bond_element.items():
        namespaced_bond.set(CML_NAME + attribute, value)
    return namespaced_bond


def namespace_atomarray(atomarray_element):
    """Namespace qualify cml:atomArray element including children.

    Parameters
    ----------
    atomarray_element : etree._Element
        cml:atomArray element with unqualified attributes.

    Returns
    -------
    namespaced_atom_array : etree._Element
        Fully namespaced cml:atomarray element.

    """
    namespaced_atom_array = etree.Element(
        CML_NAME + "atomArray", nsmap=CML_NAMESPACE_DICT
    )
    atom_elements = atomarray_element.xpath("atom")
    for atom_element in atom_elements:
        namespaced_atom_array.append(namespace_atom(atom_element))
    return namespaced_atom_array


def namespace_atom(atom_element):
    """Namespace qualify cml:atom element.

    Parameters
    ----------
    atom_element : etree._Element
        cml:atom element with unqualified attributes.

    Returns
    -------
    namespaced_atom : etree._Element
        Fully namespace qualified Atom element.

    """
    namespaced_atom = etree.Element(CML_NAME + "atom", nsmap=CML_NAMESPACE_DICT)
    for attribute, value in atom_element.items():
        namespaced_atom.set(CML_NAME + attribute, value)
    return namespaced_atom


def read_molecule_string(cml_string):
    """Read CML string to etree._Element.

    Parameters
    ----------
    cml_string : str
        Unnamespaced CML for molecule.

    Returns
    -------
    etree._Element
        Element object representation of CML.

    """
    return etree.XML(cml_string)

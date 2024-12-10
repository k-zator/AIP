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
Script to do all the outputting of cml from calculations.

:Authors:
    Mark Driver <mdd31>
"""

import logging
import copy
from lxml import etree

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.WARN)
CML_NAMESPACE_STRING = "{http://www.xml-cml.org/schema}"


def setAtomX3(x3_value, atom_cml):
    """Set atom x3 coordinate.

    Parameters
    ----------
    x3_value : str
        x3 value.
    atom_cml : etree._Element
        atom CML element.

    Returns
    -------
    None.

    """
    LOGGER.debug("Overwriting x3 value")
    atom_cml.attrib[CML_NAMESPACE_STRING + "x3"] = x3_value


def setAtomY3(y3_value, atom_cml):
    """Function takes a value and sets the y3 attribute of the atom_cml to
    this.
    """
    LOGGER.debug("Overwriting y3 value")
    atom_cml.attrib[CML_NAMESPACE_STRING + "y3"] = y3_value


def setAtomZ3(z3_value, atom_cml):
    """Function takes a value and sets the z3 attribute of the atom_cml to
    this.
    """
    LOGGER.debug("Overwriting z3 value")
    atom_cml.attrib[CML_NAMESPACE_STRING + "z3"] = z3_value


def setAtomCoordinates(coordinate_array, atom_cml):
    """Function takes a numpy array of coordinates, and then sets each
    coordinate in the atom cml.
    """
    LOGGER.debug("Setting all atom coordinataes")
    x3_value = str(coordinate_array[0])
    y3_value = str(coordinate_array[1])
    z3_value = str(coordinate_array[2])
    setAtomX3(x3_value, atom_cml)
    setAtomY3(y3_value, atom_cml)
    setAtomZ3(z3_value, atom_cml)


def setAtomCoordinatesForAtomArray(atom_dictionary, atom_array_cml):
    """Function takes an atom_dictionary and a cml atom array. It then sets the
    coordinates in the CML to match those in the atom dictionary.
    """
    LOGGER.info("Setting new coordinates")
    for atom_index in atom_dictionary.keys():
        LOGGER.debug("at index %i", atom_index)
        atom_entry = atom_dictionary[atom_index]
        LOGGER.debug("atom_entry is %s", atom_entry)
        coordinate_array = atom_entry["coordinates"]
        atom_cml = atom_array_cml[atom_index]
        setAtomCoordinates(coordinate_array, atom_cml)
    LOGGER.info("Finished setting coordinates")


def writeMoleculeCMLToFile(molecule_cml, filename):
    """function to write out the cml to a file.
    """
    molecule_cml_element_tree = etree.ElementTree(molecule_cml)
    LOGGER.debug("Writing molecule_cml_element_tree to file")
    molecule_cml_element_tree.write(
        filename, encoding="UTF-8", xml_declaration=True, pretty_print=True
    )
    return 0

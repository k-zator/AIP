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

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.WARN)


def validateXMLFile(file_etree, schema):
    """This checks to see if the file conforms to the given schema.
    """
    return schema.assertValid(file_etree)


def createCMLSchema(filename):
    """Reads in file to etree representation, then tries to create schema.
    """
    xml_file = readCMLFile(filename)
    return etree.XMLSchema(xml_file)


def readCMLFile(filename):
    """Function reads in cml file to ElementTree representation.
    """
    LOGGER.info("Reading in file: %s", filename)
    cml_file = etree.parse(filename)
    return cml_file


def extractMoleculeList(cml_file):
    """Function return a list of molecule cml elements.
    """
    LOGGER.info("extracting moecules")
    molecule_list = cml_file.xpath(
        "//cml:molecule", namespaces={"cml": "http://www.xml-cml.org/schema"}
    )
    return molecule_list


def extractStdInChiKey(cml_molecule):
    """returns the standard inchikey of the molecule.
    """
    LOGGER.info("Extracting inchikey")
    std_inchikey = cml_molecule.xpath(
        "@ssip:stdInChIKey", namespaces={"ssip": "http://www-hunter.ch.cam.ac.uk/SSIP"}
    )
    if len(std_inchikey) == 1:
        std_inchikey = std_inchikey[0]
        return std_inchikey
    else:
        LOGGER.warn("Not one inchikey in molecule.")
        LOGGER.info("Molecule:")
        LOGGER.debug(etree.tounicode(cml_molecule))
        return


def writeMoleculeCMLForNwinFile(cml_molecule):
    """function takes a cml molecule as an etree.Element and outputs a string
    representation. any line breaks: '\n' are replaced by '\n  ' so that
    when the python is parsed by nwchem no errors should be created- it should
    have the same indent as the rest of the file.
    """
    molecule_cml_string = etree.tounicode(cml_molecule)
    molecule_cml_string = molecule_cml_string.replace("\n", "\n  ")
    return molecule_cml_string

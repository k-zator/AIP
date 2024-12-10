import logging
import numpy as np
from mendeleev import element
from filereader.cml_reader import CmlReader, Atom

CML_NS = "http://www.xml-cml.org/schema"
logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.WARN)


class ATReader(CmlReader):
    def __init__(self, cml_file, scale=1):
        """Daughter class inheriting from CmlReader to also read sybyl and AIP atom types"""

        try:
            super().__init__(cml_file, from_file=True, scale=1, ns=CML_NS)
        except:
            LOGGER.critical("Error reading CML file, exiting")
            exit(1)
        self.sybyl = np.array([a.sybyl for a in self.list_atoms])
        self.aipAtomType = np.array([a.aipAtomType for a in self.list_atoms])

    def _get_atom_(self, cml_elem, scale):
        atom = Atom()
        name = cml_elem.attrib['{{{}}}id'.format(self.ns)]
        elem = cml_elem.attrib['{{{}}}elementType'.format(self.ns)]
        x, y, z = scale * float(cml_elem.attrib['{{{}}}x3'.format(self.ns)]),\
            scale * float(cml_elem.attrib['{{{}}}y3'.format(self.ns)]),\
            scale * float(cml_elem.attrib['{{{}}}z3'.format(self.ns)])
        try:
            sybyl = cml_elem.attrib['{{{}}}sybyl'.format(self.ns)]
        except:
            LOGGER.critical("No sybyl atom type in CML, exiting")
            exit(1)
        try:
            aipAtomType = cml_elem.attrib['{{{}}}aipAtomType'.format(self.ns)]
        except:
            LOGGER.critical("No aip atom type in CML, exiting")
            exit(1)

        atom.set_element(elem)
        atom.set_aname(name)
        atom.set_color(elem)
        atom.set_xyz(x, y, z)
        atom.set_sybyl(sybyl)
        atom.set_aipAtomType(aipAtomType)
        return atom

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
Script for the Atom class

@author: mdd31
"""

import logging


logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.WARN)


class Atom(object):
    """Class contains information about an atom.
    """

    def __init__(self, proton_number, cartesian_3d, id_num=0, mass=None):
        """initialise atom object.
        """
        self.id = id_num
        self.element_number = proton_number
        self.cartesian_3d = cartesian_3d
        self.mass = mass

    def __repr__(self):
        """Overload the __repr__ of the atom.
        """
        repr_string = "Atom({self.element_number},{self.cartesian_3d!r},id_num={self.id!r}, mass={self.mass})".format(
            self=self
        )
        return repr_string

    def __str__(self):
        """Overload the string representation, making use of the str of
        Cartesian3D
        """
        to_string = """Atom: element: {self.element_number}, mass {self.mass}
{self.cartesian_3d}
""".format(
            self=self
        )
        return to_string

    def __eq__(self, other_atom):
        """Overload __eq__ so that all attributes must be equal for objects to
        be equal.
        """
        return (
            self.element_number == other_atom.element_number
            and self.mass == other_atom.mass
            and self.cartesian_3d == self.cartesian_3d
        )

    def __hash__(self):
        """Overload so the hash representation of the object is a hash of the
        tuple of the atom's attributes.
        """
        md5_hash = self._md5_hash().hexdigest()
        return hash(md5_hash)

    def _md5_hash(self):
        """md5 hash representation of atom.
        """
        md5_hash = self.cartesian_3d._md5_hash()
        md5_hash.update(str(self.element_number).encode("utf-8"))
        # md5_hash.update(str(self.id))
        md5_hash.update(str(self.mass).encode("utf-8"))
        return md5_hash

    def translateAtom(self, cartesian_3d):
        """this translates the atom by the cartesian_3d coordinate.
        """
        try:
            LOGGER.debug("Trying to add Cartesian3D positions.")
            new_cartesian_3d = self.cartesian_3d.add(cartesian_3d)
        except ValueError as error:
            raise error
        except AttributeError as error:
            raise error
        else:
            new_atom = Atom(
                self.element_number, new_cartesian_3d, id_num=self.id, mass=self.mass
            )
            return new_atom

    def rotateAtom(self, rotation_matrix):
        """Rotate the atom, using the given rotation matrix.
        """
        LOGGER.debug("Rotating Atom.")
        new_cartesian_3d = self.cartesian_3d.rotate(rotation_matrix)
        new_atom = Atom(
            self.element_number, new_cartesian_3d, id_num=self.id, mass=self.mass
        )
        return new_atom

    def convertCoordinateUnitsToAU(self):
        """returns a new atom, with self.cartesian_3d in atomic units.
        """
        LOGGER.debug("Converting Cartesian3D to atomic units.")
        new_cartesian_3d = self.cartesian_3d.convertToAU()
        LOGGER.debug("Creating new Atom")
        new_atom = Atom(
            self.element_number, new_cartesian_3d, id_num=self.id, mass=self.mass
        )
        return new_atom

    def convertCoordinateUnitsToAng(self):
        """returns a new atom, with self.cartesian_3d in angstroms.
        """
        LOGGER.debug("Converting Cartesian3D to Angstroms")
        new_cartesian_3d = self.cartesian_3d.convertToAng()
        LOGGER.debug("Creating new Atom")
        new_atom = Atom(
            self.element_number, new_cartesian_3d, id_num=self.id, mass=self.mass
        )
        return new_atom

    def writeCubeFileLine(self):
        """Write out the atom as a cube file line. 
        """
        atom_fortran_format = "{:5d}{:12.6f}{:12.6f}{:12.6f}{:12.6f}"
        values = (
            self.element_number,
            0.0,
            self.cartesian_3d.vector_3d[0][0],
            self.cartesian_3d.vector_3d[1][0],
            self.cartesian_3d.vector_3d[2][0],
        )
        write_atom_line = atom_fortran_format.format(*values)
        return write_atom_line

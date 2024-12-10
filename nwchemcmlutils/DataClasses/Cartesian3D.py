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
Contains class for 3d cartesian coordinates, This inherits from ndarray.

@author: mark
"""

import logging
import numpy as np
import hashlib
from .Vector3D import Vector3D
from .DistanceUnits import DistanceUnits

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.WARN)


class Cartesian3D(object):
    """3D cartesian vector.
    """

    def __init__(self, coord_lst, unit=DistanceUnits.atomic_units):
        """initialise the vector
        """
        self.vector_3d = Vector3D(coord_lst)
        self.unit = unit

    def __repr__(self):
        """Overload the representation of the Cartesian3D
        """
        repr_string = (
            "Cartesian3D({self.vector_3d!r}, unit={se" + "lf.unit!r})"
        ).format(self=self)
        return repr_string

    def __str__(self):
        """The string representation of the object
        """
        to_string = """Coordinates are:
x={self.vector_3d[0][0]:.5f} y={self.vector_3d[1][0]:.5f} z={self.vector_3d[2][0]:.5f}
in units of: {self.unit}
""".format(
            self=self
        )
        return to_string

    def __eq__(self, other_cartesian_3d):
        """Redefine '==' operator. Cartesian3D.vector_3d and Cartesian3D.unit
        must both be equal for self == other_cartesian_3d.
        """
        return (
            self.vector_3d == other_cartesian_3d.vector_3d
            and self.unit == other_cartesian_3d.unit
        )

    def __hash__(self):
        """Overload the __hash__ method so that it returns the hash of the tuple:
        (rounded_vector_3d, self.unit) where rounded_vector_3d is 
        self.vector_3d.roundValuesVector3D()
        """
        md5_hash = self._md5_hash().hexdigest()
        return hash(md5_hash)

    def _md5_hash(self):
        """md5 hash of the object.
        """
        md5_hash = self.vector_3d._md5_hash()
        md5_hash.update(self.unit.name.encode("utf-8"))
        return md5_hash

    def length(self):
        """This returns the length of the 3d cartesian vector, and the unit.
        """
        length = self.vector_3d.length()
        return length, self.unit

    def add(self, cartesian_3d):
        """returns the sum of the cartesian vectors, only if they have the same
        units.
        """
        if self.unit == cartesian_3d.unit:
            LOGGER.debug("Units of both cartesian vectors match")
            new_vector_3d = self.vector_3d + cartesian_3d.vector_3d
            new_unit = self.unit
            new_cartesian_3d = Cartesian3D(new_vector_3d, new_unit)
            return new_cartesian_3d
        elif self.unit != cartesian_3d.unit:
            LOGGER.warning("Units of both cartesian vectors are different.")
            LOGGER.debug("self Cartesian3D is in: %s", self.unit)
            raise ValueError("Incorrect Units")

    def convertToAU(self):
        """returns the Cartesian3D of the Cartesian3D in atomic units.
        """
        conversion_factor = self.unit.conversionFactorToAU()
        new_vector_3d = np.multiply(self.vector_3d, conversion_factor)
        new_unit = DistanceUnits.atomic_units
        new_cartesian_3d = Cartesian3D(new_vector_3d, unit=new_unit)
        return new_cartesian_3d

    def convertToAng(self):
        """returns the Cartesian3D of the Cartesian3D in angstroms.
        """
        conversion_factor = self.unit.conversionFactorToAng()
        new_vector_3d = np.multiply(self.vector_3d, conversion_factor)
        new_unit = DistanceUnits.angstroms
        new_cartesian_3d = Cartesian3D(new_vector_3d, unit=new_unit)
        return new_cartesian_3d

    def rotate(self, rotation_matrix):
        """returns a new Cartesian3D after it has been rotated.
        """
        new_vector_3d = self.vector_3d.rotate(rotation_matrix)
        new_cartesian_3d = Cartesian3D(new_vector_3d, unit=self.unit)
        return new_cartesian_3d

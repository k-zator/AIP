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
Script contains class for storing a property value.

@author: mdd31
"""

import logging
import numpy as np

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.WARN)


class PropertyValue(object):
    """Class contains information about a property value.
    """

    def __init__(self, property_value, cartesian_3d):
        """This contains the property value at a point.
        """
        self.value = property_value
        self.cartesian_3d = cartesian_3d

    def __repr__(self):
        """Overload the __repr__ of the PropertyValue.
        """
        repr_string = "PropertyValue({self.value!r}, {self.cartesian_3d!r})".format(
            self=self
        )
        return repr_string

    def __str__(self):
        """Overloading of the str method.
        """
        to_string = """Value = {self.value}
Coordinate = {self.cartesian_3d}""".format(
            self=self
        )
        return to_string

    def __eq__(self, other_property_value):
        """Overload __eq__ so that they are only equal if all atributes are
        equal.
        """
        return (
            self.value == other_property_value.value
            and self.cartesian_3d == other_property_value.cartesian_3d
        )

    def __lt__(self, other_property_value):
        """Overloading less than so it compares the self.value attribute.
        """
        if self.cartesian_3d != other_property_value.cartesian_3d:
            if self.cartesian_3d.unit == other_property_value.cartesian_3d.unit:
                return self.value < other_property_value.value
            else:
                raise TypeError("Cartesian3D need same unit for comparison.")
        elif (
            self.cartesian_3d == other_property_value.cartesian_3d
            and self.value != other_property_value.value
        ):
            raise ValueError("Can't have 2 different Values in same position.")
        else:
            return False

    def __gt__(self, other_property_value):
        """Overloading greater than so it compares the self.value attribute.
        """
        if self.cartesian_3d != other_property_value.cartesian_3d:
            if self.cartesian_3d.unit == other_property_value.cartesian_3d.unit:
                return self.value > other_property_value.value
            else:
                raise TypeError("Cartesian3D need same unit for comparison.")
        elif (
            self.cartesian_3d == other_property_value.cartesian_3d
            and self.value != other_property_value.value
        ):
            raise ValueError("Can't have 2 different Values in same position.")
        else:
            return False

    def __hash__(self):
        """Overload to return hash of tuple of attributes.
        """
        md5_hash = self._md5_hash().hexdigest()
        return hash(md5_hash)

    def _md5_hash(self):
        """Creates the md5 hash of the PropertyValue.
        """
        md5_hash = self.cartesian_3d._md5_hash()
        md5_hash.update(str(self.value).encode("utf-8"))
        return md5_hash

    def translatePropertyValue(self, other_cartesian_3d):
        """Translate the PropertyValue to a new position, given by the Cartesian3D
        given.
        """
        try:
            LOGGER.debug("Trying to add Cartesian3D positions.")
            new_cartesian_3d = self.cartesian_3d.add(other_cartesian_3d)
        except ValueError as error:
            raise error
        except AttributeError as error:
            raise error
        else:
            new_property_value = PropertyValue(self.value, new_cartesian_3d)
            return new_property_value

    def rotatePropertyValue(self, rotation_matrix):
        """Rotate the PropertyValue, using the given rotation matrix.
        """
        LOGGER.debug("Rotating PropertyValue.")
        new_cartesian_3d = self.cartesian_3d.rotate(rotation_matrix)
        new_property_value = PropertyValue(self.value, new_cartesian_3d)
        return new_property_value

    def convertCoordinateUnitsToAU(self):
        """returns a new PropertyValue, with self.cartesian_3d in atomic units.
        """
        LOGGER.debug("Converting Cartesian3D to atomic units.")
        new_cartesian_3d = self.cartesian_3d.convertToAU()
        LOGGER.debug("Creating new PropertyValue")
        new_property_value = PropertyValue(self.value, new_cartesian_3d)
        return new_property_value

    def convertCoordinateUnitsToAng(self):
        """returns a new PropertyValue, with self.cartesian_3d in angstroms.
        """
        LOGGER.debug("Converting Cartesian3D to Angstroms")
        new_cartesian_3d = self.cartesian_3d.convertToAng()
        LOGGER.debug("Creating new PropertyValue")
        new_property_value = PropertyValue(self.value, new_cartesian_3d)
        return new_property_value

    def writeCubeFileLine(self):
        """Write out the PropertyValue as a cube file line.
        """
        value_fortran_format = "{:13.6f}{:13.6f}{:13.6f}{:13.6f}"
        values = (
            self.cartesian_3d.vector_3d[0][0],
            self.cartesian_3d.vector_3d[1][0],
            self.cartesian_3d.vector_3d[2][0],
            self.value,
        )
        write_prop_val_line = value_fortran_format.format(*values)
        return write_prop_val_line

    def returnValuesForPlotting(self):
        """Returns the values for plotting as a numpy array, in the order Cartesian3D
        x, y, z, value.
        """
        plot_values_array = np.array(
            [
                self.cartesian_3d.vector_3d[0][0],
                self.cartesian_3d.vector_3d[1][0],
                self.cartesian_3d.vector_3d[2][0],
                self.value,
            ]
        )
        return plot_values_array

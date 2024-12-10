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
Script for PropertyValueList class. this contains a list of PropertyValue
objects.

@author: mdd31
"""


import logging
import copy
import numpy as np
from .PropertyValue import PropertyValue
from .DistanceUnits import DistanceUnits
from .Cartesian3D import Cartesian3D

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.WARN)


class PropertyValueList(object):
    """Class contains a list of property values.
    """

    def __init__(self, property_value_list, unit=DistanceUnits.atomic_units):
        """Initialise PropertyValueList with a unit and list of PropertyValue
        objects.
        """
        self.unit = unit
        self.property_value_list = []
        self.addPropertyValueList(property_value_list)

    def __repr__(self):
        """Overload repr of PropertyValueList
        """
        repr_string = "PropertyValueList({self.property_value_list!r}, unit={self.unit!r})".format(
            self=self
        )
        return repr_string

    def __str__(self):
        """Overload PropertyValueList string method.
        """
        to_string = "PropertyValueList: {self.property_value_list}".format(self=self)
        return to_string

    def __eq__(self, other_property_value_list):
        """Overload equals.
        """
        return (
            self.property_value_list == other_property_value_list.property_value_list
            and self.unit == other_property_value_list.unit
        )

    def addPropertyValue(self, property_value):
        """Function to add an PropertyValue to the PropertyValueList.
        """
        if property_value.cartesian_3d.unit == self.unit:
            self.property_value_list.append(property_value)
        else:
            raise ValueError(
                "PropertyValue and PropertyValueList distance units do not match"
            )

    def addPropertyValueList(self, property_value_list):
        """Function adds a list of PropertyValues to the property_value_list
        attribute.
        """
        for property_value in property_value_list:
            self.addPropertyValue(property_value)

    def appendPropertyValueList(self, property_value_list):
        """function appends the PropertyValue entries from
        property_value_list.property_value_list to self.property_value_list
        """
        if self.unit == property_value_list.unit:
            merged_prop_val_list = copy.deepcopy(self)
            LOGGER.debug("Appending property vlaue list entries.")
            merged_prop_val_list.addPropertyValueList(
                property_value_list.property_value_list
            )
            return merged_prop_val_list
        else:
            raise ValueError("PropertyValueLists have different units.")

    def translatePropertyValueList(self, cartesian_3d):
        """Function translates molecule by given vector, translating all
        atoms, then returns a new molecule.
        """
        if cartesian_3d.unit == self.unit:
            translated_prop_val_list = []
            for property_value in self.property_value_list:
                translated_prop_val = property_value.translatePropertyValue(
                    cartesian_3d
                )
                translated_prop_val_list.append(translated_prop_val)
            translated_property_value_list = PropertyValueList(
                translated_prop_val_list, unit=self.unit
            )
            return translated_property_value_list
        elif cartesian_3d.unit != self.unit:
            raise ValueError("Cartesian3D has incorrect units for translation.")

    def rotatePropertyValueList(self, rotation_matrix):
        """Rotates the PropertyValueList by the given rotation_matrix,
        and returns the new position.
        """
        rotated_prop_val_list = []
        LOGGER.debug("Rotating PropertyValues and adding them to list")
        for prop_val in self.property_value_list:
            rotated_prop_val = prop_val.rotatePropertyValue(rotation_matrix)
            rotated_prop_val_list.append(rotated_prop_val)
        LOGGER.debug("Rotated PropertyValues, now creating new PropertyValueList.")
        rotated_property_value_list = PropertyValueList(
            rotated_prop_val_list, unit=self.unit
        )
        return rotated_property_value_list

    def rotatePropertyValueListAboutVector(self, rotation_matrix, translation_vector):
        """Rotates the PropertyValueList about a specified origin:

        PropertyValueList is first translated by the adding the negative
        translation vector. This is to position the grid at the new origin.

        It is then rotated by the rotation matrix.

        Finally it is translated back to it's original frame of reference (by
        adding the positive translation vector.)
        """
        LOGGER.info("Creating negative translation vector")
        neg_translation_vector = Cartesian3D(
            -translation_vector.vector_3d, unit=translation_vector.unit
        )
        LOGGER.info("Translating PropertyValueList to new origin")
        prop_val_list_new_origin = self.translatePropertyValueList(
            neg_translation_vector
        )
        LOGGER.info("rotating PropertyValueList about new origin")
        rot_prop_val_list = prop_val_list_new_origin.rotatePropertyValueList(
            rotation_matrix
        )
        LOGGER.info("Rotating PropertyValueList back to initial origin")
        rot_prop_val_list_old_origin = rot_prop_val_list.translatePropertyValueList(
            translation_vector
        )
        return rot_prop_val_list_old_origin

    def convertPropertyValueListToAU(self):
        """Convert to atomic units, and return new molecule.
        """
        LOGGER.debug("Converting PropertyValues to atomic units.")
        converted_prop_val_list = []
        for prop_val in self.property_value_list:
            prop_val_in_au = prop_val.convertCoordinateUnitsToAU()
            converted_prop_val_list.append(prop_val_in_au)
        LOGGER.debug("Creating new PropertyValueList")
        prop_val_list_in_au = PropertyValueList(
            converted_prop_val_list, unit=DistanceUnits.atomic_units
        )
        return prop_val_list_in_au

    def convertPropertyValueListToAng(self):
        """Convert to atomic units, and return new molecule.
        """
        LOGGER.debug("Converting PropertyValues to angstroms.")
        converted_prop_val_list = []
        for prop_val in self.property_value_list:
            prop_val_in_ang = prop_val.convertCoordinateUnitsToAng()
            converted_prop_val_list.append(prop_val_in_ang)
        LOGGER.debug("Creating new molecule")
        prop_val_list_in_ang = PropertyValueList(
            converted_prop_val_list, unit=DistanceUnits.angstroms
        )
        return prop_val_list_in_ang

    def returnMinValue(self):
        """Returns the PropertyValue with the lowest value attribute.
        """
        min_value = self.property_value_list[0]
        for prop_value in self.property_value_list[1:]:
            if prop_value < min_value:
                min_value = prop_value
        return min_value

    def returnMaxValue(self):
        """Returns the PropertyValue with the lowest value attribute.
        """
        max_value = self.property_value_list[0]
        for prop_value in self.property_value_list[1:]:
            if prop_value > max_value:
                max_value = prop_value
        return max_value

    def writecubeFileLines(self):
        """Function writes the lines for the cube file.
        """
        write_lines = []
        for prop_val in self.property_value_list:
            prop_line = prop_val.writeCubeFileLine()
            write_lines.append(prop_line)
        return write_lines

    def returnPlotValues(self):
        """Returns the values for plotting as a numpy array (4,n), in the order Cartesian3D
        x, y, z, value. 
        """
        plot_values = np.array([[]]).reshape(0, 4)
        for prop_val in self.property_value_list:
            prop_val_values = prop_val.returnValuesForPlotting()
            plot_values = np.append(plot_values, [prop_val_values], axis=0)
        plot_values_transposed = plot_values.T
        return plot_values_transposed

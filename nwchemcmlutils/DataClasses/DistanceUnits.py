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
Contains class which enumerates distance units, and also provides conversion
between them.

@author: mdd31
"""

import logging
import enum

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.WARN)

ATOMIC_UNITS_TO_ANGSTROMS = 0.529177249
ANGSTROMS_TO_ATOMIC_UNITS = 1.0 / ATOMIC_UNITS_TO_ANGSTROMS


class DistanceUnits(enum.Enum):
    """Class contains an enumeration of the different distance units
    """

    atomic_units = 1
    angstroms = 2

    def __repr__(self):
        """Overload the __repr__ of the enumeration, so it returns the string
        DistanceUnits.{self.name}
        """
        repr_string = "DistanceUnits.{self.name}".format(self=self)
        return repr_string

    def describe(self):
        """This returns the name and value of the enumeration.
        """
        return self.name, self.value

    def conversionFactorToAU(self):
        """This returns the conversion factor for the conversion between the
        inputted unit type and atomic units.
        """
        if self.value == 1:
            conversion_factor = 1.0
        elif self.value == 2:
            conversion_factor = ANGSTROMS_TO_ATOMIC_UNITS
        return conversion_factor

    def conversionFactorToAng(self):
        """This returns the conversion factor for the conversion between the
        inputted unti type and angstroms.
        """
        if self.value == 2:
            conversion_factor = 1.0
        elif self.value == 1:
            conversion_factor = ATOMIC_UNITS_TO_ANGSTROMS
        return conversion_factor

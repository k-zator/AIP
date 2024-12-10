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
This contains the Vector3D class. This shall be used to represent a 3D position
in space, without any units associated with it. It inherits from ndarray.

@author: mark
"""

import logging
import math
import decimal
import struct
import hashlib
import numpy as np

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.WARN)


class Vector3D(np.ndarray):
    """3D Vector array.
    """

    def __new__(cls, coord_lst):
        """overload the base so that the shape of the Vector3d is an ndarray
        of shape(3,1).
        """
        vector_3d = np.ndarray.__new__(cls, (3, 1), dtype="float64")
        if isinstance(coord_lst, list):
            if coord_lst:
                vector_3d.put((0, 1, 2), coord_lst)
            else:
                vector_3d.put((0, 1, 2), (0, 0, 0))
        elif isinstance(coord_lst, np.ndarray):
            vector_3d.put((0, 1, 2), coord_lst)
        return vector_3d

    def __hash__(self):
        """Define the hash to use the binary tuple of the numbers.
        """
        md5_hex = self._md5_hash().hexdigest()
        LOGGER.debug("Vector3D to hash, %s", self)
        LOGGER.debug("Value to hash: %s", md5_hex)
        return hash(md5_hex)

    def __mul__(self, other_vector_3d):
        """redefine multiplication of Vector3D, such that x*y==x.y for Vector3D.
        """
        return np.dot(self.T, other_vector_3d)

    def __eq__(self, other_vector_3d):
        """redefine equal to operator such that 'x==y' returns np.allclose(x,y)
        This uses absolute tolerance of 1e-05 (or 0.00001) for the comparison,
        as well as a relative tolerance of 1e-20 (this should be much smaller
        than the absolute tolerance unless dealing with units on the order of
        10e15/10e16 or greater- this is unlikely.)
        """
        if isinstance(other_vector_3d, Vector3D):
            return np.allclose(self, other_vector_3d, rtol=1e-20, atol=1e-03)
        else:
            return False

    def _binaryIthValue(self, index):
        """This returns the 32 bit binary representation of the float in the ith
        position as a string.
        solution from:

        http://stackoverflow.com/questions/16444726/binary-representation-of-float-in-python-bits-not-hex

        Print statements have been replaced by Logger messages.
        """
        ith_value = self[index][0]
        # this is to deal with the fact that -0. != 0. in the binary represntation.
        if ith_value == 0.0:
            ith_value = 0.0
        # Struct can provide us with the float packed into bytes. The '!' ensures that
        # it's in network byte order (big-endian) and the 'f' says that it should be
        # packed as a float. Alternatively, for double-precision, you could use 'd'.
        packed = struct.pack("!f", ith_value)
        LOGGER.debug("Packed: %s", repr(packed))

        # For each character in the returned string, we'll turn it into its corresponding
        # integer code point
        #
        # [62, 163, 215, 10] = [ord(c) for c in '>\xa3\xd7\n']
        # In python 3 conversion with ord is not required.
        # For each integer, we'll convert it to its binary representation.
        binaries = [bin(i) for i in packed]
        LOGGER.debug("Binaries: %s", binaries)

        # Now strip off the '0b' from each of these
        stripped_binaries = [s.replace("0b", "") for s in binaries]
        LOGGER.debug("Stripped: %s", stripped_binaries)

        # Pad each byte's binary representation's with 0's to make sure it has all 8 bits:
        #
        # ['00111110', '10100011', '11010111', '00001010']
        padded = [s.rjust(8, "0") for s in stripped_binaries]
        LOGGER.debug("Padded: %s", padded)
        # At this point, we have each of the bytes for the network byte ordered float
        # in an array as binary strings. Now we just concatenate them to get the total
        # representation of the float:
        return "".join(padded)

    def binaryTuple(self):
        """Returns a tuple of the binary representations of the vlaues in the
        Vector3D, with the index corresponding to that in the vector.
        """
        binary_values = []
        LOGGER.debug("Calculating binary values")
        for i in range(3):
            binary_value = self._binaryIthValue(i)
            binary_values.append(binary_value)
        return tuple(binary_values)

    def _md5_hash(self):
        """Function returns the hashlib md5 object. 
        """
        md5_hash = hashlib.md5()
        binary_tuple = self.roundValuesVector3D().binaryTuple()
        for binary_num in binary_tuple:
            md5_hash.update(binary_num.encode("utf-8"))
        return md5_hash

    def _roundIthValue(self, index, decimal_val, rounding):
        """round down the value at the ith position. returns the float.
        """
        LOGGER.debug(
            "Converting ith value to a Decimal. value is: %.7f", self[index][0]
        )
        ith_value = decimal.Decimal(str(self[index][0]))
        LOGGER.debug("Rounding decimal value.")
        rounded_ith_value = float(ith_value.quantize(decimal_val, rounding=rounding))
        LOGGER.debug("Rounded float is: %f", rounded_ith_value)
        return rounded_ith_value

    def roundValuesTuple(
        self, decimal_val=decimal.Decimal("0.00001"), rounding=decimal.ROUND_DOWN
    ):
        """Rounds all the values in the Vector3D. returns a tuple of the values.

        The convention chosen is to always round down values, and do this at
        the same tolerance as the __eq__ method, of 0.00001 or 1e-05 decimal
        places.
        """
        value_list = []
        LOGGER.debug("Rounding each value in Vector3D and appending to list.")
        for i in range(3):
            rounded_value = self._roundIthValue(i, decimal_val, rounding)
            value_list.append(rounded_value)
        LOGGER.debug("Returning Tuple of rounded values.")
        return tuple(value_list)

    def roundValuesVector3D(
        self, decimal_val=decimal.Decimal("0.01"), rounding=decimal.ROUND_DOWN
    ):
        """Returns the rounded Vector3D. Note the lower tolerance on the
        rounding.
        """
        rounded_values_tuple = self.roundValuesTuple(
            decimal_val=decimal_val, rounding=rounding
        )
        LOGGER.debug("Creating Vector3D from tuple.")
        rounded_vector_3d = Vector3D(list(rounded_values_tuple))
        return rounded_vector_3d

    def length(self):
        """This returns the length of the Vector3d.
        """
        length = math.sqrt(sum(self * self))
        return length

    def normalise(self):
        """Return the normalised vector.
        """
        return self / self.length()

    def cross(self, other_vector_3d):
        """return the cross product of self and other_vector_3d.
        """
        return Vector3D(np.cross(self.T, other_vector_3d.T))

    def rotate(self, rotation_matrix):
        """rotate the vec_3d and return the rotated Vector3D.
        """
        return Vector3D(np.dot(rotation_matrix, self))

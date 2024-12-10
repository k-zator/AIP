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
Script containing tests for the Atom Class.

@author: mdd31
"""

import unittest
import logging
import numpy as np
from nwchemcmlutils.DataClasses.Cartesian3D import Cartesian3D
from nwchemcmlutils.DataClasses.DistanceUnits import DistanceUnits
from nwchemcmlutils.DataClasses.Vector3D import Vector3D
from nwchemcmlutils.DataClasses.Atom import Atom

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.WARN)


class AtomTestCase(unittest.TestCase):
    """Test case for the Atom class.
    """

    def setUp(self):
        """Set up for tests.
        """
        cartesian_3d_ang = Cartesian3D([1, 1, 1], unit=DistanceUnits.angstroms)
        self.atom_o_ang = Atom(8, cartesian_3d_ang)
        cartesian_3d_au = Cartesian3D(
            [1.88972599, 1.88972599, 1.88972599], unit=DistanceUnits.atomic_units
        )
        self.atom_o_au = Atom(8, cartesian_3d_au)
        self.rot_matrix_45_100 = np.array(
            [[1, 0, 0], [0.0, 0.70710678, 0.70710678], [0.0, -0.70710678, 0.70710678]]
        )

    def tearDown(self):
        """Tear down for tests.
        """
        del self.atom_o_ang
        del self.atom_o_au
        del self.rot_matrix_45_100

    def test___repr__(self):
        """Test to see if correct repr is produced.
        """
        expected_repr = """Atom(8,Cartesian3D(Vector3D([[1.],
          [1.],
          [1.]]), unit=DistanceUnits.angstroms),id_num=0, mass=None)"""
        actual_repr = repr(self.atom_o_ang)
        self.assertMultiLineEqual(actual_repr, expected_repr)

    def test___str__(self):
        """Test to see if correct string is returned.
        """
        expected_string = """Atom: element: 8, mass None
Coordinates are:
x=1.88973 y=1.88973 z=1.88973
in units of: DistanceUnits.atomic_units

"""
        actual_string = str(self.atom_o_au)
        self.assertMultiLineEqual(actual_string, expected_string)

    def test___hash__(self):
        """Test to see that the expected hash is produced.
        """
        actual_hash = hash(self.atom_o_au)
        self.assertTrue(isinstance(actual_hash, int))

    def test__md5_hash(self):
        """Test to see if expected md5_hash is produced.
        """
        expected_md5 = "bad4dd978c3d213e79e10580e0496f4a"
        actual_md5 = self.atom_o_au._md5_hash().hexdigest()
        self.assertEqual(actual_md5, expected_md5)

    def test_translate_atom(self):
        """Test to see if molecule is correctly translated.
        """
        expected_atom = Atom(8, Cartesian3D([2, 2, 2], unit=DistanceUnits.angstroms))
        actual_atom = self.atom_o_ang.translateAtom(
            Cartesian3D([1, 1, 1], unit=DistanceUnits.angstroms)
        )
        self.assertEqual(actual_atom, expected_atom)
        # now we test to see if we get the expected errors
        with self.assertRaises(ValueError) as err_1:
            self.atom_o_ang.translateAtom(
                Cartesian3D([1, 1, 1], unit=DistanceUnits.atomic_units)
            )
        expected_args_1 = "Incorrect Units"
        actual_args_1 = err_1.exception.args[0]
        self.assertEqual(actual_args_1, expected_args_1)
        with self.assertRaises(AttributeError) as err_2:
            self.atom_o_ang.translateAtom(Vector3D([1, 1, 1]))
        expected_args_2 = "'Vector3D' object has no attribute 'unit'"
        actual_args_2 = err_2.exception.args[0]
        self.assertEqual(actual_args_2, expected_args_2)

    def test_rotate_atom(self):
        """Test to see if correct atom is returned after rotation.
        """
        expected_atom = Atom(
            8, Cartesian3D([1.0, 1.41421356, 0.0], unit=DistanceUnits.angstroms)
        )
        actual_atom = self.atom_o_ang.rotateAtom(self.rot_matrix_45_100)
        self.assertEqual(actual_atom, expected_atom)

    def test_convert_to_au_from_ang(self):
        """Test to see if conversion from an atom with a Cartesian3D position
        in angstroms is converted to one in atomic units.
        """
        expected_atom = self.atom_o_au
        actual_atom = self.atom_o_ang.convertCoordinateUnitsToAU()
        self.assertEqual(actual_atom, expected_atom)

    def test_convert_to_au_from_au(self):
        """Test to see if conversion from an atom with a Cartesian3D position
        in atomic units is converted to one in atomic units.
        """
        expected_atom = self.atom_o_au
        actual_atom = self.atom_o_au.convertCoordinateUnitsToAU()
        self.assertEqual(actual_atom, expected_atom)

    def test_convert_to_ang_from_au(self):
        """Test to see if conversion from an atom with a Cartesian3D position
        in atomic units is converted to one in angstroms.
        """
        expected_atom = self.atom_o_ang
        actual_atom = self.atom_o_au.convertCoordinateUnitsToAng()
        self.assertEqual(actual_atom, expected_atom)

    def test_convert_to_ang_from_ang(self):
        """Test to see if conversion from an atom with a Cartesian3D position
        in angstroms is converted to one in angstroms.
        """
        expected_atom = self.atom_o_ang
        actual_atom = self.atom_o_ang.convertCoordinateUnitsToAng()
        self.assertEqual(actual_atom, expected_atom)

    def test_write_atom_line(self):
        """Test to see if expected line is written out.
        """
        expected_line = "    8    0.000000    1.889726    1.889726    1.889726"
        actual_line = self.atom_o_au.writeCubeFileLine()
        self.assertEqual(actual_line, expected_line)

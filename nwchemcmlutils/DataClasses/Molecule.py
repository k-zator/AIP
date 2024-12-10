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
Created on Wed Apr 20 14:38:27 2016

@author: mark
"""

import logging
import numpy as np
from .DistanceUnits import DistanceUnits
from .Cartesian3D import Cartesian3D

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.WARN)


class Molecule(object):
    """Molcule class. Contains atom list and also unit.
    """

    def __init__(self, atom_list, unit=DistanceUnits.atomic_units):
        """Initialisation of Molecule. Set default unit to atomic units.
        """
        # want to make sure all atoms in the molecule have the same unit.
        self.unit = unit
        self.atom_list = set()
        self.addAtomList(atom_list)
        self.number_of_atoms = len(self.atom_list)

    def __repr__(self):
        """Overload __repr__ so that we can recreate the molecule.
        """
        repr_string = "Molecule({self.atom_list!r}, unit={self.unit!r})".format(
            self=self
        )
        return repr_string

    def __str__(self):
        """Overload the __str__ method so it lists the atoms in the molecule.
        """
        to_string = "Molecule: {self.atom_list}".format(self=self)
        return to_string

    def __eq__(self, other_molecule):
        """Overload the __eq__ method.
        """
        return (
            self.unit == other_molecule.unit
            and self.atom_list == other_molecule.atom_list
            and self.number_of_atoms == other_molecule.number_of_atoms
        )

    def addAtom(self, atom):
        """Function to add an atom to the molecule.
        """
        if atom.cartesian_3d.unit == self.unit:
            self.atom_list.add(atom)
        else:
            raise ValueError("Atom and Molecule distance units do not match")

    def addAtomList(self, atom_list):
        """Function adds a list of atoms to the atom_list attribute.
        """
        for atom in atom_list:
            self.addAtom(atom)

    def translateMolecule(self, cartesian_3d):
        """Function translates molecule by given vector, translating all
        atoms, then returns a new molecule.
        """
        if cartesian_3d.unit == self.unit:
            translated_atom_list = []
            for atom in self.atom_list:
                translated_atom = atom.translateAtom(cartesian_3d)
                translated_atom_list.append(translated_atom)
            translated_molecule = Molecule(translated_atom_list, unit=self.unit)
            return translated_molecule
        elif cartesian_3d.unit != self.unit:
            raise ValueError("Cartesian3D has incorrect units for translation.")

    def calculateCentroidVector(self):
        """Calculates the Centroid Vector, and returns the Cartesian3D for this.
        returns Cartesian3D in units=self.unit
        """
        coord_cartesian_vector_sum = Cartesian3D([0, 0, 0], unit=self.unit)
        for atom in self.atom_list:
            coord_cartesian_vector_sum = coord_cartesian_vector_sum.add(
                atom.cartesian_3d
            )
        LOGGER.debug("Total coordinate sum = %s", coord_cartesian_vector_sum)
        centroid_vector_3d = coord_cartesian_vector_sum.vector_3d / self.number_of_atoms
        centroid_cartesian_3d = Cartesian3D(centroid_vector_3d, unit=self.unit)
        return centroid_cartesian_3d

    def rotateMolecule(self, rotation_matrix):
        """Rotates the molecule by the given rotation_matrix, and returns the
        new position
        """
        rotated_atom_list = []
        LOGGER.debug("Rotating atoms and adding them to list")
        for atom in self.atom_list:
            rotated_atom = atom.rotateAtom(rotation_matrix)
            rotated_atom_list.append(rotated_atom)
        LOGGER.debug("Rotated Atoms, now creating new Molecule.")
        rotated_molecule = Molecule(rotated_atom_list, unit=self.unit)
        return rotated_molecule

    def rotateMoleculeAboutCentroid(self, rotation_matrix):
        """Calculates centroid vector. then translates centroid to origin.
        Rotates molecule, then translates molecule back to initial position.
        """
        LOGGER.info("Calculating Centroid vector.")
        centroid_cartesian_3d = self.calculateCentroidVector()
        centroid_cartesian_3d_neg = Cartesian3D(
            -centroid_cartesian_3d.vector_3d, unit=self.unit
        )
        LOGGER.info("Translating by centroid_cartesian_3d_neg")
        mol_at_origin = self.translateMolecule(centroid_cartesian_3d_neg)
        LOGGER.info("Rotating Molecule about origin.")
        rotated_mol_at_origin = mol_at_origin.rotateMolecule(rotation_matrix)
        LOGGER.info("Translating molecule back to initial centroid_cartesian_3d")
        rotated_mol_at_centroid = rotated_mol_at_origin.translateMolecule(
            centroid_cartesian_3d
        )
        return rotated_mol_at_centroid

    def convertMoleculeToAU(self):
        """Convert to atomic units, and return new molecule.
        """
        LOGGER.debug("Converting Atoms to atomic units.")
        converted_atom_list = []
        for atom in self.atom_list:
            atom_in_au = atom.convertCoordinateUnitsToAU()
            converted_atom_list.append(atom_in_au)
        LOGGER.debug("Creating new molecule")
        molecule_in_au = Molecule(converted_atom_list, unit=DistanceUnits.atomic_units)
        return molecule_in_au

    def convertMoleculeToAng(self):
        """Convert to atomic units, and return new molecule.
        """
        LOGGER.debug("Converting Atoms to atomic units.")
        converted_atom_list = []
        for atom in self.atom_list:
            atom_in_ang = atom.convertCoordinateUnitsToAng()
            converted_atom_list.append(atom_in_ang)
        LOGGER.debug("Creating new molecule")
        molecule_in_ang = Molecule(converted_atom_list, unit=DistanceUnits.angstroms)
        return molecule_in_ang

    def writeCubeFileLines(self):
        """Write the atom entries for a cube file. Returns a list, where the
        order is determined by the atom id number- this is to determine the
        order.
        """
        LOGGER.debug("Creating list with None values.")
        write_lines = [None] * self.number_of_atoms
        LOGGER.debug("Writing line for each atom")
        for atom in self.atom_list:
            atom_line = atom.writeCubeFileLine()
            if write_lines[atom.id] == None:
                write_lines[atom.id] = atom_line
            else:
                raise ValueError("Multiple atoms have the same id.")
        LOGGER.debug("returning completed list.")
        return write_lines

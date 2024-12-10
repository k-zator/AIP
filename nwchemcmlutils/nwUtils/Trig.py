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
Script for doing geometry operations on molecules in NWChem.


:Authors:
    Mark Driver <mdd31>
    Mark Williamson <mjw@mjw.name>
"""
import logging
import math
import numpy as np


logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.WARN)


def calculateCentroidVector(atoms):
    """Calculate the centroid position of the atom coordinates given.

    Parameters
    ----------
    atoms : np.array
        Array of atom coordinates, shape (natom,3).

    Returns
    -------
    centroid_vector : np.array
        Centroid vector position. Array shape (3,).

    """
    LOGGER.info("Setting sums of x, y,z coordinates and atom count == 0.")
    coord_vector_sum = np.zeros((3))

    atom_count = 0

    LOGGER.info("Going to iterate over all atoms.")
    for atom in atoms:
        LOGGER.debug("Adding atom coordinates to total sums. Atom: %s", atom)
        coord_vector_sum = np.add(coord_vector_sum, atom)

        atom_count += 1

    LOGGER.info("Calculating mean x, y, z coordinates.")

    centroid_vector = coord_vector_sum / atom_count
    LOGGER.debug("Calculated centroid_vector is: %s", centroid_vector)

    return centroid_vector


def generateRotationMatrix(axis, theta_in_radians):
    """Generate rotation matrix for axis and theta in radians.

    Uses the Euler-Rodrigues formula to generate the matrix. Based on question
    answer: `here http://stackoverflow.com/questions/6802577/python-rotation-of-3d-vector`

    Parameters
    ----------
    axis : np.array
        (3,) array for Axis.
    theta_in_radians : float
        Angle in radians.

    Raises
    ------
    ValueError
        raised if rotation about [0,0,0] is attempted.

    Returns
    -------
    np.array
        Rotation matrix.

    See Also
    --------
    Euler-Rodrigeus formula: `here http://en.wikipedia.org/wiki/Euler%E2%80%93Rodrigues_parameters`_
    """
    if not (axis - np.array([0, 0, 0])).any():
        LOGGER.debug("Can't rotate about [0, 0, 0]")
        raise ValueError("Can't rotate about [0, 0, 0]")
    else:
        axis = axis / math.sqrt(np.dot(axis, axis))

        a = math.cos(theta_in_radians / 2)
        b, c, d = -axis * math.sin(theta_in_radians / 2)

        return np.array(
            [
                [
                    a * a + b * b - c * c - d * d,
                    2 * (b * c - a * d),
                    2 * (b * d + a * c),
                ],
                [
                    2 * (b * c + a * d),
                    a * a + c * c - b * b - d * d,
                    2 * (c * d - a * b),
                ],
                [
                    2 * (b * d - a * c),
                    2 * (c * d + a * b),
                    a * a + d * d - b * b - c * c,
                ],
            ]
        )


def translateAtom(atom, translation_vector):
    """Translate atom by given vector.

    Parameters
    ----------
    atom : np.array
        Position of atom.
    translation_vector : np.array
        Translation vector.

    Returns
    -------
    np.array
        New Atom coordinate.

    """
    LOGGER.debug("Adding translation vector to atom coordinates.")
    translated_atom = np.add(atom, translation_vector)
    return translated_atom


def translateMolecule(atoms, translation_vector):
    """Translate all atoms in molecule by given translation vector.

    Parameters
    ----------
    atoms : np.array
        Array of atom coordinates, shape (natom,3).
    translation_vector : np.array
        Translation vector.

    Returns
    -------
    translated_atoms : np.array
        Array of translated atom coordinates.

    """
    LOGGER.debug("translateMolecule(): input translation_vector is: ")
    LOGGER.debug(translation_vector)
    LOGGER.debug("")
    LOGGER.debug("translateMolecule(): input atoms are: ")
    LOGGER.debug(atoms)
    LOGGER.debug("")

    translated_atoms = np.array([[]]).reshape(0, 3)

    for atom in atoms:
        LOGGER.debug("translateMolecule(): atom          " + str(atom))
        translated_atom = translateAtom(atom, translation_vector)
        LOGGER.debug("translateMolecule(): translated_atom" + str(translated_atom))
        translated_atoms = np.append(translated_atoms, [translated_atom], axis=0)

    LOGGER.debug("translateMolecule(): output atoms are:")
    LOGGER.debug(translated_atoms)
    LOGGER.debug("")

    return translated_atoms


def rotateAtomByMatrix(atom, rotation_matrix):
    """Rotate atom coordinates using matrix.

    Parameters
    ----------
    atom : np.array
        Position of atom.
    rotation_matrix : np.array
        Rotation matrix.

    Returns
    -------
    rotated_atom : np.array
        Position of rotated atom.

    """
    LOGGER.debug("Rotating atom.")
    LOGGER.debug("Initial coordinates: %s", atom)
    LOGGER.debug("Rotation matrix: %s", rotation_matrix)
    rotated_atom = np.dot(rotation_matrix, atom)
    LOGGER.debug("Rotated Atom is: %s", rotated_atom)
    return rotated_atom


def rotateMoleculeByMatrix(molecule_atoms, rotation_matrix):
    """Rotate molecule using matrix.

    Parameters
    ----------
    molecule_atoms : np.array
        Array of atom coordinates, shape (natom,3).
    rotation_matrix : np.array
        Rotation matrix shape (3,3).

    Returns
    -------
    rotated_molecule_atoms : np.array
        Array of rotated atom coordinates, shape (natom,3).

    """
    LOGGER.info("Rotating molecule: %s", molecule_atoms)
    LOGGER.debug("Rotation Matrix: %s", rotation_matrix)
    rotated_molecule_atoms = np.array([[]]).reshape(0, 3)
    for atom in molecule_atoms:
        rotated_atom = rotateAtomByMatrix(atom, rotation_matrix)
        rotated_molecule_atoms = np.append(
            rotated_molecule_atoms, [rotated_atom], axis=0
        )
    LOGGER.debug("Rotated molecule: %s", rotated_molecule_atoms)
    return rotated_molecule_atoms


def rotateMolecule(atoms, angle_in_degrees, axis_of_rotation):
    """Rotate molecule from input angle and axis.

    Parameters
    ----------
    atoms : np.array
        Array of atom coordinates, shape (natom,3).
    angle_in_degrees : float
        Angle of rotation in degrees.
    axis_of_rotation : np.array
        Axis of rotation.

    Raises
    ------
    ValueError
        Raised if rotation about (0,0,0) is attempted.

    Returns
    -------
    rotated_atoms : np.array
        Array of rotated atom coordinates, shape (natom,3).

    """
    LOGGER.debug("rotateMolecule(): input atoms are ")
    LOGGER.debug(atoms)
    LOGGER.debug("")
    try:
        rotation_matrix = generateRotationMatrix(
            axis_of_rotation, math.radians(angle_in_degrees)
        )
    except ValueError as err:
        raise err
    LOGGER.debug("rotateMolecule(): rotation_matrix is ")
    LOGGER.debug(rotation_matrix)
    LOGGER.debug("")

    rotated_atoms = rotateMoleculeByMatrix(atoms, rotation_matrix)

    LOGGER.debug("rotateMolecule(): rotated atoms are now:")
    LOGGER.debug(rotated_atoms)
    LOGGER.debug("")

    return rotated_atoms


def rotateMoleculeRelativeToTranslationVector(
    atoms, angle_in_degrees, axis_of_rotation, translation_vector
):
    """Translate molecule, rotate molecule then translate molecule back.
    
    Combination of operations is designed to lead to molecule being rotated
    about the centroid, if translation_vector is the distance of the origin
    to the molecular centroid.

    Parameters
    ----------
    atoms : np.array
        Array of atom coordinates, shape (natom,3).
    angle_in_degrees : float
        Angle of rotation in degrees.
    axis_of_rotation : np.array
        Axis of rotation.
    translation_vector : np.array
        Translation vector.

    Raises
    ------
    ValueError
        Raised if rotation about (0,0,0) is attempted.

    Returns
    -------
    translate_rotated_mol_back : np.array
        Array of transformed atom coordinates, shape (natom,3).

    """
    translate_mol_by_neg_vector = translateMolecule(
        atoms, np.negative(translation_vector)
    )
    try:
        rotate_mol_at_new_location = rotateMolecule(
            translate_mol_by_neg_vector, angle_in_degrees, axis_of_rotation
        )
    except ValueError as err:
        raise err
    translate_rotated_mol_back = translateMolecule(
        rotate_mol_at_new_location, translation_vector
    )

    return translate_rotated_mol_back


def rotateMoleculeAboutCentroid(atoms, angle_in_degrees, axis_of_rotation):
    """Rotate molecule about centroid.
    
    Molecule centroid is calculated and used as translation vector before
    rotation is computed.

    Parameters
    ----------
    atoms :  np.array
        Array of atom coordinates, shape (natom,3).
    angle_in_degrees : float
        Angle of rotation in degrees.
    axis_of_rotation : np.array
        Axis of rotation.

    Raises
    ------
    ValueError
        Raised if rotation about (0,0,0) is attempted.


    Returns
    -------
    rotated_molecule : np.array
        Array of transformed atom coordinates, shape (natom,3).

    """
    LOGGER.debug("Calculating centroid vector")
    centroid_vector = calculateCentroidVector(atoms)

    LOGGER.debug("Rotating molecule about given axis, relative to centroid.")
    try:
        rotated_molecule = rotateMoleculeRelativeToTranslationVector(
            atoms, angle_in_degrees, axis_of_rotation, centroid_vector
        )
    except ValueError as err:
        raise err
    return rotated_molecule

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
"""Methods for interaction with the RTDB of NWChem.

Contains functions for accessing and parsing the rtdb of an active
nwchem calculation.

:Authors:
    Mark Driver <mdd31>
"""

import logging
import copy
import numpy as np
import nwchemcmlutils.nwUtils.Trig as Trig
from nwchemcmlutils.nwUtils.nwchemerror import NWChemError

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.WARN)

try:
    import nwchem
    import nwgeom
except ImportError as err:
    LOGGER.warning(err)
    import unittest.mock as mock

    nwchem = mock.Mock()
    nwgeom = mock.Mock()




def geomGetCoordsInAngToNumpyArray(name):
    """Get atom coordinates as numpy array.

    Parameters
    ----------
    name : str
        name of property.

    Raises
    ------
    NWChemError
        Exception raised due to error accessing RTDB.

    Returns
    -------
    atoms : np.array
        Array of atom coordinates, shape (natom,3).

    """
    try:
        actualname = nwchem.rtdb_get(name)
    except NWChemError:
        actualname = name
    if actualname is None:
        actualname = name
    LOGGER.debug("Actual name: %s", actualname)
    
    coords = nwchem.rtdb_get("geometry:" + actualname + ":coords")
    units = nwchem.rtdb_get("geometry:" + actualname + ":user units")
    if units == "a.u.":
        factor = 1.0
    elif units == "angstroms":
        factor = nwchem.rtdb_get("geometry:" + actualname + ":angstrom_to_au")
    else:
        raise NWChemError("unknown units")

    # coords iterator
    i = 0
    coord = [0.0, 0.0, 0.0]

    # Create the numpy.array of the final atoms
    atoms = np.array([[]]).reshape(0, 3)

    while i < len(coords):
        atom_x = coords[i + 0] / factor
        atom_y = coords[i + 1] / factor
        atom_z = coords[i + 2] / factor
        # Debug
        # print str(x) + " " + str(y) + " " + str(z)

        coord = np.array([[atom_x, atom_y, atom_z]])
        atoms = np.append(atoms, coord, axis=0)

        i = i + 3

    LOGGER.debug(str("geomGetCoordsInAngToNpArray(): atoms"))
    LOGGER.debug(atoms)
    LOGGER.debug("")

    return atoms


def geomSetCoordsFromNumpyArray(name, numpy_array):
    """Set geometry in NWChem calculation to those read from input array.

    Parameters
    ----------
    name : str
        name of property.
    numpy_array : np.array
        Array of atom coordinates, shape (natom,3).

    Raises
    ------
    NWChemError
        Exception raised due to error accessing RTDB.

    Returns
    -------
    None.

    """
    LOGGER.debug("name: %s", name)
    try:
        actualname = nwchem.rtdb_get(name)
    except NWChemError:
        actualname = name
    if actualname is None:
        actualname = name
    LOGGER.debug("actualname: %s", actualname)
    units = nwchem.rtdb_get("geometry:" + actualname + ":user units")
    if units == "a.u.":
        factor = 1.0
    elif units == "angstroms":
        factor = nwchem.rtdb_get("geometry:" + actualname + ":angstrom_to_au")
    else:
        raise NWChemError("unknown units")

    coordinate_list = []

    for item in numpy_array:
        coordinate_list.append(item[0])
        coordinate_list.append(item[1])
        coordinate_list.append(item[2])

    nwgeom.geom_set_coords(name, coordinate_list)


def geomGetAtomTags(name):
    """Get the element names for each atom.

    Parameters
    ----------
    name : str
        name of property.

    Returns
    -------
    atom_tags : list
        Atom element names.

    """
    try:
        actual_name = nwchem.rtdb_get(name)
    except NWChemError:
        LOGGER.debug("Error with getting: %s",name)
        actual_name = name
    if actual_name is None:
        actual_name = name
    LOGGER.debug("name to use: %s", actual_name)
    atom_tags = nwchem.rtdb_get("geometry:" + actual_name + ":tags")
    return atom_tags


def mergeTagAndCoords(tags, coordinate_array):
    """Match coordinate and element type information, returning dictionary.

    Parameters
    ----------
    tags : list
        Atom element names.
    coordinate_array : np.array
        Array of atom coordinates, shape (natom,3).

    Returns
    -------
    atom_dictionary : dict
        Atom information key is index of atom, and values are the atom element name and position.

    """
    atom_dictionary = {}
    for i, element in enumerate(tags):
        ith_dict_entry = {"element": element, "coordinates": coordinate_array[i]}
        atom_dictionary[i] = copy.deepcopy(ith_dict_entry)
    return atom_dictionary


def geomGetAndMergeTagsAndCoords(tag_name, geom_name):
    """Get coordinates and atom element tags, assembles them into a dictionary.

    Parameters
    ----------
    tag_name : str
        name of property.
    geom_name : str
        name of geometry property.

    Returns
    -------
    atom_dictionary : dict
        Atom information key is index of atom, and values are the atom element name and position.

    """
    LOGGER.info("Getting atom tags")
    element_tags = geomGetAtomTags(tag_name)
    # To prevent a molecule made of a single two letter long atom ("Ca") from
    # being split, by the iterator in mergeTagAndCoords, into "C" and "a".
    if type(element_tags) == str:
        element_tags = [element_tags]
    LOGGER.info("Getting coordinates")
    geom_coords = geomGetCoordsInAngToNumpyArray(geom_name)
    LOGGER.info("Creating dictionary")
    atom_dictionary = mergeTagAndCoords(element_tags, geom_coords)
    LOGGER.debug("Dictionary:")
    LOGGER.debug(atom_dictionary)
    return atom_dictionary


def rotateMoleculeAboutCentroid(name, angle_in_degrees, axis_of_rotation):
    """Rotate the atoms about the centroid and set the new orientation.

    Parameters
    ----------
    name : str
        name of property.
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
    None.

    """
    LOGGER.info("Reading the atom coordinates.")
    atoms = geomGetCoordsInAngToNumpyArray(name)
    try:
        rotated_atoms = Trig.rotateMoleculeAboutCentroid(
            atoms, angle_in_degrees, axis_of_rotation
        )
    except ValueError as err:
        raise err
    LOGGER.info("Set the new atom coordinates.")
    LOGGER.debug("Rotated atoms:")
    LOGGER.debug(rotated_atoms)
    geomSetCoordsFromNumpyArray(name, rotated_atoms)


def get_number_of_voxels_in_surface():
    """Get number of voxels in surface.
    
    Number of grid points contained within the specificed isosurface.

    Returns
    -------
    n_voxels : int
        number of voxels.

    """
    try:
        n_voxels = nwchem.rtdb_get("prop:epsiso:nenc")
    except NWChemError:
        n_voxels = None
    return n_voxels


def setCubeFileName(cube_file_name):
    """Set cube filename for output

    Parameters
    ----------
    cube_file_name : str
        Output cube filename.

    Returns
    -------
    None.

    """
    LOGGER.info("Setting cube file name to %s", cube_file_name)
    nwchem.rtdb_put("prop:grid:output", cube_file_name)

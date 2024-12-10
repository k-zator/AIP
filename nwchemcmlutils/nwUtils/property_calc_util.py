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
"""Methods for calculation of MEPS properties.

Script containing the functions to carry out the the property calculation
during the python section of the job. Includes functionality to rotate the
molecule and recalculate the property, including changing the output file name.

:Authors:
    Mark Driver <mdd31>
"""

import logging
import copy
import numpy as np
import nwchemcmlutils.nwUtils.rtdb_util as rtdb_util
import nwchemcmlutils.nwinFileGen.NWChemCalcSetUp as NWChemCalcSetUp
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


def createCubeFilenameStem(cube_stem, iso_surf, angle_in_degrees, axis_of_rotation):
    """Create cube filename stem wihch includes molecule information.
    
    The cube stem is based (normally the InChIKey) on the molecule, with extra
    infomration added based on the molecule orientation and isosurface to calculate.
    This allows the production of unique file names.
    
    Example: water (XLYOFNOQVPJJNP-UHFFFAOYSA-N), on the 0.002 isosurface,
    in the initial rotation (angle=0, axis=0,0,0) the string should be
    'XLYOFNOQVPJJNP-UHFFFAOYSA-N_0.0020_0_000'

    Parameters
    ----------
    cube_stem : str
        stem for cube information- based on InChIKey.
    iso_surf : float
        isosurface value.
    angle_in_degrees : int
        angle of rotation in degrees.
    axis_of_rotation : array
        axis of rotation.

    Returns
    -------
    cube_filename_stem : str
        Cube filename stem with information about surface and molecule rotation.

    """
    cube_filename_stem = "{cube_stem}_{iso_surf:.4f}_{angle_in_degrees}_".format(
        cube_stem=cube_stem, iso_surf=iso_surf, angle_in_degrees=angle_in_degrees
    )
    cube_filename_stem += "{axis_of_rotation[0]}{axis_of_rotation[1]}{axis_of_rotation[2]}".format(
        axis_of_rotation=axis_of_rotation
    )
    return cube_filename_stem


def calculateProperty(cube_filename_stem, method_type):
    """Calculate Property, setting output filename first.
    
    A single point calculation is run beofre the property calculation to
    ensure electron density is computed.

    Parameters
    ----------
    cube_filename_stem : str
        Cube filename stem with information about surface and molecule rotation.
    method_type : str
        QM method type.

    Returns
    -------
    None.

    """
    cube_filename = cube_filename_stem + ".cube"
    rtdb_util.setCubeFileName(cube_filename)
    nwchem.task_energy(method_type)
    nwchem.task_property(method_type)


def setPropertyEPSISO(padding, step_size, iso_surf, tol):
    """Set Property section for EPS isosurface calculation.
    
    Parameters required for grid (padding and step size) and isosurface (value
    and tolerance).

    Parameters
    ----------
    padding : float
        padding around molecule.
    step_size : float
        Grid spacing step size.
    iso_surf : float
        electron isosurface density in atomic units.
    tol : float
        Allowed tolerance for a grid point to have correct density.

    Returns
    -------
    None.

    """
    property_block = NWChemCalcSetUp.setPropertyBlockEPSISO(
        padding, step_size, iso_surf, tol
    )
    nwchem.input_parse(property_block)


def setAndCalculatePropertyEPSISO(cube_filename_stem, method_type, epsiso_param_dict):
    """Set electron density properties before calcualtion of MEPS.
    
    Runs a single point energy claculation to ensure electron densities are
    already computed.

    Parameters
    ----------
    cube_filename_stem : str
        Cube filename stem with information about surface and molecule rotation.
    method_type : str
        QM method type.
    epsiso_param_dict : dict
        EPS electron density isosurface information.

    Returns
    -------
    int
        return number of voxels cointained by the surface.

    """
    padding = epsiso_param_dict["padding"]
    step_size = epsiso_param_dict["step_size"]
    iso_surf = epsiso_param_dict["iso_surf"]
    tol = epsiso_param_dict["tol"]
    setPropertyEPSISO(padding, step_size, iso_surf, tol)
    calculateProperty(cube_filename_stem, method_type)
    return rtdb_util.get_number_of_voxels_in_surface()


def setAndCalcPropEPSISOList(cube_stem, method_type, epsiso_param_dict_list, **kwargs):
    """Calculate MEPS for each given set of EPS iso paramters.

    Parameters
    ----------
    cube_stem : str
        stem for cube information- based on InChIKey.
    method_type : str
        QM method type.
    epsiso_param_dict_list :  list of dict
        EPS electron density isosurface information for all surfaces.
    angle_in_degrees : int
        angle of rotation in degrees.
    axis_of_rotation : array
        axis of rotation.

    Returns
    -------
    cube_filename_stem_dict : dict
        Dictionary of cube filename stem: number of voxels.

    """
    angle_in_degrees = kwargs.get("angle_in_degrees", 0)
    axis_of_rotation = kwargs.get("axis_of_rotation", np.array([0, 0, 0]))
    cube_filename_stem_dict = {}
    LOGGER.debug("Iterating over epsiso_param_dict_list")
    for epsiso_param_dict in epsiso_param_dict_list:
        cube_filename_stem = createCubeFilenameStem(
            cube_stem, epsiso_param_dict["iso_surf"], angle_in_degrees, axis_of_rotation
        )
        n_voxels = setAndCalculatePropertyEPSISO(
            cube_filename_stem, method_type, epsiso_param_dict
        )
        cube_filename_stem_dict[epsiso_param_dict["iso_surf"]] = {
            cube_filename_stem: n_voxels
        }
    return cube_filename_stem_dict


def rotateSetCalcPropEPSISOList(
    cube_stem, method_type, epsiso_param_dict_list, angle_in_degrees, axis_of_rotation
):
    """Calculate MEPS for each given set of EPS iso paramters in given orientation.

    Parameters
    ----------
    cube_stem : str
        stem for cube information- based on InChIKey.
    method_type : str
        QM method type.
    epsiso_param_dict_list :  list of dict
        EPS electron density isosurface information for all surfaces.
    angle_in_degrees : int
        angle of rotation in degrees.
    axis_of_rotation : array
        axis of rotation.

    Raises
    ------
    err
        Raises error if rotation is not possible.

    Returns
    -------
    cube_filename_stem_dict : dict
        Dictionary of cube filename stem: number of voxels.

    """
    try:
        rtdb_util.rotateMoleculeAboutCentroid(
            "geometry", angle_in_degrees, axis_of_rotation
        )
    except ValueError as err:
        if err.args[0] == "Can't rotate about [0, 0, 0]":
            LOGGER.info("Can't rotate about [0,0,0]")
            LOGGER.info("so calculating property in current orientation")
            cube_filename_stem_dict = setAndCalcPropEPSISOList(
                cube_stem, method_type, epsiso_param_dict_list
            )
        else:
            raise err
    else:
        if angle_in_degrees != 0:
            cube_filename_stem_dict = setAndCalcPropEPSISOList(
                cube_stem,
                method_type,
                epsiso_param_dict_list,
                angle_in_degrees=angle_in_degrees,
                axis_of_rotation=axis_of_rotation,
            )
            # we then want to rotate back to original position.
            rtdb_util.rotateMoleculeAboutCentroid(
                "geometry", -angle_in_degrees, axis_of_rotation
            )

        else:
            LOGGER.warning("rotating by 0 gives same molecule.")
            LOGGER.warning(" Calculating property in same orientation")
            cube_filename_stem_dict = setAndCalcPropEPSISOList(
                cube_stem, method_type, epsiso_param_dict_list
            )
    return cube_filename_stem_dict


def multiOrientationCalcPropEPSISOList(
    cube_stem, method_type, epsiso_param_dict_list, axes_of_rotation, angles_in_degrees
):
    """Carry out calculation of properties in all orientations for all entries.
    
    Iterates over angles as fast direction, and the axes as the slow direction:
    this is because you don't want to do multiple calculations about [0, 0, 0],
    as all produce the same result.

    Parameters
    ----------
    cube_stem :  str
        stem for cube information- based on InChIKey.
    method_type : str
        QM method type.
    epsiso_param_dict_list :  list of dict
        EPS electron density isosurface information for all surfaces.
    axes_of_rotation : list
        List of angles of rotation in degrees.
    angles_in_degrees : list
        List of axes of rotation.

    Returns
    -------
    cube_file_name_dict : dict of dicts
        Dictionary containing information organised by isosurface.
        Each isosurface entry has a dictionary containing cube filename: number
        of voxels. For combination and averaging in merged cube file.

    """
    cube_file_stem_dict = {}
    for axis_of_rotation in axes_of_rotation:
        if isinstance(axis_of_rotation, list):
            axis_of_rotation = np.array(axis_of_rotation)
        LOGGER.debug("Axis of rotation: %s", axis_of_rotation)
        if not (axis_of_rotation - np.array([0, 0, 0])).any():
            LOGGER.info("Doing calculation in inital position.")
            cube_file_dict = setAndCalcPropEPSISOList(
                cube_stem, method_type, epsiso_param_dict_list
            )
            LOGGER.debug("cube_file dict: %s", cube_file_dict)
            for iso_surf, cube_file_stem_list in cube_file_dict.items():
                if iso_surf in cube_file_stem_dict.keys():
                    cube_file_stem_dict[iso_surf].update(cube_file_stem_list)
                else:
                    cube_file_stem_dict[iso_surf] = cube_file_stem_list
            LOGGER.debug("cube_file_stem_dict: %s", cube_file_stem_dict)
        else:
            for angle_in_degrees in angles_in_degrees:
                LOGGER.debug("Angle: %f", angle_in_degrees)
                cube_file_dict = rotateSetCalcPropEPSISOList(
                    cube_stem,
                    method_type,
                    epsiso_param_dict_list,
                    angle_in_degrees,
                    axis_of_rotation,
                )
                LOGGER.debug("cube_file dict: %s", cube_file_dict)
                for iso_surf, cube_file_stem_list in cube_file_dict.items():
                    if iso_surf in cube_file_stem_dict.keys():
                        cube_file_stem_dict[iso_surf].update(cube_file_stem_list)
                    else:
                        cube_file_stem_dict[iso_surf] = cube_file_stem_list
                LOGGER.debug("cube_file_stem_dict: %s", cube_file_stem_dict)
    cube_file_name_dict = {}
    LOGGER.info("cube_stem_dict: %s", cube_file_stem_dict)
    for iso_surf, cube_file_stem_dict in cube_file_stem_dict.items():
        LOGGER.info("iso surf: %f", iso_surf)
        LOGGER.info("cube_file_stem_dict: %s", iso_surf)
        cube_file_iso_dict = {}
        for cube_file_stem, n_voxels in cube_file_stem_dict.items():
            cube_file_iso_dict[cube_file_stem + ".cube"] = n_voxels
        cube_file_name_dict[iso_surf] = cube_file_iso_dict
    return cube_file_name_dict

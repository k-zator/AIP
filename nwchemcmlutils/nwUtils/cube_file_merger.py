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
Script contains functions to rotate the produced cube files back to the initial
orientation, and then merge them together, to produce a single output file.

@author: mark
"""

import logging
import copy
import numpy as np
from nwchemcmlutils.DataClasses.GaussianCube import GaussianCube
from nwchemcmlutils.nwUtils.property_calc_util import createCubeFilenameStem
from nwchemcmlutils.nwUtils.filecleanup import clean_up_files

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.WARN)

INFO_LOGGER = logging.getLogger(__name__ + ".mergeinfo")
INFO_LOGGER.setLevel(logging.WARN)


def readCubeFile(filename):
    """Function reads in a cube file.
    """
    cube_file = GaussianCube()
    LOGGER.debug("Reading in file: %s", filename)
    cube_file.read(file_name=filename)
    return cube_file


def rotateCubeFile(cube_file, axis_of_rotation, angle_in_degrees):
    """Function rotates the given cube file and returns the result.
    """
    LOGGER.debug("Rotating cube file")
    rotated_cube_file = cube_file.rotateAboutAxis(angle_in_degrees, axis_of_rotation)
    return rotated_cube_file


def writeCubeFile(filename, cube_file):
    """This writes the cube file to file.
    """
    LOGGER.debug("Writing cube file.")
    cube_file.write(filename)


def rotateWriteCubeFile(cube_file, out_filename, axis_of_rotation, angle_in_degrees):
    """Function rotates it before writing to a new
    file.
    """
    rotated_cube_file = rotateCubeFile(cube_file, axis_of_rotation, angle_in_degrees)
    writeCubeFile(out_filename, rotated_cube_file)


def readCubeFileList(filename_list):
    """Function reads in a list of cube files, from different sources, and
    returns a dictionary. 
    """
    cube_file_dict = {}
    for filename in filename_list:
        cube_file_dict[filename] = readCubeFile(filename)
    return cube_file_dict


def rotateBackandWriteCubeFileDict(
    cube_file_dict, iso_surf, angles_in_degrees, axes_of_rotation
):
    """Function rotates and writes the cube files out.
    """
    rot_file_name_list = []
    for cube_filename, cube_file in cube_file_dict.items():
        for axis_of_rotation in axes_of_rotation:
            for angle_in_degrees in angles_in_degrees:
                orientation_description = "_{iso_surf:.4f}_{angle}_".format(
                    iso_surf=iso_surf, angle=angle_in_degrees
                )
                orientation_description += "{axis[0]}{axis[1]}{axis[2]}".format(
                    axis=axis_of_rotation
                )
                if (
                    orientation_description in cube_filename
                    and "_000" not in orientation_description
                ):
                    out_filename = cube_filename.replace(".cube", "rot.cube")
                    rotateWriteCubeFile(
                        cube_file, out_filename, axis_of_rotation, -angle_in_degrees
                    )
                    rot_file_name_list.append(out_filename)
                elif (
                    orientation_description in cube_filename
                    and "_000" in orientation_description
                ):
                    out_filename = cube_filename.replace(".cube", "rot.cube")
                    cube_file.write(out_filename)
                    rot_file_name_list.append(out_filename)
                else:
                    pass
    return rot_file_name_list


def createCubeFileNameDict(
    cube_stem, iso_surf_list, angles_in_degrees, axes_of_rotation
):
    """Function creates a list of all the cube filenames.
    """
    cube_filename_dict = {}
    for iso_surf in iso_surf_list:
        # group based on iso surface density
        iso_surf_file_list = []
        for axis_of_rotation in axes_of_rotation:
            if not (axis_of_rotation - np.array([0, 0, 0])).any():
                # molecule is not rotated in this case.
                cube_filename_stem = createCubeFilenameStem(
                    cube_stem, iso_surf, 0, axis_of_rotation
                )
                cube_filename = cube_filename_stem + ".cube"
                iso_surf_file_list.append(cube_filename)
            else:
                for angle_in_degrees in angles_in_degrees:
                    if angle_in_degrees != 0:
                        cube_filename_stem = createCubeFilenameStem(
                            cube_stem, iso_surf, angle_in_degrees, axis_of_rotation
                        )
                        cube_filename = cube_filename_stem + ".cube"
                        iso_surf_file_list.append(cube_filename)
                    else:
                        # this is also equivalent to no rotation.
                        cube_filename_stem = createCubeFilenameStem(
                            cube_stem, iso_surf, 0, np.array([0, 0, 0])
                        )
                        cube_filename = cube_filename_stem + ".cube"
                        iso_surf_file_list.append(cube_filename)
        iso_surf_file_set = set(iso_surf_file_list)
        cube_filename_dict[iso_surf] = iso_surf_file_set
    return cube_filename_dict


def mergeCubeFileList(cube_file_name_list):
    """Function reads in then merges cube files together.
    """
    cube_file_dict = readCubeFileList(cube_file_name_list)
    sorted_cube_file_name_list = sorted(cube_file_name_list)
    LOGGER.debug("Setting merged file equal to first cube file in list")
    voxel_volumes = {
        sorted_cube_file_name_list[0]
        .replace("rot.cube", ""): cube_file_dict[sorted_cube_file_name_list[0]]
        .calculateVoxelVolume()
    }
    merged_cube_file = cube_file_dict[sorted_cube_file_name_list[0]]
    for cube_file_name in sorted_cube_file_name_list[1:]:
        INFO_LOGGER.info("CUBE file name: %s", cube_file_name)
        try:
            merged_cube_file = merged_cube_file.addGausianCube(
                cube_file_dict[cube_file_name]
            )
        except ValueError as val_e:
            INFO_LOGGER.error("value error: %s", val_e)
            raise val_e
        voxel_volumes[cube_file_name.replace("rot.cube", "")] = cube_file_dict[
            sorted_cube_file_name_list[0]
        ].calculateVoxelVolume()
    return merged_cube_file, voxel_volumes


def rotateandMergeCubeFileDict(
    cube_file_dict, iso_surf, angles_in_degrees, axes_of_rotation
):
    """Function rotates the cube files back to original position, writes them
    to file, then reads them in and merges them.
    """
    rotated_cube_file_names = rotateBackandWriteCubeFileDict(
        cube_file_dict, iso_surf, angles_in_degrees, axes_of_rotation
    )
    merged_cube_file, voxel_volumes = mergeCubeFileList(rotated_cube_file_names)
    return merged_cube_file, voxel_volumes


def calculate_mean_volume(voxel_volumes, number_of_voxels):
    """This calculates the volume enclosed for each cube file, and returns the mean value
    """
    mean_mol_volume = 0.0
    for filename, n_voxels in number_of_voxels.items():
        voxel_volume = voxel_volumes[filename.replace(".cube", "")]
        molecule_volume = n_voxels * voxel_volume[0]
        mean_mol_volume += molecule_volume
    mean_mol_volume = mean_mol_volume / len(number_of_voxels)
    # note that the units should all be the same if a merge was successful on the cube files.
    return mean_mol_volume, voxel_volume[1]


def readRotateMergeWriteCubeFileList(
    cube_stem,
    filename_list,
    iso_surf,
    angles_in_degrees,
    axes_of_rotation,
    number_of_voxels,
):
    """Function reads in a list of cube_files, then rotates and merges them,
    before writing it out.
    """
    read_in_cube_file_dict = readCubeFileList(filename_list)
    try:
        merged_cube_file, voxel_volumes = rotateandMergeCubeFileDict(
            read_in_cube_file_dict, iso_surf, angles_in_degrees, axes_of_rotation
        )
    except ValueError as val_err:
        INFO_LOGGER.error("Value Error: %s", val_err)
        INFO_LOGGER.info("number of voxels: %s", number_of_voxels)
        raise val_err
    mean_mol_volume = calculate_mean_volume(voxel_volumes, number_of_voxels)
    out_filename = "{cube_stem}_{iso_surf:.4f}_merged.cube".format(
        cube_stem=cube_stem, iso_surf=iso_surf
    )
    header_line1 = (
        merged_cube_file.title_line_one
        + " Enclosed Volume: {:.6f} {:s}".format(
            mean_mol_volume[0], mean_mol_volume[1].describe()[0]
        )
    )
    merged_cube_file.title_line_one = header_line1
    merged_cube_file.write(out_filename)


def readRotateMergeWriteCubeFileDict(
    cube_stem,
    filename_dict,
    iso_surf_list,
    angles_in_degrees,
    axes_of_rotation,
    test=False,
):
    """Function reads in a list of cube_files for each isosurface,
    then rotates and merges them, before writing it out. This is done for each
    isosurface value.
    """
    for iso_surf in iso_surf_list:
        filename_list = filename_dict[iso_surf].keys()
        number_of_voxels = filename_dict[iso_surf]
        readRotateMergeWriteCubeFileList(
            cube_stem,
            filename_list,
            iso_surf,
            angles_in_degrees,
            axes_of_rotation,
            number_of_voxels,
        )
    LOGGER.debug("Finished loop.")
    LOGGER.debug("starting clean up")
    if not test:
        clean_up_files(cube_stem, len(iso_surf_list))

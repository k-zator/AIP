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
Script for creating the python section of an NWChem input file.

@author: mark
"""

import logging
import textwrap

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.WARN)


def setPythonImports(*import_args):
    """This sets up the python imports to be carried out in the calculation.
    The default_imports contains modules that are always imported.
    """
    default_imports = [
        "logging",
        "copy",
        "numpy as np",
        "nwchem",
        "nwgeom",
        "nwchemcmlutils.nwUtils.rtdb_util as rtdb_util",
        "nwchemcmlutils.nwUtils.geom_opt_cml as geom_opt_cml",
        "nwchemcmlutils.nwUtils.property_calc_util as property_calc_util",
        "nwchemcmlutils.nwUtils.cube_file_merger as cube_file_merger",
    ]
    import_lines = ""
    for default_import in default_imports:
        import_lines += "  import %s\n" % (default_import)
    # set up logging before the extra arguements.
    logging_initialisation = """
  logging.basicConfig()
  LOGGER = logging.getLogger(__name__)
  LOGGER.setLevel(logging.INFO)
"""
    import_lines += logging_initialisation
    for import_arg in import_args:
        # we add the extra imports to the default ones.
        # this is done by try except block
        try_block = "  try:\n"
        try_block += "      import %s\n" % (import_arg)
        try_block += "  except ImportError:\n"
        try_block += "      LOGGER.warn('import failed')\n"
        import_lines += try_block

    return import_lines


def setCMLStringinPython(molecule_cml_string):
    """Function takes a string representation of the cml. This is written as a
    variable in the python script in the calculation.
    """
    read_in_cml_block = ""
    wrapped_cml_string = textwrap.fill(
        molecule_cml_string, width=200, break_long_words=False, subsequent_indent="  "
    )
    mol_cml_line = "  molecule_cml_string = '''" + wrapped_cml_string + "'''\n"
    read_in_cml_block += mol_cml_line
    return read_in_cml_block


def setVariablesForOptCalc(
    method_type, tag_name, geom_name, filename, geom_opt=True, **basis_functions
):
    """Function takes the other arguments required for geometry optimisation,
    and writes section in job.
    """
    opt_variables = """  method_type = '{method_type}'
  tag_name = '{tag_name}'
  geom_name = '{geom_name}'
  filename = '{filename}'
""".format(
        method_type=method_type,
        tag_name=tag_name,
        geom_name=geom_name,
        filename=filename,
    )
    if geom_opt:
        opt_variables += "  geom_opt = True\n"
    else:
        opt_variables += "  geom_opt = False\n"
    if basis_functions:
        opt_variables += "  basis_functions = {basis_functions}\n".format(
            basis_functions=basis_functions
        )
    else:
        opt_variables += "  basis_functions = None\n"
    return opt_variables


def setGeomOptToCMLLinePython():
    """Function writes the line for the geometry optimisation assumes the
    variables already exist in the python section.
    """
    geom_opt_cml = """
  if basis_functions:
      geom_calc_out = geom_opt_cml.readCMLOptAndWriteOutputToFile(molecule_cml_string,
                                                                  method_type, tag_name,
                                                                  geom_name, filename,
                                                                  geom_opt,
                                                                  **basis_functions)
  else:
      geom_calc_out = geom_opt_cml.readCMLOptAndWriteOutputToFile(molecule_cml_string,
                                                                  method_type, tag_name,
                                                                  geom_name, filename,
                                                                  geom_opt)
"""
    return geom_opt_cml


def setGeomOptToCMLPython(
    molecule_cml_string, method_type, tag_name, geom_name, filename, **basis_functions
):
    """Function writes lines to carry out a geometry optimisation from given information.
    """
    geom_opt_sec = setCMLStringinPython(molecule_cml_string)
    geom_opt_sec += setVariablesForOptCalc(
        method_type, tag_name, geom_name, filename, **basis_functions
    )
    geom_opt_sec += setGeomOptToCMLLinePython()
    return geom_opt_sec


def setParametersEPSISOPython(
    cube_stem, method_type, epsiso_param_dict_list, axes_of_rotation, angles_in_degrees
):
    """Function writes lines, setting the variables needed for the EPSISO property
    calculation.
    """
    axes_of_rotation_str = "["
    for index, axis_of_rotation in enumerate(axes_of_rotation):
        if index != 0:
            axes_of_rotation_str += ", np.{axis_of_rotation!r}".format(
                axis_of_rotation=axis_of_rotation
            )
        elif index == 0:
            axes_of_rotation_str += "np.{axis_of_rotation!r}".format(
                axis_of_rotation=axis_of_rotation
            )
    axes_of_rotation_str += "]"
    parameter_lines = """
  cube_stem = '{cube_stem}'
  method_type = '{method_type}'
  epsiso_param_dict_list = {epsiso_param_dict_list}
  axes_of_rotation = {axes_of_rotation}
  angles_in_degrees = {angles_in_degrees}
""".format(
        cube_stem=cube_stem,
        method_type=method_type,
        epsiso_param_dict_list=epsiso_param_dict_list,
        axes_of_rotation=axes_of_rotation_str,
        angles_in_degrees=angles_in_degrees,
    )

    return parameter_lines


def setEPSISOCalcPython():
    """Function writes line for carrying out property calculation.
    """
    prop_calc_line = """
  cube_file_list = property_calc_util.multiOrientationCalcPropEPSISOList(cube_stem,
                                                                         method_type,
                                                                         epsiso_param_dict_list,
                                                                         axes_of_rotation,
                                                                         angles_in_degrees)
"""
    return prop_calc_line


def setParametersCubeRotation(epsiso_param_dict_list):
    """Function sets the extra parameters needed to rotate the cube files back.
    """
    epsiso_val_list = [
        epsiso_param_dict["iso_surf"] for epsiso_param_dict in epsiso_param_dict_list
    ]
    cube_rot_param_lines = """
  iso_surf_list = {epsiso_val_list}
""".format(
        epsiso_val_list=epsiso_val_list
    )
    return cube_rot_param_lines


def setCubeRotMergePython():
    """Function writes the lines required for the rotation and merging of the
    Cube files.
    """
    cube_merge_lines = """
  if nwchem.ga_nodeid() == 0:
      cube_file_merger.readRotateMergeWriteCubeFileDict(cube_stem, cube_file_list,
                                                        iso_surf_list,
                                                        angles_in_degrees,
                                                        axes_of_rotation)
"""
    return cube_merge_lines


def setEPSISOCalcMergeCubePython(
    cube_stem, method_type, epsiso_param_dict_list, axes_of_rotation, angles_in_degrees
):
    """Function sets all the variables for the EPSISO calculation, and runs it,
    then merges the Cube files based on iso surface density.
    """
    epsiso_calc_lines = setParametersEPSISOPython(
        cube_stem,
        method_type,
        epsiso_param_dict_list,
        axes_of_rotation,
        angles_in_degrees,
    )
    epsiso_calc_lines += setEPSISOCalcPython()
    epsiso_calc_lines += setParametersCubeRotation(epsiso_param_dict_list)
    epsiso_calc_lines += setCubeRotMergePython()
    return epsiso_calc_lines


def setPythonTaskOpt():
    """function to write the geometry optimisation in input.
    """
    optimise_line = "  nwchem.task_optimize('dft')\n"
    return optimise_line


def setPythonTaskProperty():
    """Function sets the line for carrying out the property calculation.
    """
    property_line = """  # Carry out the Property calculation
  #This is normally an EPSISO calculation
  nwchem.task_property('dft')
"""
    return property_line


def setPythonBlockStart():
    """function starts python block.
    """
    python_start = "python\n"
    return python_start


def setPythonEnd():
    """function sets the end statement of a block.
    """
    end_line = "end\n"
    return end_line


def setPythonBlock(
    molecule_cml_string,
    method_type,
    tag_name,
    geom_name,
    filename_stem,
    epsiso_calc=False,
    epsiso_data=None,
    *import_args,
    **kwargs
):
    """function writes the python block, calling the other methods to construct
    the block.
    """

    LOGGER.info("Writing python block")
    python_block_text = setPythonBlockStart()
    python_block_text += setPythonImports(*import_args)
    cml_filename = filename_stem + ".cml"
    basis_functions = kwargs.get("basis_functions")
    geom_opt = kwargs.get("geom_opt", True)
    if basis_functions:
        python_block_text += setGeomOptToCMLPython(
            molecule_cml_string,
            method_type,
            tag_name,
            geom_name,
            cml_filename,
            geom_opt=geom_opt,
            **basis_functions
        )
    else:
        python_block_text += setGeomOptToCMLPython(
            molecule_cml_string,
            method_type,
            tag_name,
            geom_name,
            cml_filename,
            geom_opt=geom_opt,
        )
    if epsiso_calc and epsiso_data:
        LOGGER.debug("Writing section for epsiso calculation.")
        python_block_text += setEPSISOCalcMergeCubePython(
            filename_stem,
            method_type,
            epsiso_data["epsiso_param_dict_list"],
            epsiso_data["axes_of_rotation"],
            epsiso_data["angles_in_degrees"],
        )
    elif epsiso_calc and not epsiso_data:
        raise TypeError("no epsiso_data present")
    else:
        pass
    python_block_text += setPythonEnd()
    return python_block_text

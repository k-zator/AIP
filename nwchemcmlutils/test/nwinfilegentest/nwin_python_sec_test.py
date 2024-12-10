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
Script contains the tests for the writing of the python section of the NWChem
input file.

@author: mark
"""


import unittest
import logging
import numpy as np
import pathlib
from lxml import etree
import nwchemcmlutils.nwUtils.cml_reading_util as cml_reading_util
import nwchemcmlutils.nwinFileGen.NWChemPythonSec as NWInPythonSec


logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.WARN)


class NWChemPythonSecTestCase(unittest.TestCase):
    """Test class for the NWChem input generator script.
    """

    def setUp(self):
        """setUp for the class. Test molecule is water to match the example
        calculation.
        """
        self.maxDiff = None
        self.memory_limit = 8
        self.molecule_name = "water_example"
        self.title = "Example input script for geometry optimisation"
        self.element_set = {"O", "H"}
        parent_directory = pathlib.Path(__file__).parents[1]
        water_cml_file = etree.parse(
            (parent_directory / "resources/examplecml.xml").absolute().as_posix()
        )
        self.cml_namespace = {"cml": "http://www.xml-cml.org/schema"}
        self.water_cml = water_cml_file.findall(
            "cml:molecule", namespaces=self.cml_namespace
        )[0]
        self.atom_array = cml_reading_util.extractAtomArray(self.water_cml)
        self.atom_line_list = cml_reading_util.createAtomLineListFromArray(
            self.atom_array
        )
        self.molecule_cml = (
            u'<cml:molecule xmlns:cml="http://www.xml-cml.o'
            + 'rg/schema" StdInChiKey="XLYOFNOQVPJJNP-UHFFFAO'
            + 'YSA-N" StuartId="13_01_20121_new" id="XLYOFNOQ'
            + 'VPJJNP-UHFFFAOYSA-N">\n  <cml:atomArray>\n  <cm'
            + 'l:atom elementType="O" id="a1" x3="-0.198000" '
            + 'y3="0.000000" z3="0.347000"/>\n  <cml:atom ele'
            + 'mentType="H" id="a2" x3="0.760000" y3="0.00000'
            + '0" z3="0.204000"/>\n  <cml:atom elementType="H'
            + '" id="a3" x3="-0.561000" y3="0.000000" z3="-0.'
            + '551000"/>\n  </cml:atomArray>\n  <cml:bondArray>'
            + '\n  <cml:bond atomRefs2="a1 a2" order="1"/>\n '
            + ' <cml:bond atomRefs2="a1 a3" order="1"/>\n  </c'
            + "ml:bondArray>\n  </cml:molecule>\n  "
        )
        self.filename_stem = "XLYOFNOQVPJJNP-UHFFFAOYSA-N"
        self.epsiso_prop_dict_list = [
            {"padding": 2.0, "step_size": 0.088, "iso_surf": 0.002, "tol": 0.00003},
            {"padding": 2.0, "step_size": 0.088, "iso_surf": 0.001, "tol": 0.000015},
        ]
        self.axes_of_rotation = [
            np.array([0, 0, 0]),
            np.array([1, 0, 0]),
            np.array([0, 1, 0]),
            np.array([0, 0, 1]),
        ]
        self.angles_in_degrees = [0, 45, 120]

    def tearDown(self):
        """tearDown for the class.
        """
        del self.memory_limit
        del self.molecule_name
        del self.title
        del self.element_set
        del self.cml_namespace
        del self.water_cml
        del self.atom_array
        del self.atom_line_list
        del self.molecule_cml
        del self.filename_stem
        del self.epsiso_prop_dict_list
        del self.axes_of_rotation
        del self.angles_in_degrees

    def test_set_python_imports(self):
        """test to see if correct python imports are written.
        """
        expected_line_1 = """  import logging
  import copy
  import numpy as np
  import nwchem
  import nwgeom
  import nwchemcmlutils.nwUtils.rtdb_util as rtdb_util
  import nwchemcmlutils.nwUtils.geom_opt_cml as geom_opt_cml
  import nwchemcmlutils.nwUtils.property_calc_util as property_calc_util
  import nwchemcmlutils.nwUtils.cube_file_merger as cube_file_merger

  logging.basicConfig()
  LOGGER = logging.getLogger(__name__)
  LOGGER.setLevel(logging.INFO)
"""
        actual_line_1 = NWInPythonSec.setPythonImports()
        self.assertMultiLineEqual(actual_line_1, expected_line_1)
        expected_line_2 = """  import logging
  import copy
  import numpy as np
  import nwchem
  import nwgeom
  import nwchemcmlutils.nwUtils.rtdb_util as rtdb_util
  import nwchemcmlutils.nwUtils.geom_opt_cml as geom_opt_cml
  import nwchemcmlutils.nwUtils.property_calc_util as property_calc_util
  import nwchemcmlutils.nwUtils.cube_file_merger as cube_file_merger

  logging.basicConfig()
  LOGGER = logging.getLogger(__name__)
  LOGGER.setLevel(logging.INFO)
  try:
      import sys
  except ImportError:
      LOGGER.warn('import failed')
"""
        actual_line_2 = NWInPythonSec.setPythonImports("sys")
        self.assertMultiLineEqual(actual_line_2, expected_line_2)

    def test_set_python_task_opt(self):
        """test to see if expected optimisation task line is written.
        """
        expected_line = "  nwchem.task_optimize('dft')\n"
        actual_line = NWInPythonSec.setPythonTaskOpt()
        self.assertEqual(actual_line, expected_line)

    def test_set_python_block_start(self):
        """test to see if expected start line for python block is written.
        """
        expected_line = "python\n"
        actual_line = NWInPythonSec.setPythonBlockStart()
        self.assertEqual(actual_line, expected_line)

    def test_set_python_end(self):
        """Test to see if expected end line is written.
        """
        expected_line = "end\n"
        actual_line = NWInPythonSec.setPythonEnd()
        self.assertEqual(actual_line, expected_line)

    def test_set_cml_string_in_python(self):
        """Test to see if correct block of python is written for the
        molecule cml in the calculation.
        """
        expected_line = u"""  molecule_cml_string = '''<cml:molecule xmlns:cml="http://www.xml-cml.org/schema" StdInChiKey="XLYOFNOQVPJJNP-UHFFFAOYSA-N" StuartId="13_01_20121_new" id="XLYOFNOQVPJJNP-UHFFFAOYSA-N">   <cml:atomArray>   <cml:atom
  elementType="O" id="a1" x3="-0.198000" y3="0.000000" z3="0.347000"/>   <cml:atom elementType="H" id="a2" x3="0.760000" y3="0.000000" z3="0.204000"/>   <cml:atom elementType="H" id="a3"
  x3="-0.561000" y3="0.000000" z3="-0.551000"/>   </cml:atomArray>   <cml:bondArray>   <cml:bond atomRefs2="a1 a2" order="1"/>   <cml:bond atomRefs2="a1 a3" order="1"/>   </cml:bondArray>
  </cml:molecule>'''
"""
        actual_line = NWInPythonSec.setCMLStringinPython(self.molecule_cml)
        self.assertMultiLineEqual(actual_line, expected_line)

    def test_set_vars_opt_calc(self):
        """Test to see if expected values are set.
        """
        expected_lines_1 = """  method_type = 'dft'
  tag_name = 'driverinitial'
  geom_name = 'geometry'
  filename = 'XLYOFNOQVPJJNP-UHFFFAOYSA-N.cml'
  geom_opt = True
  basis_functions = None
"""
        actual_lines_1 = NWInPythonSec.setVariablesForOptCalc(
            "dft", "driverinitial", "geometry", "XLYOFNOQVPJJNP-UHFFFAOYSA-N.cml"
        )
        self.assertMultiLineEqual(actual_lines_1, expected_lines_1)
        expected_lines_2 = """  method_type = 'dft'
  tag_name = 'driverinitial'
  geom_name = 'geometry'
  filename = 'XLYOFNOQVPJJNP-UHFFFAOYSA-N.cml'
  geom_opt = True
  basis_functions = {'O': '6-31G*'}
"""
        actual_lines_2 = NWInPythonSec.setVariablesForOptCalc(
            "dft",
            "driverinitial",
            "geometry",
            "XLYOFNOQVPJJNP-UHFFFAOYSA-N.cml",
            O="6-31G*",
        )
        self.assertMultiLineEqual(actual_lines_2, expected_lines_2)

    def test_set_geom_opt_cml_line(self):
        """Test to see if expected lines are written for carrying out geometry
        opt.
        """
        expected_lines = """
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
        actual_lines = NWInPythonSec.setGeomOptToCMLLinePython()
        self.assertMultiLineEqual(actual_lines, expected_lines)

    def test_set_geom_opt_cml(self):
        """test to see if complete section is written as expected for geometry
        optimisation.
        """
        expected_lines = u"""  molecule_cml_string = '''<cml:molecule xmlns:cml="http://www.xml-cml.org/schema" StdInChiKey="XLYOFNOQVPJJNP-UHFFFAOYSA-N" StuartId="13_01_20121_new" id="XLYOFNOQVPJJNP-UHFFFAOYSA-N">   <cml:atomArray>   <cml:atom
  elementType="O" id="a1" x3="-0.198000" y3="0.000000" z3="0.347000"/>   <cml:atom elementType="H" id="a2" x3="0.760000" y3="0.000000" z3="0.204000"/>   <cml:atom elementType="H" id="a3"
  x3="-0.561000" y3="0.000000" z3="-0.551000"/>   </cml:atomArray>   <cml:bondArray>   <cml:bond atomRefs2="a1 a2" order="1"/>   <cml:bond atomRefs2="a1 a3" order="1"/>   </cml:bondArray>
  </cml:molecule>'''
"""
        expected_lines += """  method_type = 'dft'
  tag_name = 'driverinitial'
  geom_name = 'geometry'
  filename = 'XLYOFNOQVPJJNP-UHFFFAOYSA-N.cml'
  geom_opt = True
  basis_functions = None

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
        actual_lines = NWInPythonSec.setGeomOptToCMLPython(
            self.molecule_cml,
            "dft",
            "driverinitial",
            "geometry",
            "XLYOFNOQVPJJNP-UHFFFAOYSA-N.cml",
        )
        self.assertMultiLineEqual(actual_lines, expected_lines)

    def test_params_epsiso_python(self):
        """Test to see if parameters are written as expected.
        """
        self.maxDiff = None
        expected_lines = """
  cube_stem = 'XLYOFNOQVPJJNP-UHFFFAOYSA-N'
  method_type = 'dft'
  epsiso_param_dict_list = [{'padding': 2.0, 'step_size': 0.088, 'iso_surf': 0.002, 'tol': 3e-05}, {'padding': 2.0, 'step_size': 0.088, 'iso_surf': 0.001, 'tol': 1.5e-05}]
  axes_of_rotation = [np.array([0, 0, 0]), np.array([1, 0, 0]), np.array([0, 1, 0]), np.array([0, 0, 1])]
  angles_in_degrees = [0, 45, 120]
"""
        actual_lines = NWInPythonSec.setParametersEPSISOPython(
            self.filename_stem,
            "dft",
            self.epsiso_prop_dict_list,
            self.axes_of_rotation,
            self.angles_in_degrees,
        )
        self.assertMultiLineEqual(actual_lines, expected_lines)

    def test_set_epsiso_prop_calc(self):
        """Test for the writing of the EPSISO calculation list.
        """
        expected_lines = """
  cube_file_list = property_calc_util.multiOrientationCalcPropEPSISOList(cube_stem,
                                                                         method_type,
                                                                         epsiso_param_dict_list,
                                                                         axes_of_rotation,
                                                                         angles_in_degrees)
"""
        actual_lines = NWInPythonSec.setEPSISOCalcPython()
        self.assertMultiLineEqual(actual_lines, expected_lines)

    def test_set_params_cube_file_rot(self):
        """Test to see if expected lines are produced.
        """
        expected_lines = """
  iso_surf_list = [0.002, 0.001]
"""
        actual_lines = NWInPythonSec.setParametersCubeRotation(
            self.epsiso_prop_dict_list
        )
        self.assertMultiLineEqual(actual_lines, expected_lines)

    def test_set_cube_rot_merge(self):
        """Test to see if expected lines are written.
        """
        expected_lines = """
  if nwchem.ga_nodeid() == 0:
      cube_file_merger.readRotateMergeWriteCubeFileDict(cube_stem, cube_file_list,
                                                        iso_surf_list,
                                                        angles_in_degrees,
                                                        axes_of_rotation)
"""
        actual_lines = NWInPythonSec.setCubeRotMergePython()
        self.assertMultiLineEqual(actual_lines, expected_lines)

    def test_set_epsiso_calc_cube_merge_python(self):
        """Test to see if expected section is written.
        """
        expected_lines = """
  cube_stem = 'XLYOFNOQVPJJNP-UHFFFAOYSA-N'
  method_type = 'dft'
  epsiso_param_dict_list = [{'padding': 2.0, 'step_size': 0.088, 'iso_surf': 0.002, 'tol': 3e-05}, {'padding': 2.0, 'step_size': 0.088, 'iso_surf': 0.001, 'tol': 1.5e-05}]
  axes_of_rotation = [np.array([0, 0, 0]), np.array([1, 0, 0]), np.array([0, 1, 0]), np.array([0, 0, 1])]
  angles_in_degrees = [0, 45, 120]

  cube_file_list = property_calc_util.multiOrientationCalcPropEPSISOList(cube_stem,
                                                                         method_type,
                                                                         epsiso_param_dict_list,
                                                                         axes_of_rotation,
                                                                         angles_in_degrees)

  iso_surf_list = [0.002, 0.001]

  if nwchem.ga_nodeid() == 0:
      cube_file_merger.readRotateMergeWriteCubeFileDict(cube_stem, cube_file_list,
                                                        iso_surf_list,
                                                        angles_in_degrees,
                                                        axes_of_rotation)
"""
        actual_lines = NWInPythonSec.setEPSISOCalcMergeCubePython(
            self.filename_stem,
            "dft",
            self.epsiso_prop_dict_list,
            self.axes_of_rotation,
            self.angles_in_degrees,
        )
        self.assertMultiLineEqual(actual_lines, expected_lines)

    def test_set_python_task_property(self):
        """Test to see if python task_property line is written as expected.
        """
        expected_task_prop = """  # Carry out the Property calculation
  #This is normally an EPSISO calculation
  nwchem.task_property('dft')
"""
        actual_task_prop = NWInPythonSec.setPythonTaskProperty()
        self.assertMultiLineEqual(actual_task_prop, expected_task_prop)

    def test_write_python_block(self):
        """test to see if the complete python block is as expected.
        """
        expected_block_1 = """python
  import logging
  import copy
  import numpy as np
  import nwchem
  import nwgeom
  import nwchemcmlutils.nwUtils.rtdb_util as rtdb_util
  import nwchemcmlutils.nwUtils.geom_opt_cml as geom_opt_cml
  import nwchemcmlutils.nwUtils.property_calc_util as property_calc_util
  import nwchemcmlutils.nwUtils.cube_file_merger as cube_file_merger

  logging.basicConfig()
  LOGGER = logging.getLogger(__name__)
  LOGGER.setLevel(logging.INFO)
"""
        expected_block_1 += u"""  molecule_cml_string = '''<cml:molecule xmlns:cml="http://www.xml-cml.org/schema" StdInChiKey="XLYOFNOQVPJJNP-UHFFFAOYSA-N" StuartId="13_01_20121_new" id="XLYOFNOQVPJJNP-UHFFFAOYSA-N">   <cml:atomArray>   <cml:atom
  elementType="O" id="a1" x3="-0.198000" y3="0.000000" z3="0.347000"/>   <cml:atom elementType="H" id="a2" x3="0.760000" y3="0.000000" z3="0.204000"/>   <cml:atom elementType="H" id="a3"
  x3="-0.561000" y3="0.000000" z3="-0.551000"/>   </cml:atomArray>   <cml:bondArray>   <cml:bond atomRefs2="a1 a2" order="1"/>   <cml:bond atomRefs2="a1 a3" order="1"/>   </cml:bondArray>
  </cml:molecule>'''
"""
        expected_block_1 += """  method_type = 'dft'
  tag_name = 'driverinitial'
  geom_name = 'geometry'
  filename = 'XLYOFNOQVPJJNP-UHFFFAOYSA-N.cml'
  geom_opt = True
  basis_functions = None

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
        expected_block_1 += "end\n"
        actual_block_1 = NWInPythonSec.setPythonBlock(
            self.molecule_cml,
            "dft",
            "driverinitial",
            "geometry",
            "XLYOFNOQVPJJNP-UHFFFAOYSA-N",
        )
        self.assertEqual(actual_block_1, expected_block_1)
        expected_block_2 = """python
  import logging
  import copy
  import numpy as np
  import nwchem
  import nwgeom
  import nwchemcmlutils.nwUtils.rtdb_util as rtdb_util
  import nwchemcmlutils.nwUtils.geom_opt_cml as geom_opt_cml
  import nwchemcmlutils.nwUtils.property_calc_util as property_calc_util
  import nwchemcmlutils.nwUtils.cube_file_merger as cube_file_merger

  logging.basicConfig()
  LOGGER = logging.getLogger(__name__)
  LOGGER.setLevel(logging.INFO)
"""
        expected_block_2 += u"""  molecule_cml_string = '''<cml:molecule xmlns:cml="http://www.xml-cml.org/schema" StdInChiKey="XLYOFNOQVPJJNP-UHFFFAOYSA-N" StuartId="13_01_20121_new" id="XLYOFNOQVPJJNP-UHFFFAOYSA-N">   <cml:atomArray>   <cml:atom
  elementType="O" id="a1" x3="-0.198000" y3="0.000000" z3="0.347000"/>   <cml:atom elementType="H" id="a2" x3="0.760000" y3="0.000000" z3="0.204000"/>   <cml:atom elementType="H" id="a3"
  x3="-0.561000" y3="0.000000" z3="-0.551000"/>   </cml:atomArray>   <cml:bondArray>   <cml:bond atomRefs2="a1 a2" order="1"/>   <cml:bond atomRefs2="a1 a3" order="1"/>   </cml:bondArray>
  </cml:molecule>'''
"""
        expected_block_2 += """  method_type = 'dft'
  tag_name = 'driverinitial'
  geom_name = 'geometry'
  filename = 'XLYOFNOQVPJJNP-UHFFFAOYSA-N.cml'
  geom_opt = True
  basis_functions = {'O': '6-311G**'}

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
        expected_block_2 += "end\n"
        actual_block_2 = NWInPythonSec.setPythonBlock(
            self.molecule_cml,
            "dft",
            "driverinitial",
            "geometry",
            "XLYOFNOQVPJJNP-UHFFFAOYSA-N",
            basis_functions={"O": "6-311G**"},
        )
        self.assertEqual(actual_block_2, expected_block_2)

    def test_write_python_block_epsiso(self):
        """Test to see if expected python block is written 
        """
        expected_block = """python
  import logging
  import copy
  import numpy as np
  import nwchem
  import nwgeom
  import nwchemcmlutils.nwUtils.rtdb_util as rtdb_util
  import nwchemcmlutils.nwUtils.geom_opt_cml as geom_opt_cml
  import nwchemcmlutils.nwUtils.property_calc_util as property_calc_util
  import nwchemcmlutils.nwUtils.cube_file_merger as cube_file_merger

  logging.basicConfig()
  LOGGER = logging.getLogger(__name__)
  LOGGER.setLevel(logging.INFO)
"""
        expected_block += u"""  molecule_cml_string = '''<cml:molecule xmlns:cml="http://www.xml-cml.org/schema" StdInChiKey="XLYOFNOQVPJJNP-UHFFFAOYSA-N" StuartId="13_01_20121_new" id="XLYOFNOQVPJJNP-UHFFFAOYSA-N">   <cml:atomArray>   <cml:atom
  elementType="O" id="a1" x3="-0.198000" y3="0.000000" z3="0.347000"/>   <cml:atom elementType="H" id="a2" x3="0.760000" y3="0.000000" z3="0.204000"/>   <cml:atom elementType="H" id="a3"
  x3="-0.561000" y3="0.000000" z3="-0.551000"/>   </cml:atomArray>   <cml:bondArray>   <cml:bond atomRefs2="a1 a2" order="1"/>   <cml:bond atomRefs2="a1 a3" order="1"/>   </cml:bondArray>
  </cml:molecule>'''
"""
        expected_block += """  method_type = 'dft'
  tag_name = 'driverinitial'
  geom_name = 'geometry'
  filename = 'XLYOFNOQVPJJNP-UHFFFAOYSA-N.cml'
  geom_opt = True
  basis_functions = None

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
        expected_block += """
  cube_stem = 'XLYOFNOQVPJJNP-UHFFFAOYSA-N'
  method_type = 'dft'
  epsiso_param_dict_list = [{'padding': 2.0, 'step_size': 0.088, 'iso_surf': 0.002, 'tol': 3e-05}, {'padding': 2.0, 'step_size': 0.088, 'iso_surf': 0.001, 'tol': 1.5e-05}]
  axes_of_rotation = [np.array([0, 0, 0]), np.array([1, 0, 0]), np.array([0, 1, 0]), np.array([0, 0, 1])]
  angles_in_degrees = [0, 45, 120]

  cube_file_list = property_calc_util.multiOrientationCalcPropEPSISOList(cube_stem,
                                                                         method_type,
                                                                         epsiso_param_dict_list,
                                                                         axes_of_rotation,
                                                                         angles_in_degrees)

  iso_surf_list = [0.002, 0.001]

  if nwchem.ga_nodeid() == 0:
      cube_file_merger.readRotateMergeWriteCubeFileDict(cube_stem, cube_file_list,
                                                        iso_surf_list,
                                                        angles_in_degrees,
                                                        axes_of_rotation)
"""
        expected_block += "end\n"
        epsiso_data = {
            "epsiso_param_dict_list": self.epsiso_prop_dict_list,
            "axes_of_rotation": self.axes_of_rotation,
            "angles_in_degrees": self.angles_in_degrees,
        }
        actual_block = NWInPythonSec.setPythonBlock(
            self.molecule_cml,
            "dft",
            "driverinitial",
            "geometry",
            "XLYOFNOQVPJJNP-UHFFFAOYSA-N",
            epsiso_calc=True,
            epsiso_data=epsiso_data,
        )
        self.assertEqual(actual_block, expected_block)
        with self.assertRaises(TypeError) as err:
            actual_block = NWInPythonSec.setPythonBlock(
                self.molecule_cml,
                "dft",
                "driverinitial",
                "geometry",
                "XLYOFNOQVPJJNP-UHFFFAOYSA-N",
                epsiso_calc=True,
            )
        expected_args = "no epsiso_data present"
        actual_args = err.exception.args[0]
        self.assertEqual(actual_args, expected_args)

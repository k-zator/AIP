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
Script to test the NWChem input generation process.

@author: mark
"""

import unittest
import logging
import numpy as np
import os
import pathlib
from lxml import etree
import nwchemcmlutils.nwUtils.cml_reading_util as cml_reading_util
import nwchemcmlutils.nwinFileGen.NWChemInputGenerator as NWInGen


logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.WARN)


class NWChemInputGeneratorTestCase(unittest.TestCase):
    """Test class for the NWChem input generator script.
    """

    def setUp(self):
        """setUp for the class. Test molecule is water to match the example
        calculation.
        """
        self.memory_limit = 8
        self.molecule_name = "water_example"
        self.title = "Example input script for geometry optimisation"
        self.element_set = {"O", "H"}
        self.parent_directory = pathlib.Path(__file__).parents[1]
        water_cml_file = etree.parse(
            (self.parent_directory / "resources/examplecml.xml").absolute().as_posix()
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
            + 'rg/schema" xmlns:ssip="http://www-hunter.ch.cam.ac.uk/SSIP" ssip:stdInChIKey="XLYOFNOQVPJJNP-UHFFFAO'
            + 'YSA-N" ssip:StuartId="13_01_20121_new" cml:id="XLYOFNOQ'
            + 'VPJJNP-UHFFFAOYSA-N">\n   <cml:atomArray>\n    <cm'
            + 'l:atom cml:elementType="O" cml:id="a1" cml:x3="-0.198000" '
            + 'cml:y3="0.000000" cml:z3="0.347000"/>\n    <cml:atom cml:ele'
            + 'mentType="H" cml:id="a2" cml:x3="0.760000" cml:y3="0.00000'
            + '0" cml:z3="0.204000"/>\n    <cml:atom cml:elementType="H'
            + '" cml:id="a3" cml:x3="-0.561000" cml:y3="0.000000" cml:z3="-0.'
            + '551000"/>\n   </cml:atomArray>\n   <cml:bondArray>'
            + '\n    <cml:bond cml:atomRefs2="a1 a2" cml:order="1"/>\n '
            + '   <cml:bond cml:atomRefs2="a1 a3" cml:order="1"/>\n   </c'
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
        if os.path.isfile("actual_input.nwin"):
            os.remove("actual_input.nwin")
        if os.path.isfile("XLYOFNOQVPJJNP-UHFFFAOYSA-N.nwin"):
            os.remove("XLYOFNOQVPJJNP-UHFFFAOYSA-N.nwin")

    def test_write_nwin_file_contents(self):
        """test to see if total file text is returned as expected.
        """
        expected_file = """memory 8 gb
scratch_dir /tmp
start XLYOFNOQVPJJNP-UHFFFAOYSA-N
title "calculation"
charge 0
"""
        expected_file += """# DFT module
# http://www.nwchem-sw.org/index.php/Density_Functional_Theory_for_Molecules
dft
   # Calculate all integrals "on-the-fly"
   direct
   
   #Exchange functional- B3LYP is default
   # B3LYP Combined Exchange and Correlation Functional
   xc b3lyp

   # Do not print the MO vector coefficients; just too much data.
   noprint "final vectors analysis"

   # Set the energy convergence to be at 1e-08, SCF=Tight in Gaussian
   convergence energy 1e-08

   #now we increase the number of iterations in case there is a bad initial guess
   iterations 100

   # Multiplicity
   mult 1

end
# Module for the geometry optimisation
# http://www.nwchem-sw.org/index.php/Geometry_Optimization#Geometry_Optimization_with_DRIVER
driver
  # Maximum number of steps allowed for the geometry optimisation
  maxiter 600
end
"""
        expected_file += """python
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
        expected_file += """  molecule_cml_string = '''<cml:molecule xmlns:cml="http://www.xml-cml.org/schema" xmlns:ssip="http://www-hunter.ch.cam.ac.uk/SSIP" ssip:stdInChIKey="XLYOFNOQVPJJNP-UHFFFAOYSA-N" ssip:StuartId="13_01_20121_new"
  cml:id="XLYOFNOQVPJJNP-UHFFFAOYSA-N">    <cml:atomArray>     <cml:atom cml:elementType="O" cml:id="a1" cml:x3="-0.198000" cml:y3="0.000000" cml:z3="0.347000"/>     <cml:atom cml:elementType="H"
  cml:id="a2" cml:x3="0.760000" cml:y3="0.000000" cml:z3="0.204000"/>     <cml:atom cml:elementType="H" cml:id="a3" cml:x3="-0.561000" cml:y3="0.000000" cml:z3="-0.551000"/>    </cml:atomArray>
  <cml:bondArray>     <cml:bond cml:atomRefs2="a1 a2" cml:order="1"/>     <cml:bond cml:atomRefs2="a1 a3" cml:order="1"/>    </cml:bondArray>   </cml:molecule>'''
"""
        expected_file += """  method_type = 'dft'
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
        expected_file += "end\n"
        expected_file += "task python\n"
        actual_file = NWInGen.writeNwinFileContents(
            self.memory_limit,
            self.molecule_cml,
            "dft",
            "driverinitial",
            "geometry",
            self.filename_stem,
            molecule_name=self.filename_stem,
        )
        self.assertMultiLineEqual(actual_file, expected_file)

    def test_write_nwin_epsiso_file_contents(self):
        """test to see if total file text is returned as expected.
        """
        expected_file = """memory 8 gb
scratch_dir /tmp
start XLYOFNOQVPJJNP-UHFFFAOYSA-N
title "calculation"
charge 0
"""
        expected_file += """# DFT module
# http://www.nwchem-sw.org/index.php/Density_Functional_Theory_for_Molecules
dft
   # Calculate all integrals "on-the-fly"
   direct
   
   #Exchange functional- B3LYP is default
   # B3LYP Combined Exchange and Correlation Functional
   xc b3lyp

   # Do not print the MO vector coefficients; just too much data.
   noprint "final vectors analysis"

   # Set the energy convergence to be at 1e-08, SCF=Tight in Gaussian
   convergence energy 1e-08

   #now we increase the number of iterations in case there is a bad initial guess
   iterations 100

   # Multiplicity
   mult 1

end
# Module for the geometry optimisation
# http://www.nwchem-sw.org/index.php/Geometry_Optimization#Geometry_Optimization_with_DRIVER
driver
  # Maximum number of steps allowed for the geometry optimisation
  maxiter 600
end
"""
        expected_file += """python
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
        expected_file += u"""  molecule_cml_string = '''<cml:molecule xmlns:cml="http://www.xml-cml.org/schema" xmlns:ssip="http://www-hunter.ch.cam.ac.uk/SSIP" ssip:stdInChIKey="XLYOFNOQVPJJNP-UHFFFAOYSA-N" ssip:StuartId="13_01_20121_new"
  cml:id="XLYOFNOQVPJJNP-UHFFFAOYSA-N">    <cml:atomArray>     <cml:atom cml:elementType="O" cml:id="a1" cml:x3="-0.198000" cml:y3="0.000000" cml:z3="0.347000"/>     <cml:atom cml:elementType="H"
  cml:id="a2" cml:x3="0.760000" cml:y3="0.000000" cml:z3="0.204000"/>     <cml:atom cml:elementType="H" cml:id="a3" cml:x3="-0.561000" cml:y3="0.000000" cml:z3="-0.551000"/>    </cml:atomArray>
  <cml:bondArray>     <cml:bond cml:atomRefs2="a1 a2" cml:order="1"/>     <cml:bond cml:atomRefs2="a1 a3" cml:order="1"/>    </cml:bondArray>   </cml:molecule>'''
"""
        expected_file += """  method_type = 'dft'
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
        expected_file += """
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
        expected_file += "end\n"
        expected_file += "task python\n"
        epsiso_data = {
            "epsiso_param_dict_list": self.epsiso_prop_dict_list,
            "axes_of_rotation": self.axes_of_rotation,
            "angles_in_degrees": self.angles_in_degrees,
        }
        actual_file = NWInGen.writeNwinFileContents(
            self.memory_limit,
            self.molecule_cml,
            "dft",
            "driverinitial",
            "geometry",
            self.filename_stem,
            epsiso_calc=True,
            epsiso_data=epsiso_data,
            molecule_name=self.filename_stem,
        )
        self.assertEqual(actual_file, expected_file)

    def test_write_nwin_contents_to_file(self):
        """Test to check that expected file is written.
        """
        expected_filename = (
            (self.parent_directory / "resources/expected_input.nwin")
            .absolute()
            .as_posix()
        )
        actual_filename = "actual_input.nwin"
        nwin_contents = NWInGen.writeNwinFileContents(
            self.memory_limit,
            self.molecule_cml,
            "dft",
            "driverinitial",
            "geometry",
            self.filename_stem,
            molecule_name=self.filename_stem,
        )
        NWInGen.writeNwinFileContentsToFile(actual_filename, nwin_contents)
        with open(expected_filename, "r") as expected_file:
            with open(actual_filename, "r") as actual_file:
                actual_file_read_in = actual_file.read()
                expected_file_read_in = expected_file.read()
                self.assertMultiLineEqual(actual_file_read_in, expected_file_read_in)

    def test_write_nwin_file(self):
        """Test to check that the expected file is written.
        """
        self.maxDiff = None
        expected_filename = (
            (self.parent_directory / "resources/expected_input.nwin")
            .absolute()
            .as_posix()
        )
        actual_filename = "XLYOFNOQVPJJNP-UHFFFAOYSA-N.nwin"
        NWInGen.writeNwinFile(
            "XLYOFNOQVPJJNP-UHFFFAOYSA-N",
            self.memory_limit,
            self.molecule_cml,
            "dft",
            "driverinitial",
            "geometry",
            directory="",
            molecule_name=self.filename_stem,
        )
        with open(expected_filename, "r") as expected_file:
            with open(actual_filename, "r") as actual_file:
                actual_file_read_in = actual_file.read()
                expected_file_read_in = expected_file.read()
                self.assertMultiLineEqual(actual_file_read_in, expected_file_read_in)

    def test_write_epsiso_file(self):
        """Test to see if a file is produced with the expected contents for
        epsiso job.
        """
        self.maxDiff = None
        expected_filename = (
            (self.parent_directory / "resources/expected_input2.nwin")
            .absolute()
            .as_posix()
        )
        actual_filename = "XLYOFNOQVPJJNP-UHFFFAOYSA-N.nwin"
        epsiso_data = {
            "epsiso_param_dict_list": self.epsiso_prop_dict_list,
            "axes_of_rotation": self.axes_of_rotation,
            "angles_in_degrees": self.angles_in_degrees,
        }
        NWInGen.writeNwinFile(
            "XLYOFNOQVPJJNP-UHFFFAOYSA-N",
            self.memory_limit,
            self.molecule_cml,
            "dft",
            "driverinitial",
            "geometry",
            epsiso_calc=True,
            epsiso_data=epsiso_data,
            directory="",
            molecule_name=self.filename_stem,
        )
        with open(expected_filename, "r") as expected_file:
            with open(actual_filename, "r") as actual_file:
                actual_file_read_in = actual_file.read()
                expected_file_read_in = expected_file.read()
                self.assertMultiLineEqual(actual_file_read_in, expected_file_read_in)

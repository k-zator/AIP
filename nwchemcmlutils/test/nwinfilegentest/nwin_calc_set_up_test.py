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
Script to test the NWChem input generation process. This tests the creation of
the sections for reading in by the NWChem input parser- this is the information
for setting up the calculation.

@author: mark
"""

import unittest
import logging
import pathlib
from lxml import etree
import nwchemcmlutils.nwUtils.cml_reading_util as cml_reading_util
import nwchemcmlutils.nwinFileGen.NWChemCalcSetUp as NWInCalc


logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.WARN)


class NWChemCalcSetUpTestCase(unittest.TestCase):
    """Test class for the testing of the NWChem calculation set up scirpt
    functions.
    """

    def setUp(self):
        """setUp for the class. Test molecule is water to match the example
        calculation.
        """
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

    def test_write_memory_line(self):
        """tests to see if the correct memory line is written.
        """
        expected_line = "memory 8 gb\n"
        actual_line = NWInCalc.writeMemoryLine(self.memory_limit)
        self.assertEqual(actual_line, expected_line)

    def test_set_scratch_dir(self):
        """test to see if the correct scratch dir line is returned.
        """
        expected_line = "scratch_dir /tmp\n"
        actual_line = NWInCalc.setScratchDir("/tmp")
        self.assertEqual(actual_line, expected_line)

    def test_set_start(self):
        """test to see if the correct start line is written.
        """
        expected_line = "start water_example\n"
        actual_line = NWInCalc.setStart(molecule_name=self.molecule_name)
        self.assertEqual(actual_line, expected_line)

    def test_set_title(self):
        """test to see if the correct title line is written.
        """
        expected_line = 'title "Example input script for geometry optimisation"\n'
        actual_line = NWInCalc.setTitle(title=self.title)
        self.assertEqual(actual_line, expected_line)

    def test_set_charge(self):
        """test to see if the correct charge line is written.
        """
        expected_line = "charge 0\n"
        actual_line = NWInCalc.setCharge(0)
        self.assertEqual(actual_line, expected_line)

    def test_set_dft(self):
        """test to see if th correct dft block is set.
        """
        expected_line = """# DFT module
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
"""
        actual_line = NWInCalc.setDFT(1, "b3lyp", 100)
        self.assertEqual(actual_line, expected_line)

    def test_set_driver(self):
        """test to see if correct driver section is written.
        """
        expected_line = """# Module for the geometry optimisation
# http://www.nwchem-sw.org/index.php/Geometry_Optimization#Geometry_Optimization_with_DRIVER
driver
  # Maximum number of steps allowed for the geometry optimisation
  maxiter 200
end
"""
        actual_line = NWInCalc.setGeometryDriver(200)
        self.assertEqual(actual_line, expected_line)

    def test_set_basis(self):
        """test to see if correct basis block is written.
        """
        expected_basis_1 = """basis
H library 6-31G*
O library 6-31G*
end
"""
        actual_basis_1 = NWInCalc.setBasis(self.element_set)
        self.assertMultiLineEqual(actual_basis_1, expected_basis_1)
        # test to see if extra definitions are handled correctly.
        expected_basis_2 = """basis
H library 6-31G*
O library 6-311G**
end
"""
        actual_basis_2 = NWInCalc.setBasis(self.element_set, **{"O": "6-311G**"})
        self.assertEqual(actual_basis_2, expected_basis_2)

    def test_set_property_epsiso(self):
        """Test to see if expected string is returned.
        """
        expected_block = """
# Property module
# http://www.nwchem-sw.org/index.php/Properties
# New property keyword in use
property
  grid pad 2.000000 step 0.088000
  espiso iso 0.002000 tol 0.000030
end
"""
        actual_block = NWInCalc.setPropertyBlockEPSISO(2.0, 0.088, 0.002, 0.00003)
        self.assertMultiLineEqual(actual_block, expected_block)

    def test_set_geometry(self):
        """test to see if correct geometry block is written
        """
        expected_line = """geometry units an autoz nocenter noautosym
        O       -0.198000        0.000000        0.347000
        H        0.760000        0.000000        0.204000
        H       -0.561000        0.000000       -0.551000
end
"""
        actual_line = NWInCalc.setGeometry(self.atom_line_list)
        self.assertEqual(actual_line, expected_line)

    def test_write_nwin_set_up_lines(self):
        """test to see that correct lines are generated for the start of the
        nwin file.
        """
        expected_lines = """memory 8 gb
scratch_dir /tmp
start molecule
title "calculation"
charge 0
"""
        actual_lines = NWInCalc.writeNwinSetUpLines(self.memory_limit)
        self.assertEqual(actual_lines, expected_lines)

    def test_write_dft_driver_lines(self):
        """test to see if expected dft and driver blocks are created.
        """
        expected_blocks = """# DFT module
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
        actual_blocks = NWInCalc.writeDFTAndDriverLines()
        self.assertEqual(actual_blocks, expected_blocks)

    def test_set_task(self):
        """test to see expected task line is written.
        """
        expected_line = "task python\n"
        actual_line = NWInCalc.setTask("python")
        self.assertEqual(actual_line, expected_line)

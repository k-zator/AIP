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
Script contains tests for geom_opt_cml.

@author: mark
"""

import unittest
import logging
import os
import pathlib
import unittest.mock as mock
from lxml import etree
import nwchemcmlutils.nwUtils.geom_opt_cml as geom_opt_cml
import nwchemcmlutils.nwUtils.rtdb_util as rtdb_util

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.WARN)


class GeomOptCMLTestCase(unittest.TestCase):
    """Test class for the testing of the Geom opt calc functions.
    """

    def setUp(self):
        """setUp for the class. Test molecule is water to match the example
        calculation.
        """
        self.maxDiff = None
        parent_directory = pathlib.Path(__file__).parents[1]
        self.expected_filename = (
            (parent_directory / "resources/water_exp_output.cml").absolute().as_posix()
        )
        self.out_filename = (
            (parent_directory / "resources/water_out_example.cml").absolute().as_posix()
        )
        self.molecule_cml = (
            u'<cml:molecule xmlns:cml="http://www.xml-cml.o'
            + 'rg/schema" xmlns:ssip="http://www-hunter.ch.cam.ac.uk/SSIP" ssip:stdInChIKey="XLYOFNOQVPJJNP-UHFFFAO'
            + 'YSA-N" StuartId="13_01_20121_new" cml:id="XLYOFNOQ'
            + 'VPJJNP-UHFFFAOYSA-N">\n  <cml:atomArray>\n  <cm'
            + 'l:atom cml:elementType="O" cml:id="a1" cml:x3="-0.198000" '
            + 'cml:y3="0.000000" cml:z3="0.347000"/>\n  <cml:atom cml:ele'
            + 'mentType="H" cml:id="a2" cml:x3="0.760000" cml:y3="0.00000'
            + '0" cml:z3="0.204000"/>\n  <cml:atom cml:elementType="H'
            + '" cml:id="a3" cml:x3="-0.561000" cml:y3="0.000000" cml:z3="-0.'
            + '551000"/>\n  </cml:atomArray>\n  <cml:bondArray>'
            + '\n  <cml:bond cml:atomRefs2="a1 a2" cml:order="1"/>\n '
            + ' <cml:bond cml:atomRefs2="a1 a3" cml:order="1"/>\n  </c'
            + "ml:bondArray>\n  </cml:molecule>\n  "
        )
        example_cml = (
            (parent_directory / "resources/examplecml.xml").absolute().as_posix()
        )
        water_cml_file = etree.parse(example_cml)
        cml_namespace = {"cml": "http://www.xml-cml.org/schema"}
        self.water_cml = water_cml_file.findall(
            "cml:molecule", namespaces=cml_namespace
        )[0]
        self.get_values_1 = [
            rtdb_util.NWChemError,
            ["O", "H", "H"],
            rtdb_util.NWChemError,
            [
                -3.00686549e-02,
                -2.32489942e-29,
                7.00333940e-01,
                1.72635129e00,
                2.61236666e-29,
                1.84332816e-01,
                -9.46061415e-01,
                1.98520447e-30,
                -8.84666757e-01,
            ],
            "angstroms",
            1.8897259885789233,
        ]

    def tearDown(self):
        """Tear down for the class.
        """
        del self.molecule_cml
        del self.water_cml
        if os.path.isfile("water_test_output.cml"):
            os.remove("water_test_output.cml")

    def test_read_set_geom_basis(self):
        """Test to see if geometry and basis are  set as expected, and expected
        molecule is returned.
        """
        # set up the input_parse mock
        input_parse_mock = mock.Mock()
        geom_opt_cml.nwchem.attach_mock(input_parse_mock, "input_parse")
        actual_mol_cml = geom_opt_cml.readCMLSetGeomAndBasis(self.molecule_cml)
        # check input_parse calls.
        actual_call_arg_list = [
            str(arg_list) for arg_list in geom_opt_cml.nwchem.input_parse.call_args_list
        ]
        expected_call_arg_list = [
            (
                "call('geometry units an autoz nocenter n"
                + "oautosym\\n        O       -0.198000       "
                + " 0.000000        0.347000\\n        H      "
                + "  0.760000        0.000000        0.204000"
                + "\\n        H       -0.561000        0.0000"
                + "00       -0.551000\\nend\\n')"
            ),
            "call('basis\\nH library 6-31G*\\nO library 6-31G*\\nend\\n')",
        ]
        self.assertListEqual(actual_call_arg_list, expected_call_arg_list)
        expected_mol_cml = (
            u'<cml:molecule xmlns:cml="http://www.xml-cml.o'
            + 'rg/schema" xmlns:ssip="http://www-hunter.ch.cam.ac.uk/SSIP" ssip:stdInChIKey="XLYOFNOQVPJJNP-UHFFFAO'
            + 'YSA-N" StuartId="13_01_20121_new" cml:id="XLYOFNOQ'
            + 'VPJJNP-UHFFFAOYSA-N">\n  <cml:atomArray>\n  <c'
            + 'ml:atom cml:elementType="O" cml:id="a1" cml:x3="-0.198000" '
            + 'cml:y3="0.000000" cml:z3="0.347000"/>\n  <cml:atom cml:elem'
            + 'entType="H" cml:id="a2" cml:x3="0.760000" cml:y3="0.000000" '
            + 'cml:z3="0.204000"/>\n  <cml:atom cml:elementType="H" cml:id='
            + '"a3" cml:x3="-0.561000" cml:y3="0.000000" cml:z3="-0.551000"'
            + "/>\n  </cml:atomArray>\n  <cml:bondArray>\n  <cm"
            + 'l:bond cml:atomRefs2="a1 a2" cml:order="1"/>\n  <cml:bond'
            + ' cml:atomRefs2="a1 a3" cml:order="1"/>\n  </cml:bondArr'
            + "ay>\n  </cml:molecule>"
        )

        self.assertEqual(etree.tounicode(actual_mol_cml), expected_mol_cml)

    def test_geom_opt_calc(self):
        """Test to see if input_parse is called with the expected arguements.
        """
        task_opt_mock = mock.Mock()
        geom_opt_cml.nwchem.attach_mock(task_opt_mock, "task_optimize")
        geom_opt_cml.geomOptCalc("dft")
        actual_call_arg_list_1 = [
            str(arg_list)
            for arg_list in geom_opt_cml.nwchem.task_optimize.call_args_list
        ]
        expected_call_arg_list_1 = ["call('dft')"]
        self.assertListEqual(actual_call_arg_list_1, expected_call_arg_list_1)

    def test_geom_extract_to_cml(self):
        """Test to see if geometry is correctly extracted, and written.
        """
        self.maxDiff = None
        rtdb_get_mock = mock.Mock(side_effect=self.get_values_1)
        geom_opt_cml.rtdb_util.nwchem.attach_mock(rtdb_get_mock, "rtdb_get")
        ga_nodeid_mock = mock.Mock(return_value=0)
        geom_opt_cml.nwchem.attach_mock(ga_nodeid_mock, "ga_nodeid")
        actual_mol_cml = geom_opt_cml.extractCoordWriteToCML(
            self.water_cml, "driverinitial", "geometry"
        )
        actual_call_arg_list_1 = [
            str(arg_list)
            for arg_list in geom_opt_cml.rtdb_util.nwchem.rtdb_get.call_args_list
        ]
        expected_call_arg_list_1 = [
            "call('driverinitial')",
            "call('geometry:driverinitial:tags')",
            "call('geometry')",
            "call('geometry:geometry:coords')",
            "call('geometry:geometry:user units')",
            "call('geometry:geometry:angstrom_to_au')",
        ]
        self.assertListEqual(actual_call_arg_list_1, expected_call_arg_list_1)
        actual_call_arg_list_2 = [
            str(arg_list) for arg_list in geom_opt_cml.nwchem.ga_nodeid.call_args_list
        ]
        expected_call_arg_list_2 = ["call()"]
        self.assertListEqual(actual_call_arg_list_2, expected_call_arg_list_2)
        expected_water_cml = u"""<cml:molecule xmlns:cml="http://www.xml-cml.org/schema" xmlns:ssip="http://www-hunter.ch.cam.ac.uk/SSIP" ssip:stdInChIKey="XLYOFNOQVPJJNP-UHFFFAOYSA-N" StuartId="13_01_20121_new" cml:id="XLYOFNOQVPJJNP-UHFFFAOYSA-N">
 <cml:atomArray>
  <cml:atom cml:elementType="O" cml:id="a1" cml:x3="-0.015911648081112367" cml:y3="-1.2302838792772956e-29" cml:z3="0.37060078775053107"/>
  <cml:atom cml:elementType="H" cml:id="a2" cml:x3="0.9135458264498012" cml:y3="1.3824050025181182e-29" cml:z3="0.09754473247130319"/>
  <cml:atom cml:elementType="H" cml:id="a3" cml:x3="-0.5006341769747473" cml:y3="1.050525040137103e-30" cml:z3="-0.46814552075101146"/>
 </cml:atomArray>
 <cml:bondArray>
  <cml:bond cml:atomRefs2="a1 a2" cml:order="1"/>
  <cml:bond cml:atomRefs2="a1 a3" cml:order="1"/>
 </cml:bondArray>
</cml:molecule>
"""
        self.assertMultiLineEqual(etree.tounicode(actual_mol_cml), expected_water_cml)

    def test_write_cml_to_file(self):
        """Test to see if geometry is written to file as expected.
        """
        ga_nodeid_mock = mock.Mock(return_value=0)
        geom_opt_cml.nwchem.attach_mock(ga_nodeid_mock, "ga_nodeid")
        actual_filename = "water_test_output.cml"
        file_out = geom_opt_cml.writeCMLToFile(actual_filename, self.water_cml)
        actual_call_arg_list_1 = [
            str(arg_list) for arg_list in geom_opt_cml.nwchem.ga_nodeid.call_args_list
        ]
        expected_call_arg_list_1 = ["call()"]
        self.assertListEqual(actual_call_arg_list_1, expected_call_arg_list_1)
        expected_out = 0
        self.assertEqual(file_out, expected_out)
        expected_filename = self.out_filename
        with open(expected_filename, "r") as expected_file:
            with open(actual_filename, "r") as actual_file:
                actual_file_read_in = actual_file.read()
                expected_file_read_in = expected_file.read()
                self.assertEqual(actual_file_read_in, expected_file_read_in)

    def test_read_cml_opt_write_file(self):
        """Test for function which reads in molecule_cml_string, extracts the
        information and puts it in the rtdb, carries out a geometry optimisation,
        then extracts the output coordinates and writes those to file.
        """
        LOGGER.debug("Setting up mocks")
        input_parse_mock = mock.Mock()
        geom_opt_cml.nwchem.attach_mock(input_parse_mock, "input_parse")
        task_opt_mock = mock.Mock()
        geom_opt_cml.nwchem.attach_mock(task_opt_mock, "task_optimize")
        rtdb_get_mock = mock.Mock(side_effect=self.get_values_1)
        geom_opt_cml.rtdb_util.nwchem.attach_mock(rtdb_get_mock, "rtdb_get")
        ga_nodeid_mock = mock.Mock(return_value=0)
        geom_opt_cml.nwchem.attach_mock(ga_nodeid_mock, "ga_nodeid")

        actual_filename = "water_test_output.cml"
        file_out = geom_opt_cml.readCMLOptAndWriteOutputToFile(
            self.molecule_cml, "dft", "driverinitial", "geometry", actual_filename
        )
        actual_call_arg_list_1 = [
            str(arg_list) for arg_list in geom_opt_cml.nwchem.input_parse.call_args_list
        ]
        expected_call_arg_list_1 = [
            (
                "call('geometry units an autoz nocenter n"
                + "oautosym\\n        O       -0.198000       "
                + " 0.000000        0.347000\\n        H      "
                + "  0.760000        0.000000        0.204000"
                + "\\n        H       -0.561000        0.0000"
                + "00       -0.551000\\nend\\n')"
            ),
            "call('basis\\nH library 6-31G*\\nO library 6-31G*\\nend\\n')",
        ]
        self.assertListEqual(actual_call_arg_list_1, expected_call_arg_list_1)
        actual_call_arg_list_2 = [
            str(arg_list)
            for arg_list in geom_opt_cml.nwchem.task_optimize.call_args_list
        ]
        expected_call_arg_list_2 = ["call('dft')"]
        self.assertListEqual(actual_call_arg_list_2, expected_call_arg_list_2)
        actual_call_arg_list_3 = [
            str(arg_list)
            for arg_list in geom_opt_cml.rtdb_util.nwchem.rtdb_get.call_args_list
        ]
        expected_call_arg_list_3 = [
            "call('driverinitial')",
            "call('geometry:driverinitial:tags')",
            "call('geometry')",
            "call('geometry:geometry:coords')",
            "call('geometry:geometry:user units')",
            "call('geometry:geometry:angstrom_to_au')",
        ]
        self.assertListEqual(actual_call_arg_list_3, expected_call_arg_list_3)
        actual_call_arg_list_4 = [
            str(arg_list) for arg_list in geom_opt_cml.nwchem.ga_nodeid.call_args_list
        ]
        expected_call_arg_list_4 = ["call()", "call()"]
        self.assertListEqual(actual_call_arg_list_4, expected_call_arg_list_4)
        expected_out = 0
        self.assertEqual(file_out, expected_out)
        expected_filename = self.expected_filename
        with open(expected_filename, "r") as expected_file:
            with open(actual_filename, "r") as actual_file:
                actual_file_read_in = actual_file.read()
                expected_file_read_in = expected_file.read()
                self.assertMultiLineEqual(actual_file_read_in, expected_file_read_in)

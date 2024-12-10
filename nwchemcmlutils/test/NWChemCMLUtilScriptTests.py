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
Script to operate all tests.

@author: mark
"""

import unittest
import logging
import sys
from nwchemcmlutils.test.nwutilstest.cml_reading_util_test import CMLReadingUtilTestCase
from nwchemcmlutils.test.nwutilstest.cml_writing_util_test import CMLWritingUtilTestCase
from nwchemcmlutils.test.nwutilstest.trig_test import TrigTestCase
from nwchemcmlutils.test.nwutilstest.rtdb_util_test import RTDBUtilTestCase
from nwchemcmlutils.test.nwutilstest.property_calc_util_test import PropertyCalcUtilTestCase
from nwchemcmlutils.test.nwutilstest.cube_file_merger_test import CubeFileMergerTestCase
from nwchemcmlutils.test.nwutilstest.geom_opt_cml_test import GeomOptCMLTestCase
from nwchemcmlutils.test.dataclassestest.gaussian_cube_overload_methods_test import (
    GaussianCubeOverloadMethodsTestCase,
)
from nwchemcmlutils.test.dataclassestest.gaussian_cube_read_write_test import GaussianCubeReadWriteTestCase
from nwchemcmlutils.test.dataclassestest.gaussian_cube_rotate_add_plot_test import (
    GaussianCubeRotateAddPlotTestCase,
)
from nwchemcmlutils.test.nwinfilegentest.nwin_calc_set_up_test import NWChemCalcSetUpTestCase
from nwchemcmlutils.test.nwinfilegentest.nwin_python_sec_test import NWChemPythonSecTestCase
from nwchemcmlutils.test.nwinfilegentest.nwin_input_generator_test import NWChemInputGeneratorTestCase
from nwchemcmlutils.test.nwinfilegentest.ziggy_submit_generator_test import ZiggySubmitGeneratorTestCase
from nwchemcmlutils.test.nwinfilegentest.slurmsubmitgeneratortest import SLURMSubmitGeneratorTestCase
from nwchemcmlutils.test.dataclassestest.distance_unit_test import DistanceUnitTestCase
from nwchemcmlutils.test.dataclassestest.cartesian_3d_test import Cartesian3DTestCase
from nwchemcmlutils.test.dataclassestest.vector_3d_test import Vector3DTestCase
from nwchemcmlutils.test.dataclassestest.atom_test import AtomTestCase
from nwchemcmlutils.test.dataclassestest.property_value_test import PropertyValueTestCase
from nwchemcmlutils.test.dataclassestest.molecule_test import MoleculeTestCase
from nwchemcmlutils.test.dataclassestest.property_value_list_test import PropertyValueListTestCase
from nwchemcmlutils.test.calccreatortest.cml_file_reader_test import CMLFileReaderTestCase
from nwchemcmlutils.test.calccreatortest.job_file_creator_test import JobFileCreatorTestCase
from nwchemcmlutils.test.calccreatortest.directory_creator_test import DirectoryCreatorTestCase
from nwchemcmlutils.test.calccreatortest.submit_script_creator_test import SubmitScriptCreatorTestCase
from nwchemcmlutils.test.epsisocalccreatortest import EPSISOCalcCreatorTestCase

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.WARN)

DATACLASSESTESTCASES = [
    Vector3DTestCase,
    DistanceUnitTestCase,
    Cartesian3DTestCase,
    AtomTestCase,
    PropertyValueTestCase,
    MoleculeTestCase,
    PropertyValueListTestCase,
    GaussianCubeOverloadMethodsTestCase,
    GaussianCubeReadWriteTestCase,
    GaussianCubeRotateAddPlotTestCase,
]

NWUTILSTESTCASES = [
    TrigTestCase,
    RTDBUtilTestCase,
    CMLReadingUtilTestCase,
    CMLWritingUtilTestCase,
    GeomOptCMLTestCase,
    PropertyCalcUtilTestCase,
    CubeFileMergerTestCase,
]

NWINFILEGENTESTCASES = [
    NWChemCalcSetUpTestCase,
    NWChemPythonSecTestCase,
    NWChemInputGeneratorTestCase,
    ZiggySubmitGeneratorTestCase,
    SLURMSubmitGeneratorTestCase,
]
CALCCREATORTESTCASES = [
    CMLFileReaderTestCase,
    JobFileCreatorTestCase,
    DirectoryCreatorTestCase,
    SubmitScriptCreatorTestCase,
]

MODULETESTCASES = [
    EPSISOCalcCreatorTestCase,
]

def create_test_suite():
    """Function creates a test suite and then loads all the tests from the
    different test cases.
    """
    LOGGER.info("setting up loader and test suite")
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    for test_case in DATACLASSESTESTCASES:
        LOGGER.debug("Adding %s", test_case)
        suite.addTests(loader.loadTestsFromTestCase(test_case))
    for test_case in NWUTILSTESTCASES:
        LOGGER.debug("Adding %s", test_case)
        suite.addTests(loader.loadTestsFromTestCase(test_case))
    for test_case in NWINFILEGENTESTCASES:
        LOGGER.debug("Adding %s", test_case)
        suite.addTests(loader.loadTestsFromTestCase(test_case))
    for test_case in CALCCREATORTESTCASES:
        LOGGER.debug("Adding %s", test_case)
        suite.addTests(loader.loadTestsFromTestCase(test_case))
    for test_case in MODULETESTCASES:
        LOGGER.debug("Adding %s", test_case)
        suite.addTests(loader.loadTestsFromTestCase(test_case))
    return suite


def run_tests():
    """Run test suite. Exits if there is a failure.

    Returns
    -------
    None
    """
    LOGGER.info("calling test suite method")
    suite = create_test_suite()
    LOGGER.info("running test suite")
    ret = (
        not unittest.TextTestRunner(verbosity=2, stream=sys.stderr)
        .run(suite)
        .wasSuccessful()
    )
    if ret:
        sys.exit(ret)


if __name__ == "__main__":
    run_tests()

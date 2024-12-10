# -*- coding: utf-8 -*-
#    cmlgenerator creates fully namespaced CML for molecules from input structures.
#    Copyright (C) 2019  Mark D. Driver
#
#    cmlgenerator is free software: you can redistribute it and/or modify
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
"""
Script for running tests.

@author: mark
"""

import unittest
import logging
import sys
from tests.structuretocmltest import StructureToCMLTestCase
from tests.smilestocmltest import SmilesToCMLTestCase
from tests.cmlnamespacingtest import CMLNamespacingTestCase
from tests.cmlmakertest import CMLMakerTestCase

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.WARN)

TESTCASES = [
    StructureToCMLTestCase,
    SmilesToCMLTestCase,
    CMLNamespacingTestCase,
    CMLMakerTestCase,
]


def create_test_suite():
    """Function creates a test suite and then loads all the tests from the
    different test cases.
    """
    LOGGER.info("setting up loader and test suite")
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    for test_case in TESTCASES:
        LOGGER.debug("Adding %s", test_case)
        suite.addTests(loader.loadTestsFromTestCase(test_case))
    return suite


def run_tests():
    """Runs test suite. Exits if there is a failure.

    Returns
    --------
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

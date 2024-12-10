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
Script containing test case for epsisocalccreator methods.

@author: mark
"""

import logging
import unittest
import pathlib
import nwchemcmlutils.epsisocalccreator as epscalc

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.WARN)

class EPSISOCalcCreatorTestCase(unittest.TestCase):
    """Test case for epsisocalccreator methods"""
    def setUp(self):
        """Set up before tests.

        Returns
        -------
        None.

        """
        self.maxDiff=None
        parent_directory = pathlib.Path(__file__).parents[0]
        self.json_filename = (parent_directory / "resources/exampleparameters.json").absolute().as_posix()
    def test_get_parameters_from_json_file(self):
        """Test expected file contents are read in.

        Returns
        -------
        None.

        """
        expected_dict = {
    "account": "HUNTER-SL3-CPU",
    "partition": "skylake-himem",
    "walltime": "12:00:00",
    "processors": 128,
    "nodes": 4,
    "conda_version": "miniconda2-4.3.14-gcc-5.4.0-xjtq53h",
    "conda_env": "nwchempy",
    "modules": ["rhel7/default-peta4",
             "gcc-5.4.0-gcc-4.8.5-fis24gg",
             "intel/impi/2017.4/gnu"],
    "memory_limit": 12,  # this should be per CPU.
    "paths":["NWCHEM_TOP=/home/mdd31/rds/hpc-work/nwchem-6.8.1-release",
             "PATH=$PATH:$NWCHEM_TOP/bin/LINUX64",
             "PYTHONHOME=~/.conda/envs/nwchempy",
             "LD_LIBRARY_PATH=~/.conda/envs/nwchempy/lib:$LD_LIBRARY_PATH"],
    "npernodesec":"-ppn $mpi_tasks_per_node",
    "machinefile": "export NODEFILE=`generate_pbs_nodefile`\n        cat $NODEFILE | uniq > machine.file.$JOBID",
    "utils_path":"PYTHONPATH=$PYTHONPATH:$NWCHEM_TOP/contrib/python:${CONDA_PREFIX}/lib/python3.11/site-packages",
    "slurm_type": "camhpc",
}
        actual_dict = epscalc.get_parameters_from_json_file(self.json_filename)
        self.assertDictEqual(expected_dict, actual_dict)
        self.assertDictEqual(expected_dict, epscalc.CAM_HPC_DEFAULTS)
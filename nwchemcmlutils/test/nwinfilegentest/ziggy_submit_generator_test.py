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
script containing the tests for the ZiggySubmitGenerator script functions.

@author: mark
"""

import unittest
import logging
import os
import pathlib
import nwchemcmlutils.nwinFileGen.ZiggySubmitGenerator as ZiggyGen

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.WARN)


class ZiggySubmitGeneratorTestCase(unittest.TestCase):
    """Test case for the ZiggySubmitGenerator functions.
    """

    def setUp(self):
        """set up for tests.
        """
        self.parent_directory = pathlib.Path(__file__).parents[1]

    def tearDown(self):
        """tear down for tests.
        """
        if os.path.isfile("actual_ziggy_submit.sh"):
            os.remove("actual_ziggy_submit.sh")
        if os.path.isfile("actual_input2_submit.sh"):
            os.remove("actual_input2_submit.sh")

    def test_write_submit_file_start(self):
        """test to see if start of file is written as expected.
        """
        expected_line = """#!/bin/bash -x
# Start of submission script


/opt/torque/bin/qsub << EOF
"""
        actual_line = ZiggyGen.writeStartOfSubmitScript()
        self.assertEqual(actual_line, expected_line)

    def test_write_pbs_start(self):
        """test to see if start of pbs directives is written as expected.
        """
        expected_line = """
########################
# Start PBS directives #
########################

"""
        actual_line = ZiggyGen.writePBSDirectiveStart()
        self.assertEqual(actual_line, expected_line)

    def test_set_pbs_name(self):
        """test to see if the correct lines are written.
        """
        expected_line = """
#          Set the name of the job (up to 15 characters,
#          no blank spaces, start with alphanumeric character)

#PBS -N default_name
"""
        actual_line = ZiggyGen.setPBSJobName("default_name")
        self.assertEqual(actual_line, expected_line)

    def test_set_pbs_queue(self):
        """test to see if the correct lines are written
        """
        expected_line = """#          Specify the queue.

#PBS -q s16
"""
        actual_line = ZiggyGen.setPBSQueue("s16")
        self.assertEqual(actual_line, expected_line)

    def test_set_wall_node_proc(self):
        """
        """
        expected_line = """
#          Specify the maximum wall clock time.
#          Specify the number of nodes requested and the
#          number of processors per node.


#PBS -l nodes=1:ppn=16,walltime=24:00:00
"""
        actual_line = ZiggyGen.setPBSWalltimeNodesProcessors(1, 16, "24:00:00")
        self.assertEqual(actual_line, expected_line)

    def test_set_pbs_job_output(self):
        """test to see correct job output is returned.
        """
        expected_line = """
#          The directive below directs that the standard output and
#          error streams are to be merged, intermixed, as standard
#          output.

#PBS -j oe
"""
        actual_line = ZiggyGen.setPBSJobOutput()
        self.assertEqual(actual_line, expected_line)

    def test_set_pbs_directive_end(self):
        """test to see end of directives is written as expected.
        """
        expected_line = """
########################
# End PBS directives   #
########################
"""
        actual_line = ZiggyGen.writePBSDirectiveEnd()
        self.assertEqual(actual_line, expected_line)

    def test_write_pbs_section(self):
        """test to see pbs directive section is written as expected.
        """
        expected_section = """
########################
# Start PBS directives #
########################


#          Set the name of the job (up to 15 characters,
#          no blank spaces, start with alphanumeric character)

#PBS -N default_name
#          Specify the queue.

#PBS -q s16

#          Specify the maximum wall clock time.
#          Specify the number of nodes requested and the
#          number of processors per node.


#PBS -l nodes=1:ppn=16,walltime=24:00:00

#          The directive below directs that the standard output and
#          error streams are to be merged, intermixed, as standard
#          output.

#PBS -j oe

########################
# End PBS directives   #
########################
"""
        actual_section = ZiggyGen.writePBSSection()
        self.assertEqual(actual_section, expected_section)

    def test_write_module_purge(self):
        """test to check that correct line is returned.
        """
        expected_line = "\nmodule purge\n"
        actual_line = ZiggyGen.writeModulePurge()
        self.assertEqual(actual_line, expected_line)

    def test_anaconda_set_up(self):
        """test to see expected anaconda version and environment is set.
        """
        expected_line = """
module load anaconda/python2/2.4.1
source activate nwchempy

"""
        actual_line = ZiggyGen.setAnacondaAndEnvironment("nwchempy")
        self.assertEqual(actual_line, expected_line)

    def test_write_module_load(self):
        """test to see expected modules are loaded. 
        """
        expected_line = """
module add gcc/4.8.3
module add mpi/openmpi/gnu/1.8.1
module add nwchem/6.7_mw529_anaconda_2.4.1
"""
        actual_line = ZiggyGen.setModulesToLoad()
        self.assertEqual(actual_line, expected_line)

    def test_set_paths(self):
        """test to see if paths are correctly set
        """
        expected_line = """
export SCRATCH_DIR=/scratch/\${LOGNAME}

# Change to directory where you submitted the job from
cd \${PBS_O_WORKDIR}
export PYTHONPATH=\${PYTHONPATH}:\${PBS_O_WORKDIR}
"""
        actual_line = ZiggyGen.setPaths()
        self.assertEqual(actual_line, expected_line)

    def test_add_cml_util(self):
        """Test to see if correct line is added to append nwchemcmlutils to
        pythonpath.
        """
        expected_line = """
#add nwchemcmlutils directory to python path- crude way
export PYTHONPATH=\${PYTHONPATH}:/home/mdd31/nwchemcmlutils
"""
        actual_line = ZiggyGen.addCMLUtilsToPythonPath()
        self.assertEqual(actual_line, expected_line)

    def test_write_calculation_line(self):
        """Test to see if expected calculation line is written.
        """
        expected_line = """
mpiexec -np 16 nwchem expected_input.nwin &> expected_input.nwlog && wait $!
"""
        actual_line = ZiggyGen.writeCalculationLine("expected_input", 16)
        self.assertEqual(actual_line, expected_line)

    def test_write_calc_clean_up_lines(self):
        """test to see if expected clean up lines are written.
        """
        expected_line = """
# Clean ups
rm -f nwchem.py*
rm -f molecule*.hess
rm -f molecule.movecs
rm -f molecule.db
"""
        actual_line = ZiggyGen.writeCalcCleanUp("molecule")
        self.assertEqual(actual_line, expected_line)

    def test_write_cube_clean_up(self):
        """Test to see if expected lines are written for cube file clean up.
        """
        expected_line = """
#Look to see if merged cube files exist.
export merged_files=`ls molecule_*_merged.cube`


echo "merged_files is:"
echo \${merged_files}

for merged_file in \${merged_files}; do
    #for each merged cube file in the list check they are a normal file.
    echo \${merged_file}
    if [ -f \${merged_file} ]
        then
        #then if it exists, remove all the intermediate files for that surface
        #this uses search and replacement of merged.cube
        rm *_*_*_*.cube
    fi
done
"""
        actual_line = ZiggyGen.writeCubeFileCleanUp("molecule")
        self.assertMultiLineEqual(actual_line, expected_line)

    def test_write_clean_up(self):
        """Test to see if expected clean up section is written for the 
        """
        expected_clean_up_1 = """
# Clean ups
rm -f nwchem.py*
rm -f molecule*.hess
rm -f molecule.movecs
rm -f molecule.db
"""
        actual_clean_up_1 = ZiggyGen.writeCleanUp("molecule", False)
        self.assertMultiLineEqual(actual_clean_up_1, expected_clean_up_1)
        expected_clean_up_2 = """
# Clean ups
rm -f nwchem.py*
rm -f molecule*.hess
rm -f molecule.movecs
rm -f molecule.db

#Look to see if merged cube files exist.
export merged_files=`ls molecule_*_merged.cube`


echo "merged_files is:"
echo \${merged_files}

for merged_file in \${merged_files}; do
    #for each merged cube file in the list check they are a normal file.
    echo \${merged_file}
    if [ -f \${merged_file} ]
        then
        #then if it exists, remove all the intermediate files for that surface
        #this uses search and replacement of merged.cube
        rm *_*_*_*.cube
    fi
done
"""
        actual_clean_up_2 = ZiggyGen.writeCleanUp("molecule", True)
        self.assertMultiLineEqual(actual_clean_up_2, expected_clean_up_2)

    def test_write_end_of_file(self):
        """test to see that expected end of file is written.
        """
        expected_line = "EOF\n"
        actual_line = ZiggyGen.writeFileEnd()
        self.assertEqual(actual_line, expected_line)

    def test_write_submit_file_contents(self):
        """test to see if expected file contents are returned.
        """
        self.maxDiff = None
        expected_line_1 = """#!/bin/bash -x
# Start of submission script


/opt/torque/bin/qsub << EOF

########################
# Start PBS directives #
########################


#          Set the name of the job (up to 15 characters,
#          no blank spaces, start with alphanumeric character)

#PBS -N default_name
#          Specify the queue.

#PBS -q s16

#          Specify the maximum wall clock time.
#          Specify the number of nodes requested and the
#          number of processors per node.


#PBS -l nodes=1:ppn=16,walltime=24:00:00

#          The directive below directs that the standard output and
#          error streams are to be merged, intermixed, as standard
#          output.

#PBS -j oe

########################
# End PBS directives   #
########################

module purge

module load anaconda/python2/2.4.1
source activate nwchempy


module add gcc/4.8.3
module add mpi/openmpi/gnu/1.8.1
module add nwchem/6.7_mw529_anaconda_2.4.1

export SCRATCH_DIR=/scratch/\${LOGNAME}

# Change to directory where you submitted the job from
cd \${PBS_O_WORKDIR}
export PYTHONPATH=\${PYTHONPATH}:\${PBS_O_WORKDIR}

#add nwchemcmlutils directory to python path- crude way
export PYTHONPATH=\${PYTHONPATH}:/home/mdd31/nwchemcmlutils

mpiexec -np 16 nwchem expected_input.nwin &> expected_input.nwlog && wait $!

# Clean ups
rm -f nwchem.py*
rm -f molecule*.hess
rm -f molecule.movecs
rm -f molecule.db
EOF
"""
        actual_line_1 = ZiggyGen.writeSubmitFileContents("expected_input", "molecule")
        self.assertMultiLineEqual(actual_line_1, expected_line_1)
        expected_line_2 = """#!/bin/bash -x
# Start of submission script


/opt/torque/bin/qsub << EOF

########################
# Start PBS directives #
########################


#          Set the name of the job (up to 15 characters,
#          no blank spaces, start with alphanumeric character)

#PBS -N default_name
#          Specify the queue.

#PBS -q s16

#          Specify the maximum wall clock time.
#          Specify the number of nodes requested and the
#          number of processors per node.


#PBS -l nodes=1:ppn=16,walltime=24:00:00

#          The directive below directs that the standard output and
#          error streams are to be merged, intermixed, as standard
#          output.

#PBS -j oe

########################
# End PBS directives   #
########################

module purge

module load anaconda/python2/2.4.1
source activate nwchempy


module add gcc/4.8.3
module add mpi/openmpi/gnu/1.8.1
module add nwchem/6.7_mw529_anaconda_2.4.1

export SCRATCH_DIR=/scratch/\${LOGNAME}

# Change to directory where you submitted the job from
cd \${PBS_O_WORKDIR}
export PYTHONPATH=\${PYTHONPATH}:\${PBS_O_WORKDIR}

#add nwchemcmlutils directory to python path- crude way
export PYTHONPATH=\${PYTHONPATH}:/home/mdd31/nwchemcmlutils

mpiexec -np 16 nwchem expected_input.nwin &> expected_input.nwlog && wait $!

# Clean ups
rm -f nwchem.py*
rm -f molecule*.hess
rm -f molecule.movecs
rm -f molecule.db

#Look to see if merged cube files exist.
export merged_files=`ls molecule_*_merged.cube`


echo "merged_files is:"
echo \${merged_files}

for merged_file in \${merged_files}; do
    #for each merged cube file in the list check they are a normal file.
    echo \${merged_file}
    if [ -f \${merged_file} ]
        then
        #then if it exists, remove all the intermediate files for that surface
        #this uses search and replacement of merged.cube
        rm *_*_*_*.cube
    fi
done
EOF
"""
        actual_line_2 = ZiggyGen.writeSubmitFileContents(
            "expected_input", "molecule", cube_clean=True
        )
        self.assertMultiLineEqual(actual_line_2, expected_line_2)

    def test_write_submit_file_contents_to_file(self):
        """test to see if expected file is written out.
        """
        expected_filename = (
            (self.parent_directory / "resources/expected_ziggy_submit.sh")
            .absolute()
            .as_posix()
        )
        actual_filename = "actual_ziggy_submit.sh"
        submit_file_contents = ZiggyGen.writeSubmitFileContents(
            "expected_input", "molecule"
        )
        ZiggyGen.writeSubmitFileContentsToFile(actual_filename, submit_file_contents)
        with open(expected_filename, "r") as expected_file:
            with open(actual_filename, "r") as actual_file:
                actual_file_read_in = actual_file.read()
                expected_file_read_in = expected_file.read()
                self.assertMultiLineEqual(actual_file_read_in, expected_file_read_in)

    def test_write_submit_file(self):
        """Test to see if expected file is created.
        """
        expected_filename = (
            (self.parent_directory / "resources/expected_input2_submit.sh")
            .absolute()
            .as_posix()
        )
        actual_filename = "actual_input2_submit.sh"
        submit_file = ZiggyGen.writeSubmitFile(
            "actual_input2", "molecule", directory=""
        )
        with open(expected_filename, "r") as expected_file:
            with open(actual_filename, "r") as actual_file:
                actual_file_read_in = actual_file.read()
                expected_file_read_in = expected_file.read()
                self.assertMultiLineEqual(actual_file_read_in, expected_file_read_in)

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
script containing the tests for the slurmsubmitgenerator.

@author: mark
"""

import unittest
import logging
import os
import pathlib
import nwchemcmlutils.nwinFileGen.slurmsubmitgenerator as slurmgen

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.WARN)


class SLURMSubmitGeneratorTestCase(unittest.TestCase):
    """Class for the unit tests for the SLURM submit script.
    """

    def setUp(self):
        """set up before tests.
        """
        self.maxDiff = None
        self.parent_directory = pathlib.Path(__file__).parents[1]

    def tearDown(self):
        """clean up after tests.
        """
        if os.path.isfile("actual_slurm_submit.sh"):
            os.remove("actual_slurm_submit.sh")
        if os.path.isfile("actual_input2_submit.sh"):
            os.remove("actual_input2_submit.sh")

    def test_set_slurm_job_name(self):
        """Test to see if job name is set as expected.
        """
        expected_line = """#!/bin/bash
# Name of the job
#SBATCH -J default_name
"""
        actual_line = slurmgen.set_slurm_job_name("default_name")
        self.assertEqual(expected_line, actual_line)

    def test_set_account(self):
        """Test to see if expected account line is written.
        """
        expected_line = "\n# account to charge\n#SBATCH -A HUNTER-SL2\n"
        actual_line = slurmgen.set_account("HUNTER-SL2")
        self.assertEqual(expected_line, actual_line)

    def test_set_partition(self):
        """Test to see if partition is correctly set.
        """
        expected_line = """
# Partition to use
#SBATCH -p HUNTER
"""
        actual_line = slurmgen.set_partition("HUNTER")
        self.assertEqual(expected_line, actual_line)

    def test_set_walltime(self):
        """Test to see if walltime is set correctly.
        """
        expected_line = """
# Time limit. Often not needed as there are defaults, 
# but you might have to specify it to get the maximum allowed
# time
#SBATCH --time=24:00:00
"""
        actual_line = slurmgen.set_walltime("24:00:00")
        self.assertEqual(expected_line, actual_line)

    def test_set_number_of_processors(self):
        """Test to see if number of processors is set as expected.
        """
        expected_line = """
# Number of processors
#SBATCH -n16
"""
        actual_line = slurmgen.set_number_of_processors(16)
        self.assertEqual(expected_line, actual_line)

    def test_set_number_of_nodes(self):
        """Test to see if expected lines are written.
        """
        expected_line = """# Number of processors
#SBATCH -N1
"""
        actual_line = slurmgen.set_number_of_nodes(1)
        self.assertEqual(expected_line, actual_line)

    def test_write_slurm_section(self):
        """Test to see if SLURM section is as expected.
        """
        expected_line = """#!/bin/bash
# Name of the job
#SBATCH -J default_name

# Partition to use
#SBATCH -p HUNTER

# Time limit. Often not needed as there are defaults, 
# but you might have to specify it to get the maximum allowed
# time
#SBATCH --time=24:00:00

# Number of processors
#SBATCH -n16
# Number of processors
#SBATCH -N1
"""
        actual_line = slurmgen.write_slurm_section()
        self.assertEqual(expected_line, actual_line)

    def test_writeModuleSetUp(self):
        """Test to see if expected lines are written
        """
        expected_line = "\n. /etc/profile.d/modules.sh\n"
        actual_line = slurmgen.writeModuleSetUp()
        self.assertEqual(expected_line, actual_line)

    def test_write_hpc_section(self):
        """Test to check expected section was written.
        """
        expected_section = """#!/bin/bash
# Name of the job
#SBATCH -J default_name

# Partition to use
#SBATCH -p HUNTER

# Time limit. Often not needed as there are defaults, 
# but you might have to specify it to get the maximum allowed
# time
#SBATCH --time=24:00:00

# Number of processors
#SBATCH -n16
# Number of processors
#SBATCH -N1

. /etc/profile.d/modules.sh

module purge
module add rhel7/default-peta4
module add gcc-5.4.0-gcc-4.8.5-fis24gg
module add intel/impi/2017.4/gnu
module add miniconda2-4.3.14-gcc-5.4.0-xjtq53h
source activate nwchempy


dir=/tmp/$SLURM_JOB_ID; mkdir -p $dir; trap "rm -r $dir" EXIT

ln -snf $dir ./scratch
export NWCHEM_TOP=/home/mdd31/rds/hpc-work/nwchem-6.8.1-release
export PATH=$PATH:$NWCHEM_TOP/bin/LINUX64
export PYTHONHOME=~/.conda/envs/nwchempy
export LD_LIBRARY_PATH=~/.conda/envs/nwchempy/lib:$LD_LIBRARY_PATH

#add nwchemcmlutils directory to python path- crude way
export PYTHONPATH=$PYTHONPATH:$NWCHEM_TOP/contrib/python:${CONDA_PREFIX}/lib/python3.11/site-packages


#! Full path to application executable:
application="nwchem"

#! Run options for the application:
options="molecule.nwin &> molecule.nwlog"

#! Work directory (i.e. where the job will run):
workdir="$SLURM_SUBMIT_DIR"  # The value of SLURM_SUBMIT_DIR sets workdir to the directory
                             # in which sbatch is run.

numnodes=$SLURM_JOB_NUM_NODES
mpi_tasks_per_node=$SLURM_CPUS_ON_NODE

#! Are you using OpenMP (NB this is unrelated to OpenMPI)? If so increase this
#! safe value to no more than 32:
export OMP_NUM_THREADS=1

#! Number of MPI tasks to be started by the application per node and in total (do not change):
np=$(($numnodes * $mpi_tasks_per_node))

#! The following variables define a sensible pinning strategy for Intel MPI tasks -
#! this should be suitable for both pure MPI and hybrid MPI/OpenMP jobs:
export I_MPI_PIN_DOMAIN=omp:compact # Domains are $OMP_NUM_THREADS cores in size
export I_MPI_PIN_ORDER=scatter # Adjacent domains have minimal sharing of caches/sockets
#! Notes:
#! 1. These variables influence Intel MPI only.
#! 2. Domains are non-overlapping sets of cores which map 1-1 to MPI tasks.
#! 3. I_MPI_PIN_PROCESSOR_LIST is ignored if I_MPI_PIN_DOMAIN is set.
#! 4. If MPI tasks perform better when sharing caches/sockets, try I_MPI_PIN_ORDER=compact.

#! Uncomment one choice for CMD below (add mpirun/mpiexec options if necessary):

#! Choose this for a MPI code (possibly using OpenMP) using Intel MPI.
CMD="mpirun  -np $np $application $options"

#! Choose this for a pure shared-memory OpenMP parallel program on a single node:
#! (OMP_NUM_THREADS threads will be created):
#CMD="$application $options"

#! Choose this for a MPI code (possibly using OpenMP) using OpenMPI:
#CMD="mpirun -npernode $mpi_tasks_per_node -np $np $application $options"

###############################################################
### You should not have to change anything below this line ####
###############################################################

cd $workdir
echo -e "Changed directory to `pwd`.
"
JOBID=$SLURM_JOB_ID
echo -e "JobID: $JOBID
======"
echo "Time: `date`"
echo "Running on master node: `hostname`"
echo "Current directory: `pwd`"

if [ "$SLURM_JOB_NODELIST" ]; then
        #! Create a machine file:
        touch machine.file.$JOBID
        echo -e "
Nodes allocated:
================"
        echo `cat machine.file.$JOBID | sed -e 's/\..*$//g'`
fi
echo -e "
numtasks=$numtasks, numnodes=$numnodes, mpi_tasks_per_node=$mpi_tasks_per_node (OMP_NUM_THREADS=$OMP_NUM_THREADS)"
echo -e "
Executing command:
==================
$CMD
"
eval $CMD

"""
        actual_section = slurmgen.write_hpc_section("molecule", "molecule")
        self.maxDiff = None
        self.assertMultiLineEqual(expected_section, actual_section)

    def test_write_module_purge(self):
        """test to check that correct line is returned.
        """
        expected_line = "\nmodule purge\n"
        actual_line = slurmgen.writeModulePurge()
        self.assertEqual(expected_line, actual_line)

    def test_anaconda_set_up(self):
        """test to see expected anaconda version and environment is set.
        """
        expected_line = """
source activate nwchempy

"""
        actual_line = slurmgen.set_anaconda_environment("nwchempy")
        self.assertEqual(expected_line, actual_line)

    def test_set_modules(self):
        """test to see expected modules are loaded. 
        """
        expected_line = """
. /etc/profile.d/modules.sh

module purge
module add gcc/4.8.3
module add mpi/openmpi/gnu/2.1.2
module add nwchem/6.9_febae52_anaconda_4.4.0
module add anaconda/python3/4.4.0"""
        actual_line = slurmgen.set_modules(
            ["gcc/4.8.3", "mpi/openmpi/gnu/2.1.2", "nwchem/6.9_febae52_anaconda_4.4.0"],
            "anaconda/python3/4.4.0"
        )
        self.assertMultiLineEqual(expected_line, actual_line)

    def test_add_cml_util(self):
        """Test to see if correct line is added to append nwchemcmlutils to
        pythonpath.
        """
        expected_line = """
#add nwchemcmlutils directory to python path- crude way
export PYTHONPATH=${PYTHONPATH}:${CONDA_PREFIX}/lib/python3.11/site-packages
"""
        actual_line = slurmgen.addCMLUtilsToPythonPath("PYTHONPATH=${PYTHONPATH}:${CONDA_PREFIX}/lib/python3.11/site-packages")
        self.assertEqual(actual_line, expected_line)

    def test_write_calculation_line(self):
        """Test to see if expected calculation line is written.
        """
        expected_line = """

#! Full path to application executable:
application="nwchem"

#! Run options for the application:
options="expected_input.nwin &> expected_input.nwlog"

#! Work directory (i.e. where the job will run):
workdir="$SLURM_SUBMIT_DIR"  # The value of SLURM_SUBMIT_DIR sets workdir to the directory
                             # in which sbatch is run.

numnodes=$SLURM_JOB_NUM_NODES
mpi_tasks_per_node=$SLURM_CPUS_ON_NODE

#! Are you using OpenMP (NB this is unrelated to OpenMPI)? If so increase this
#! safe value to no more than 32:
export OMP_NUM_THREADS=1

#! Number of MPI tasks to be started by the application per node and in total (do not change):
np=$(($numnodes * $mpi_tasks_per_node))

#! The following variables define a sensible pinning strategy for Intel MPI tasks -
#! this should be suitable for both pure MPI and hybrid MPI/OpenMP jobs:
export I_MPI_PIN_DOMAIN=omp:compact # Domains are $OMP_NUM_THREADS cores in size
export I_MPI_PIN_ORDER=scatter # Adjacent domains have minimal sharing of caches/sockets
#! Notes:
#! 1. These variables influence Intel MPI only.
#! 2. Domains are non-overlapping sets of cores which map 1-1 to MPI tasks.
#! 3. I_MPI_PIN_PROCESSOR_LIST is ignored if I_MPI_PIN_DOMAIN is set.
#! 4. If MPI tasks perform better when sharing caches/sockets, try I_MPI_PIN_ORDER=compact.

#! Uncomment one choice for CMD below (add mpirun/mpiexec options if necessary):

#! Choose this for a MPI code (possibly using OpenMP) using Intel MPI.
CMD="mpirun -ppn $mpi_tasks_per_node -np $np $application $options"

#! Choose this for a pure shared-memory OpenMP parallel program on a single node:
#! (OMP_NUM_THREADS threads will be created):
#CMD="$application $options"

#! Choose this for a MPI code (possibly using OpenMP) using OpenMPI:
#CMD="mpirun -npernode $mpi_tasks_per_node -np $np $application $options"

###############################################################
### You should not have to change anything below this line ####
###############################################################

cd $workdir
echo -e "Changed directory to `pwd`.\n"
JOBID=$SLURM_JOB_ID
echo -e "JobID: $JOBID\n======"
echo "Time: `date`"
echo "Running on master node: `hostname`"
echo "Current directory: `pwd`"

if [ "$SLURM_JOB_NODELIST" ]; then
        #! Create a machine file:
        export NODEFILE=`generate_pbs_nodefile`
        cat $NODEFILE | uniq > machine.file.$JOBID
        echo -e "\nNodes allocated:\n================"
        echo `cat machine.file.$JOBID | sed -e 's/\..*$//g'`
fi
echo -e "\nnumtasks=$numtasks, numnodes=$numnodes, mpi_tasks_per_node=$mpi_tasks_per_node (OMP_NUM_THREADS=$OMP_NUM_THREADS)"
echo -e "\nExecuting command:\n==================\n$CMD\n"
eval $CMD

"""
        actual_line = slurmgen.write_calculation_section("expected_input",
                                                         npernodesec="-ppn $mpi_tasks_per_node",
        machinefile="""export NODEFILE=`generate_pbs_nodefile`
        cat $NODEFILE | uniq > machine.file.$JOBID""")
        self.assertEqual(expected_line, actual_line)

    def test_write_submit_file_contents(self):
        """test to see if expected file contents are returned.
        """
        self.maxDiff = None
        expected_line_1 = """#!/bin/bash
# Name of the job
#SBATCH -J default_name

# Partition to use
#SBATCH -p HUNTER

# Time limit. Often not needed as there are defaults, 
# but you might have to specify it to get the maximum allowed
# time
#SBATCH --time=24:00:00

# Number of processors
#SBATCH -n16
# Number of processors
#SBATCH -N1

. /etc/profile.d/modules.sh

module purge
module add gcc/4.8.3
module add mpi/openmpi/gnu/2.1.2
module add nwchem/6.9_febae52_anaconda_4.4.0
module add anaconda/python3/4.4.0
source activate nwchempy


dir=/tmp/$SLURM_JOB_ID; mkdir -p $dir; trap "rm -r $dir" EXIT

ln -snf $dir ./scratch

#add nwchemcmlutils directory to python path- crude way
export PYTHONPATH=${PYTHONPATH}:${CONDA_PREFIX}/lib/python3.11/site-packages


#! Full path to application executable:
application="nwchem"

#! Run options for the application:
options="expected_input.nwin &> expected_input.nwlog"

#! Work directory (i.e. where the job will run):
workdir="$SLURM_SUBMIT_DIR"  # The value of SLURM_SUBMIT_DIR sets workdir to the directory
                             # in which sbatch is run.

numnodes=$SLURM_JOB_NUM_NODES
mpi_tasks_per_node=$SLURM_CPUS_ON_NODE

#! Are you using OpenMP (NB this is unrelated to OpenMPI)? If so increase this
#! safe value to no more than 32:
export OMP_NUM_THREADS=1

#! Number of MPI tasks to be started by the application per node and in total (do not change):
np=$(($numnodes * $mpi_tasks_per_node))

#! The following variables define a sensible pinning strategy for Intel MPI tasks -
#! this should be suitable for both pure MPI and hybrid MPI/OpenMP jobs:
export I_MPI_PIN_DOMAIN=omp:compact # Domains are $OMP_NUM_THREADS cores in size
export I_MPI_PIN_ORDER=scatter # Adjacent domains have minimal sharing of caches/sockets
#! Notes:
#! 1. These variables influence Intel MPI only.
#! 2. Domains are non-overlapping sets of cores which map 1-1 to MPI tasks.
#! 3. I_MPI_PIN_PROCESSOR_LIST is ignored if I_MPI_PIN_DOMAIN is set.
#! 4. If MPI tasks perform better when sharing caches/sockets, try I_MPI_PIN_ORDER=compact.

#! Uncomment one choice for CMD below (add mpirun/mpiexec options if necessary):

#! Choose this for a MPI code (possibly using OpenMP) using Intel MPI.
CMD="mpirun  -np $np $application $options"

#! Choose this for a pure shared-memory OpenMP parallel program on a single node:
#! (OMP_NUM_THREADS threads will be created):
#CMD="$application $options"

#! Choose this for a MPI code (possibly using OpenMP) using OpenMPI:
#CMD="mpirun -npernode $mpi_tasks_per_node -np $np $application $options"

###############################################################
### You should not have to change anything below this line ####
###############################################################

cd $workdir
echo -e "Changed directory to `pwd`.\n"
JOBID=$SLURM_JOB_ID
echo -e "JobID: $JOBID\n======"
echo "Time: `date`"
echo "Running on master node: `hostname`"
echo "Current directory: `pwd`"

if [ "$SLURM_JOB_NODELIST" ]; then
        #! Create a machine file:
        touch machine.file.$JOBID
        echo -e "\nNodes allocated:\n================"
        echo `cat machine.file.$JOBID | sed -e 's/\..*$//g'`
fi
echo -e "\nnumtasks=$numtasks, numnodes=$numnodes, mpi_tasks_per_node=$mpi_tasks_per_node (OMP_NUM_THREADS=$OMP_NUM_THREADS)"
echo -e "\nExecuting command:\n==================\n$CMD\n"
eval $CMD

"""
        actual_line_1 = slurmgen.writeSubmitFileContents("expected_input", "molecule")
        self.assertMultiLineEqual(actual_line_1, expected_line_1)
        expected_line_2 = """#!/bin/bash
# Name of the job
#SBATCH -J default_name

# Partition to use
#SBATCH -p HUNTER

# Time limit. Often not needed as there are defaults, 
# but you might have to specify it to get the maximum allowed
# time
#SBATCH --time=24:00:00

# Number of processors
#SBATCH -n16
# Number of processors
#SBATCH -N1

. /etc/profile.d/modules.sh

module purge
module add gcc/4.8.3
module add mpi/openmpi/gnu/2.1.2
module add nwchem/6.9_febae52_anaconda_4.4.0
module add anaconda/python3/4.4.0
source activate nwchempy


dir=/tmp/$SLURM_JOB_ID; mkdir -p $dir; trap "rm -r $dir" EXIT

ln -snf $dir ./scratch

#add nwchemcmlutils directory to python path- crude way
export PYTHONPATH=${PYTHONPATH}:${CONDA_PREFIX}/lib/python3.11/site-packages


#! Full path to application executable:
application="nwchem"

#! Run options for the application:
options="expected_input.nwin &> expected_input.nwlog"

#! Work directory (i.e. where the job will run):
workdir="$SLURM_SUBMIT_DIR"  # The value of SLURM_SUBMIT_DIR sets workdir to the directory
                             # in which sbatch is run.

numnodes=$SLURM_JOB_NUM_NODES
mpi_tasks_per_node=$SLURM_CPUS_ON_NODE

#! Are you using OpenMP (NB this is unrelated to OpenMPI)? If so increase this
#! safe value to no more than 32:
export OMP_NUM_THREADS=1

#! Number of MPI tasks to be started by the application per node and in total (do not change):
np=$(($numnodes * $mpi_tasks_per_node))

#! The following variables define a sensible pinning strategy for Intel MPI tasks -
#! this should be suitable for both pure MPI and hybrid MPI/OpenMP jobs:
export I_MPI_PIN_DOMAIN=omp:compact # Domains are $OMP_NUM_THREADS cores in size
export I_MPI_PIN_ORDER=scatter # Adjacent domains have minimal sharing of caches/sockets
#! Notes:
#! 1. These variables influence Intel MPI only.
#! 2. Domains are non-overlapping sets of cores which map 1-1 to MPI tasks.
#! 3. I_MPI_PIN_PROCESSOR_LIST is ignored if I_MPI_PIN_DOMAIN is set.
#! 4. If MPI tasks perform better when sharing caches/sockets, try I_MPI_PIN_ORDER=compact.

#! Uncomment one choice for CMD below (add mpirun/mpiexec options if necessary):

#! Choose this for a MPI code (possibly using OpenMP) using Intel MPI.
CMD="mpirun  -np $np $application $options"

#! Choose this for a pure shared-memory OpenMP parallel program on a single node:
#! (OMP_NUM_THREADS threads will be created):
#CMD="$application $options"

#! Choose this for a MPI code (possibly using OpenMP) using OpenMPI:
#CMD="mpirun -npernode $mpi_tasks_per_node -np $np $application $options"

###############################################################
### You should not have to change anything below this line ####
###############################################################

cd $workdir
echo -e "Changed directory to `pwd`.\n"
JOBID=$SLURM_JOB_ID
echo -e "JobID: $JOBID\n======"
echo "Time: `date`"
echo "Running on master node: `hostname`"
echo "Current directory: `pwd`"

if [ "$SLURM_JOB_NODELIST" ]; then
        #! Create a machine file:
        touch machine.file.$JOBID
        echo -e "\nNodes allocated:\n================"
        echo `cat machine.file.$JOBID | sed -e 's/\..*$//g'`
fi
echo -e "\nnumtasks=$numtasks, numnodes=$numnodes, mpi_tasks_per_node=$mpi_tasks_per_node (OMP_NUM_THREADS=$OMP_NUM_THREADS)"
echo -e "\nExecuting command:\n==================\n$CMD\n"
eval $CMD

"""
        actual_line_2 = slurmgen.writeSubmitFileContents(
            "expected_input", "molecule", cube_clean=True
        )
        self.assertMultiLineEqual(actual_line_2, expected_line_2)

    def test_write_submit_file_contents_to_file(self):
        """test to see if expected file is written out.
        """
        expected_filename = (
            (self.parent_directory / "resources/expected_slurm_submit.sh")
            .absolute()
            .as_posix()
        )
        actual_filename = "actual_slurm_submit.sh"
        submit_file_contents = slurmgen.writeSubmitFileContents(
            "expected_input", "molecule"
        )
        slurmgen.writeSubmitFileContentsToFile(actual_filename, submit_file_contents)
        with open(expected_filename, "r") as expected_file:
            with open(actual_filename, "r") as actual_file:
                actual_file_read_in = actual_file.read()
                expected_file_read_in = expected_file.read()
                self.assertMultiLineEqual(actual_file_read_in, expected_file_read_in)

    def test_write_submit_file(self):
        """Test to see if expected file is created.
        """
        expected_filename = (
            (self.parent_directory / "resources/expected_slurm2_submit.sh")
            .absolute()
            .as_posix()
        )
        actual_filename = "actual_input2_submit.sh"
        submit_file = slurmgen.writeSubmitFile(
            "actual_input2", "molecule", directory=""
        )
        with open(expected_filename, "r") as expected_file:
            with open(actual_filename, "r") as actual_file:
                actual_file_read_in = actual_file.read()
                expected_file_read_in = expected_file.read()
                self.assertMultiLineEqual(actual_file_read_in, expected_file_read_in)

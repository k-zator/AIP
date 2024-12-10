#!/bin/bash
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

#! Optionally modify the environment seen by the application
#! (note that SLURM reproduces the environment at submission irrespective of ~/.bashrc):
. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load rhel7/default-peta4            # REQUIRED - loads the basic environment
module add gcc-5.4.0-gcc-4.8.5-fis24gg
module add intel/impi/2017.4/gnu
module add miniconda2-4.3.14-gcc-5.4.0-xjtq53h
source activate nwchempy
export NWCHEM_TOP=/home/mdd31/rds/hpc-work/nwchem-6.8.1-release
export PATH=$PATH:$NWCHEM_TOP/bin/LINUX64
export PYTHONPATH=$PYTHONPATH:$NWCHEM_TOP/contrib/python:${CONDA_PREFIX}/lib/python2.7/site-packages
export PYTHONHOME=~/.conda/envs/nwchempy
export LD_LIBRARY_PATH=~/.conda/envs/nwchempy/lib:$LD_LIBRARY_PATH
#! Insert additional module load commands after this line if needed:

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
        export NODEFILE=`generate_pbs_nodefile`
        cat $NODEFILE | uniq > machine.file.$JOBID
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


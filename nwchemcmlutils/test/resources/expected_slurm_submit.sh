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


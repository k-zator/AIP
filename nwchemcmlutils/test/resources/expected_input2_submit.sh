#!/bin/bash -x
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

mpiexec -np 16 nwchem actual_input2.nwin &> actual_input2.nwlog && wait $!

# Clean ups
rm -f nwchem.py*
rm -f molecule*.hess
rm -f molecule.movecs
rm -f molecule.db
EOF

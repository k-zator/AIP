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
Script for generating input files for SLURM queuing systems.
"""

import logging

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.WARN)


def set_slurm_job_name(job_name):
    """Set the job name for the slurm script.

    Parameters
    ----------
    job_name : TYPE
        DESCRIPTION.

    Returns
    -------
    str
        job name section.

    """

    slurm_job_name = """#!/bin/bash
# Name of the job
#SBATCH -J {job_name}
""".format(
        job_name=job_name
    )
    return slurm_job_name


def set_account(account):
    """This sets the account to charge, if applicable for resource.
    """
    if account != None:
        return """
# account to charge
#SBATCH -A {account}
""".format(
            account=account
        )
    else:
        return ""


def set_partition(partition):
    """This sets the partition to submit the job to.
    """
    slurm_partition = """
# Partition to use
#SBATCH -p {partition}
""".format(
        partition=partition
    )
    return slurm_partition


def set_walltime(walltime):
    """This sets the walltime for the calculation.
    """
    slurm_walltime = """
# Time limit. Often not needed as there are defaults, 
# but you might have to specify it to get the maximum allowed
# time
#SBATCH --time={walltime}
""".format(
        walltime=walltime
    )
    return slurm_walltime


def set_number_of_processors(processors):
    """This sets the number of processors.
    """
    slurm_processors = """
# Number of processors
#SBATCH -n{processors}
""".format(
        processors=processors
    )
    return slurm_processors


def set_number_of_nodes(nodes):
    """This sets the number of nodes.
    """
    slurm_nodes = """# Number of processors
#SBATCH -N{nodes}
""".format(
        nodes=nodes
    )
    return slurm_nodes


def write_slurm_section(**kwargs):
    """This writes the SLURM set up section.
    """
    slurm_section = set_slurm_job_name(
        job_name=kwargs.pop("job_name", "default_name"))
    slurm_section += set_account(account=kwargs.pop("account", None))
    slurm_section += set_partition(partition=kwargs.pop("partition", "HUNTER"))
    slurm_section += set_walltime(walltime=kwargs.pop("walltime", "24:00:00"))
    slurm_section += set_number_of_processors(
        processors=kwargs.pop("processors", 16))
    slurm_section += set_number_of_nodes(nodes=kwargs.pop("nodes", 1))
    return slurm_section


def writeModuleSetUp():
    """This writes the line to set up the module command, so that it can be called.
    """
    return "\n. /etc/profile.d/modules.sh\n"


def write_hpc_section(filename_stem, molecule_name, **kwargs):
    """This writes contents for CSD3 HPC file.
    """
    file_contents = write_slurm_section(**kwargs)
    file_contents += set_modules(
        kwargs.pop(
            "modules",
            ["rhel7/default-peta4",
             "gcc-5.4.0-gcc-4.8.5-fis24gg",
             "intel/impi/2017.4/gnu"]
        ),
        conda_version=kwargs.pop(
            "conda_version", "miniconda2-4.3.14-gcc-5.4.0-xjtq53h"),
    )
    file_contents += set_anaconda_environment(
        conda_env=kwargs.pop("conda_env", "nwchempy"),
    )
    file_contents += add_scratchdir()
    file_contents += write_path_updates(kwargs.get("paths", ["NWCHEM_TOP=/home/mdd31/rds/hpc-work/nwchem-6.8.1-release",
                                                             "PATH=$PATH:$NWCHEM_TOP/bin/LINUX64",
                                                             "PYTHONHOME=~/.conda/envs/nwchempy",
                                                             "LD_LIBRARY_PATH=~/.conda/envs/nwchempy/lib:$LD_LIBRARY_PATH"]),
                                        cml_utils_path=kwargs.get("utils_path", "PYTHONPATH=$PYTHONPATH:$NWCHEM_TOP/contrib/python:${CONDA_PREFIX}/lib/python3.11/site-packages"))
    file_contents += write_calculation_section(filename_stem)
    return file_contents


def writeModulePurge():
    """function writes 'module purge' to clean up any pre loaded modules.
    """
    return "\nmodule purge\n"


def set_anaconda_environment(conda_env):
    """function writes lines to load the anaconda module, and then set up the
    the anaconda environment specified.
    """
    conda_lines = """
source activate {0}

""".format(
        conda_env
    )
    return conda_lines


def set_modules(modules, conda_version):
    """Write the modules to set for an NWChem calculation.

    Note- depending on the compliation of NWChem by sysadmins, the anaconda
    module may need to be loaded with the other modules, before NWWChem. In this case,
    include it in the list of modules and set conda_version=None.

    Parameters
    ----------
    modules : list
        list of modules to add.
    conda_version : str
        conda version.

    Returns
    -------
    str
        module block for SLURM submission file.

    """
    module_strings = ["module add {0}".format(module) for module in modules]
    if conda_version is not None:
        module_strings.append("module add {0}".format(conda_version))
    module_start = writeModuleSetUp() + writeModulePurge()
    return module_start + "\n".join(module_strings)


def add_scratchdir():
    """This creates a temporary directory.
    """
    scratch_dir = """
dir=/tmp/$SLURM_JOB_ID; mkdir -p $dir; trap "rm -r $dir" EXIT

ln -snf $dir ./scratch
"""
    return scratch_dir


def addCMLUtilsToPythonPath(cml_utils_path):
    """function adds the nwchemcmlutils directory to the pythonpath.
    """
    add_cml_utils = """
#add nwchemcmlutils directory to python path- crude way
export {0}
""".format(cml_utils_path)
    return add_cml_utils


def write_path_updates(path_list, cml_utils_path):
    """Write section conatining information to update all paths.

    Parameters
    ----------
    path_list : list
        List of path variable to update.
    cml_utils_path: str
        Path update required to append NWChemCMLUtils to pythonpath.

    Returns
    -------
    str
        Block of export statements for paths variable to update.

    """
    path_set_list = ["export {0}".format(path) for path in path_list]
    path_set_list.append(addCMLUtilsToPythonPath(cml_utils_path))
    return "\n".join(path_set_list)


def write_calculation_section(filename_stem, **kwargs):
    """Write section to run calculation and set up output.

    Parameters
    ----------
    filename_stem : str
        filename stem.

    Returns
    -------
    str
        calcuation section for a job.

    """
    npernodesec = kwargs.pop("npernodesec", "")
    # -ppn $mpi_tasks_per_node
    machinefile = kwargs.pop("machinefile", "touch machine.file.$JOBID")
    #        export NODEFILE=`generate_pbs_nodefile`
#        cat $NODEFILE | uniq > machine.file.$JOBID
    calculation_section = """

#! Full path to application executable:
application="nwchem"

#! Run options for the application:
options="{filename_stem}.nwin &> {filename_stem}.nwlog"

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
CMD="mpirun {npernodesec} -np $np $application $options"

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
        {machinefile}
        echo -e "\nNodes allocated:\n================"
        echo `cat machine.file.$JOBID | sed -e 's/\..*$//g'`
fi
echo -e "\nnumtasks=$numtasks, numnodes=$numnodes, mpi_tasks_per_node=$mpi_tasks_per_node (OMP_NUM_THREADS=$OMP_NUM_THREADS)"
echo -e "\nExecuting command:\n==================\n$CMD\n"
eval $CMD

""".format(
        filename_stem=filename_stem,
        npernodesec=npernodesec,
        machinefile=machinefile,

    )
    return calculation_section


def writeSubmitFileContents(filename_stem, molecule_name, **kwargs):
    """Function writes the submit file contents
    """
    file_contents = write_slurm_section(**kwargs)
    file_contents += set_modules(
        kwargs.pop(
            "modules",
            ["gcc/7.5.0",
             "mpi/openmpi/gnu7/4.1.0",
             "nwchem/7.0.2"]
        ),
        conda_version=kwargs.pop(
            "conda_version", "hunter_nwchem_anaconda/python3/4.4.0"),
    )
    file_contents += set_anaconda_environment(
        conda_env=kwargs.pop("conda_env", "nwchempy"),
    )
    file_contents += add_scratchdir()
    file_contents += write_path_updates(kwargs.get("paths", []),
                                        cml_utils_path=kwargs.get("utils_path", "PYTHONPATH=${PYTHONPATH}:${CONDA_PREFIX}/lib/python3.11/site-packages"))
    file_contents += write_calculation_section(
        filename_stem=filename_stem,
        npernodesec=kwargs.pop("npernodesec", ""),
        machinefile=kwargs.pop("machinefile", "touch machine.file.$JOBID")
    )
    return file_contents


def writeSubmitFileContentsToFile(filename, submit_file_contents):
    """Function writes submit filecontents to a file.
    """
    with open(filename, "w") as submit_file:
        submit_file.write(submit_file_contents)
    return 0


def writeSubmitFile(
    filename_stem, molecule_name, directory="", slurm_type="ziggy", **kwargs
):
    """Function writes a submit file.
    """
    file_contents = ""
    if slurm_type == "ziggy":
        file_contents = writeSubmitFileContents(
            filename_stem=filename_stem, molecule_name=molecule_name, **kwargs
        )
    elif slurm_type == "camhpc":
        file_contents = write_hpc_section(
            filename_stem=filename_stem, molecule_name=molecule_name, **kwargs
        )
    else:
        file_contents = writeSubmitFileContents(
            filename_stem=filename_stem, molecule_name=molecule_name, **kwargs
        )
    filename = directory + filename_stem + "_submit.sh"
    file_out = writeSubmitFileContentsToFile(
        filename=filename, submit_file_contents=file_contents
    )
    return file_out

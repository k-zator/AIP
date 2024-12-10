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
script containing functions for generating ziggy submit scripts for nwchem
jobs.

@author: mark
"""

import logging

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.WARN)


def writeStartOfSubmitScript():
    """writes start of submission script.
    """
    LOGGER.debug("Writing start of submission script.")
    submission_script_start = """#!/bin/bash -x
# Start of submission script


/opt/torque/bin/qsub << EOF
"""
    return submission_script_start


def writePBSDirectiveStart():
    """Function writes the start of the PBS directive section
    """
    pbs_start = """
########################
# Start PBS directives #
########################

"""
    return pbs_start


def setPBSJobName(job_name):
    """function sets the name of the job in the PBS script.
    name is up to 15 characters, no blank spaces, starting with alphanumeric
    character.
    """
    pbs_set_job_name = """
#          Set the name of the job (up to 15 characters,
#          no blank spaces, start with alphanumeric character)

#PBS -N %s
""" % (
        job_name
    )
    return pbs_set_job_name


def setPBSQueue(queue):
    """function sets the queue for the job to be submitted to. Default is s16.
    """
    pbs_set_queue = """#          Specify the queue.

#PBS -q %s
""" % (
        queue
    )
    return pbs_set_queue


def setPBSWalltimeNodesProcessors(nodes, processors, walltime):
    """function sets the number of nodes, processor number and walltime.
    """
    pbs_walltime_node_proc = """
#          Specify the maximum wall clock time.
#          Specify the number of nodes requested and the
#          number of processors per node.


#PBS -l nodes=%i:ppn=%i,walltime=%s
""" % (
        nodes,
        processors,
        walltime,
    )
    return pbs_walltime_node_proc


def setPBSJobOutput():
    """function sets the output and error streams to be merged and intermixed
    as standard output.
    """
    pbs_job_output = """
#          The directive below directs that the standard output and
#          error streams are to be merged, intermixed, as standard
#          output.

#PBS -j oe
"""
    return pbs_job_output


def writePBSDirectiveEnd():
    """Function sets the end of the pbs directive
    """
    pbs_directive_end = """
########################
# End PBS directives   #
########################
"""
    return pbs_directive_end


def writePBSSection(**kwargs):
    """Function writes the PBS section of the file.
    """
    pbs_section = writePBSDirectiveStart()
    pbs_section += setPBSJobName(job_name=kwargs.get("job_name",
                                 "default_name"))
    pbs_section += setPBSQueue(queue=kwargs.get("queue", "s16"))
    pbs_section += setPBSWalltimeNodesProcessors(
        nodes=kwargs.get("nodes", 1),
        processors=kwargs.get("processors", 16),
        walltime=kwargs.get("walltime", "24:00:00"),
    )
    pbs_section += setPBSJobOutput()
    pbs_section += writePBSDirectiveEnd()
    return pbs_section


def writeModulePurge():
    """function writes 'module purge' to clean up any pre loaded modules.
    """
    return "\nmodule purge\n"


def setAnacondaAndEnvironment(environment):
    """function writes lines to load the anaconda module, and then set up the
    the anaconda environment specified.
    """
    conda_lines = """
module load hunter_nwchem_anaconda/python3/4.4.0
source activate %s

""" % (
        environment
    )
    return conda_lines


def setModulesToLoad():
    """function writes themodules to load
    """
    modules_to_load = """


module add gcc/7.5.0
module add mpi/openmpi/gnu7/4.1.0
module add hunter_nwchem_anaconda/python3/4.4.0
module add nwchem/7.0.2
"""
    return modules_to_load


def setPaths():
    """function writes the lines setting the paths for the scratch directory,
    change to the working directory, and adding working directory to the
    pythonpath.
    """
    set_paths_lines = """
export SCRATCH_DIR=/scratch/\${LOGNAME}

# Change to directory where you submitted the job from
cd \${PBS_O_WORKDIR}
export PYTHONPATH=\${PYTHONPATH}:\${PBS_O_WORKDIR}
"""
    return set_paths_lines


def addCMLUtilsToPythonPath():
    """function adds the nwchemcmlutils directory to the pythonpath.
    """
    add_cml_utils = """
#add nwchemcmlutils directory to python path- crude way
export PYTHONPATH=\${PYTHONPATH}:/home/workspace/nwchemcmlutils
"""
    return add_cml_utils


def writeCalculationLine(filename_stem, processors):
    """function writes the nwchem execution line. It takes a filename stem and
    substitutes it in.
    """
    calculation_line = """
mpiexec -np %i nwchem %s.nwin &> %s.nwlog && wait $!
""" % (
        processors,
        filename_stem,
        filename_stem,
    )
    return calculation_line


def writeCalcCleanUp(molecule_name):
    """function writes the clean up lines for a molecule.
    """
    clean_up_lines = """
# Clean ups
rm -f nwchem.py*
rm -f %s*.hess
rm -f %s.movecs
rm -f %s.db
""" % (
        molecule_name,
        molecule_name,
        molecule_name,
    )
    return clean_up_lines


def writeCubeFileCleanUp(molecule_name):
    """This writes the bash lines for clean up of the cube files.
    """
    cube_clean_up = """
#Look to see if merged cube files exist.
export merged_files=`ls %s_*_merged.cube`


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
""" % (
        molecule_name
    )
    return cube_clean_up


def writeCleanUp(molecule_name, cube_clean):
    """Function writes the clean up lines to remove files after a job.
    """
    clean_up_lines = writeCalcCleanUp(molecule_name)
    if cube_clean:
        clean_up_lines += writeCubeFileCleanUp(molecule_name)
    return clean_up_lines


def writeFileEnd():
    """function writes end of file statement.
    """
    return "EOF\n"


def writeSubmitFileContents(filename_stem, molecule_name, **kwargs):
    """Function writes the submit file contents
    """
    file_contents = writeStartOfSubmitScript()
    file_contents += writePBSSection(**kwargs)
    file_contents += writeModulePurge()
    file_contents += setAnacondaAndEnvironment(
        environment=kwargs.get("environment", "nwchempy")
    )
    file_contents += setModulesToLoad()
    file_contents += setPaths()
    file_contents += addCMLUtilsToPythonPath()
    file_contents += writeCalculationLine(
        filename_stem=filename_stem, processors=kwargs.get("processors", 16)
    )
    file_contents += writeCleanUp(
        molecule_name=molecule_name, cube_clean=kwargs.get("cube_clean", False)
    )
    file_contents += writeFileEnd()
    return file_contents


def writeSubmitFileContentsToFile(filename, submit_file_contents):
    """Function writes submit filecontents to a file.
    """
    with open(filename, "w") as submit_file:
        submit_file.write(submit_file_contents)
    return 0


def writeSubmitFile(filename_stem, molecule_name, directory="", **kwargs):
    """Function writes a submit file.
    """
    file_contents = writeSubmitFileContents(
        filename_stem=filename_stem, molecule_name=molecule_name, **kwargs
    )
    filename = directory + filename_stem + "_submit.sh"
    file_out = writeSubmitFileContentsToFile(
        filename=filename, submit_file_contents=file_contents
    )
    return file_out

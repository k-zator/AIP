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
Script for creating the set up lines required for an NWChem input file.

@author: mark
"""

import logging

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.WARN)


def writeMemoryLine(memory_limit):
    """function sets the memory limit for the calculation to memory_limit in GB.
    """
    LOGGER.info("Setting memory")
    memory_line = "memory %i gb\n" % (memory_limit)
    return memory_line


def setScratchDir(scratch_dir):
    """sets the scratch directory for the calculation- defaults to '/tmp'
    """
    LOGGER.info("Setting scratch directory")
    scratch_line = "scratch_dir %s\n" % scratch_dir
    return scratch_line


def setStart(molecule_name):
    """Start molecule section. defaults to 'molecule'.
    """
    LOGGER.info("Setting Start")
    start_line = "start %s\n" % (molecule_name)
    return start_line


def setTitle(title):
    """sets the titlefor the calculation- defaults to 'calculation'.
    """
    LOGGER.info("Setting title")
    title_line = 'title "%s"\n' % (title)
    return title_line


def setCharge(charge):
    """sets the charge on the molecule. default is '0'.
    """
    LOGGER.info("Setting charge")
    charge_line = "charge %i\n" % (charge)
    return charge_line


def setDFT(multiplicity, exchange, iterations):
    """function sets up the DFT section of the calculation. Creates a string
    block.
    """
    LOGGER.info("Setting DFT block")
    initial_comment_lines = """# DFT module
# http://www.nwchem-sw.org/index.php/Density_Functional_Theory_for_Molecules
"""
    dft_block = """dft
   # Calculate all integrals "on-the-fly"
   direct
   
   #Exchange functional- B3LYP is default
   # B3LYP Combined Exchange and Correlation Functional
   xc %s

   # Do not print the MO vector coefficients; just too much data.
   noprint "final vectors analysis"

   # Set the energy convergence to be at 1e-08, SCF=Tight in Gaussian
   convergence energy 1e-08

   #now we increase the number of iterations in case there is a bad initial guess
   iterations %i

   # Multiplicity
   mult %i

end
""" % (
        exchange,
        iterations,
        multiplicity,
    )

    dft_section = initial_comment_lines + dft_block
    return dft_section


def setGeometryDriver(maxiter):
    """function sets up the driver block for geometry optimisation calculation.
    """
    LOGGER.info("Setting Driver block")
    driver_block = """# Module for the geometry optimisation
# http://www.nwchem-sw.org/index.php/Geometry_Optimization#Geometry_Optimization_with_DRIVER
driver
  # Maximum number of steps allowed for the geometry optimisation
  maxiter %i
end
""" % (
        maxiter
    )
    return driver_block


def setPropertyBlockEPSISO(padding, step_size, iso_surf, tol):
    """Function sets up the property block for calculating the electrostatic
    potential surface (EPS) for the specified isosurface density.

    padding refers to the extra space around the molecule in Angstroms? NEED to check
    step_size is distance between points in grid.
    iso_surf is the isosurface density to use.
    tol is the tolerance for points to lie on the iso surface density.
    """
    LOGGER.info("Setting Property block for EPSISO")
    property_block = """
# Property module
# http://www.nwchem-sw.org/index.php/Properties
# New property keyword in use
property
  grid pad {0:f} step {1:f}
  espiso iso {2:f} tol {3:f}
end
""".format(
        padding, step_size, iso_surf, tol
    )
    return property_block


def setBasis(mol_element_set, **element_basis_definitions):
    """Function writes the basis section of an NWChem input-
    requires a set of all the elements in the molecule. The basis dict supplies
    the basis set to use for each atom. If one isn't specified, 6-31G* is used
    as default, except for iodine, where 6-311G** is specified.

    Extra keyword arguements fed to this function are added to the basis dict-
    use the element name as the keyword, and basis set as the value.
    """
    basis_dict = {
        "default": "6-31G*",
        "Br": "6-311G**",
        "I": "6-311G**",
        "Se": "6-311G**",
    }
    for element, basis_set in element_basis_definitions.items():
        basis_dict[element] = basis_set
    basis_block = "basis\n"
    for element in sorted(mol_element_set):
        if element in basis_dict.keys():
            element_basis = "%s library %s\n" % (element, basis_dict[element])
        else:
            element_basis = "%s library %s\n" % (element, basis_dict["default"])
        basis_block += element_basis
    end_line = "end\n"
    basis_block += end_line
    return basis_block


def setGeometry(
    atom_list,
    units="an",
    autoz_option="autoz",
    center_molecule="nocenter",
    symmetry="noautosym",
):
    """sets the geometry of the molecule.
    header- specifies symmetry to use, units to use, whether to
    recenter molecule and whether to use internal z matrix.
    """
    LOGGER.info("Setting geometry")
    geometry_section = ""
    geometry_header_line = "geometry units %s %s %s %s\n" % (
        units,
        autoz_option,
        center_molecule,
        symmetry,
    )
    end_line = "end\n"
    geometry_section += geometry_header_line
    for atom in atom_list:
        atom_line = "    " + atom
        geometry_section += atom_line
    geometry_section += end_line
    return geometry_section


def writeNwinSetUpLines(memory_limit, **kwargs):
    """Function writes the first few lines of the nwin file.
    """
    LOGGER.info("Start writing set up lines for calculation.")
    set_up_lines = ""
    set_up_lines += writeMemoryLine(memory_limit)
    set_up_lines += setScratchDir(scratch_dir=kwargs.get("scratch_dir", "/tmp"))
    set_up_lines += setStart(molecule_name=kwargs.get("molecule_name", "molecule"))
    set_up_lines += setTitle(title=kwargs.get("title", "calculation"))
    set_up_lines += setCharge(charge=kwargs.get("charge", 0))
    return set_up_lines


def writeDFTAndDriverLines(**kwargs):
    """function to write the dft and driver blocks of code.
    """
    LOGGER.info("Reading in kwargs")
    dft_driver_lines = ""
    dft_driver_lines += setDFT(
        multiplicity=kwargs.get("multiplicity", 1),
        exchange=kwargs.get("exchange", "b3lyp"),
        iterations=kwargs.get("iterations", 100),
    )
    dft_driver_lines += setGeometryDriver(maxiter=kwargs.get("maxiter", 600))
    return dft_driver_lines


def setTask(task):
    """sets the task to be run
    """
    LOGGER.info("Setting Task")
    task_line = "task %s\n" % (task)
    return task_line


# New SSIP method
This is a comprehensive algorithm  created by multiple PhD students.

1) NWChemCMLGenerator - MD Driver's code to create CML for DFT calculation. For further details, see the nwchemcmlutils folder

2) AIP_footprinting - KJ Zator and MC Storer's code to create ssip.xml file with AIP description of the molecules 

3) ssip-visualisation - MD Driver's code to view SSIPs/AIPs

4) aip-visualisation - KJ Zator and MC Storer's code to visualise AIPs (not yet)

They will be described in succession.

## Installation

The module can be installed in python 3.x.

#### Ziggy environment setup ####

Note that if you are on ziggy you need to first load anaconda:

        module load anaconda/python3/4.4.0
        git clone git@gitlab.developers.cam.ac.uk:ch/hunter/newssiptools/newssip.git
        cd newssip

        conda env create -f environment36.yml
        source activate newssip36
        pip install .

#### Using pip on any machine ###

Clone the repository:

        git clone git@gitlab.developers.cam.ac.uk:ch/hunter/newssiptools/newssip.git
        cd newssip

        conda env create -f environment.yml
        #enter environment
        source activate newssip
        pip install .

##  NWChemCMLGenerator 

When in the python environment where it is installed, it can be called on the command line using:
    
        python -m cmlgenerator generate -f {mol type} -o {directory}

Otherwise, more aid is available with -h flag. There are two more arguments which will be necessary for proper footprinting:
  -t, --aip_atom_types  Specify if you want to add aip atom types to the cml
  -y, --sybyl           Specify if you want to add sybyl atom types to the cml

This will create a CML file for the next steps.

###  AIP_footprinting 

The algorithm requires completed merged.cube files and appropriately formatted CML file. Terminal command is:

        python -m aip_footprinting -c {}.cml -m {}_0.0300_merged.cube -b {}_0.0104_merged.cube -n {}_0.002_merged.cube -w ssip.xml

  --cml_file       -c  cml,           cml (along with aipAtomType information)
  --cube_polar     -m  cube_polar,    cube at the 0.0300 e.Bohr-3
  --cube_nonpolar  -n  cube_nonpolar, cube at the 0.0020 e.Bohr-3
  --cube_middle    -b  cube_middle,   cube at the 0.0104 e.Bohr-3
  --write          -w  PATH,          path for the output file
  --centre_surface_percentile, -csp  (default: 80) - percentile of surface MEP density for edge cut-off

It produces the ssip.xml file which can be used for further calculations or viewed.

###  SSIP visualisation  
The standard visualisation of the AIPs can be achieved with:

      ssip-vis -r ssip.xml

Here all the AIPs are red or white, negative or positive, transparent spheres with areas proportional to the value.
The code was originally written for SSIPs according to its rather strict Schema, hence sone inconsistencies will be flaged up.

#### Who do I talk to?

Any queries please contact Katarzyna Zator, kz265.

### License

&copy; Katarzyna Zator, Maria Chiara Storer,  Christopher Hunter at the University of Cambridge
This is released under an AGPLv3 license for academic use.



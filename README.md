
# AIP

Module for calculating AIP for neutral organic molecules to be used to examine non-covalent interactions (as AIP-AIP pairing) as well as examinations of solvation and co-crystal prediction.
The relevant paper: https://doi.org/10.1039/D3SC05690B

Author: Katarzyna Zator

#### Environment setup.

Clone the repository:

    git clone https://github.com/k-zator/AIP.git
    cd AIP

Then install this environment instead (based on python 3.6):

        conda env create -f environment.yml
        conda activate aip

#### Calculation.

To obtain an AIP description of a molecule, several files are necessary. Primarily, the custom MEPS cube which can be calculated using the AIP_MEPS script (https://github.com/k-zator/AIP_MEPS). It produces the necessary three isosurface cube files (0.0020, 0.0300, 0.0104) with correct formatting. Secondly, the module requires the input molecule be specified as a CML file to feature correct atom types. This is generated with
    
        python -m cmlgenerator generate -f mol2 -o . -y -t water.mol2

for an example water.mol2 molecule. The flags (-y -t) are necessary to specify the correct atom types.
The created CML file will be saved with the name formatted as an Inchikey but this need not be the case for the AIP calculation following (such was not the case for the SSIP calculations). To obtain the final aip.xml file, run

        python -m aip_footprinting -c water.cml -m water_0.0300.cube -b water_0.0104.cube -n water_0.002.cube -w aip.xml

The aip.xml file can then be used used in AIP_map module (https://github.com/k-zator/AIP_map) to look at NCIs.

###  AIP visualisation.

The standard visualisation of the AIPs can be achieved with:

      jmol ssip-visualisation-7.0.0-jar-with-dependencies -r aip.xml

Here all the AIPs are red or white, negative or positive, transparent spheres with areas proportional to the value.
The code was originally written for SSIPs according to its rather strict Schema, hence some formatting inconsistencies might be flagged up.


Who do I talk to?

Any queries please contact Katarzyna Zator, kz265.

License

&copy; Katarzyna Zator, Maria Chiara Storer,  Christopher Hunter at the University of Cambridge
This is released under an AGPLv3 license for academic use.



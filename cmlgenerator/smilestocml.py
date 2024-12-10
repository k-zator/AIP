# -*- coding: utf-8 -*-
#    cmlgenerator creates fully namespaced CML for molecules from input structures.
#    Copyright (C) 2019  Mark D. Driver
#
#    cmlgenerator is free software: you can redistribute it and/or modify
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
"""
Script for converting SMILES to 3D molecule representation.

:Authors:
    Mark Driver <mdd31>
"""

import logging
from rdkit import Chem
from rdkit.Chem import AllChem
import cmlgenerator.structuretocml as structuretocml

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.WARN)


def generate_cml_and_inchikey_from_smiles(smiles_string, **kwargs):
    """Generate unqualified CML and InChIKey from input SMILES/
    

    Parameters
    ----------
    smiles_string : str
        SMILES string.
    max_confs : int, optional
        Maximum number of conformers to search for. The default is 10.
    MaxIters : int, optional
        Maximum number of iterations in calculation. The default is 1000.


    Returns
    -------
    inchikey : str
        InChIKey.
    cml : TYPE
        CML with unqualified attributes.

    """
    cml = generate_cml_from_smiles(smiles_string, **kwargs)
    inchikey = structuretocml.generate_inchikey_for_cml(cml)
    return inchikey, cml


def generate_cml_from_smiles(smiles_string, **kwargs):
    """Generate CML (attribute unqualified) from SMILES.

    Parameters
    ----------
    smiles_string : str
        SMILES string.
    max_confs : int, optional
        Maximum number of conformers to search for. The default is 10.
    MaxIters : int, optional
        Maximum number of iterations in calculation. The default is 1000.

    Returns
    -------
    str
        unnamespaced CML string.

    """
    rd_mol = generate_mol_from_smiles(smiles_string)
    rd_mol, min_conf_id = generate_3d_structure_for_mol(rd_mol, **kwargs)
    return convert_structure_to_cml(rd_mol, min_conf_id)


def convert_structure_to_cml(rd_mol, conf_id=-1):
    """Convert rdkit.Chem.rdchem.Mol to CML format.
    
    A mol2 version of the molecule is outputted from RDKit and then converted
    to CML with openbabel.

    Parameters
    ----------
    rd_mol : rdkit.Chem.rdchem.Mol 
        Molecule to convert to CML.
    conf_id : int, optional
        conformer ID to return. The default is -1.

    Returns
    -------
    str
        CML moelcule string.

    """
    mol_block = convert_structure_to_molblock(rd_mol, confId=conf_id)
    return structuretocml.perform_conversion_from_string(mol_block, "mol")


def convert_structure_to_molblock(rd_mol, confId=-1):
    """Output mol2 representation of rd_mol molecule.
    
    Conformer selected for output is determined by confId

    Parameters
    ----------
    rd_mol : rdkit.Chem.rdchem.Mol
        RDkit molecule internal representation.
    confId : int, optional
        Conformer to select. The default is -1.

    Returns
    -------
    str
        mol2 representation of molecule.

    """
    return AllChem.MolToMolBlock(rd_mol, confId=confId)


def GenerateAndMinimiseConformers(mol, max_confs=10, MaxIters=1000):
    """Generate conformers returning mol and minimum energy conformer ID.
    
    This uses the ETKDG algorithm for conformer generation and the UFF
    force field for energy evaluation.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        RDkit molecule internal representation.
    max_confs : int, optional
        Maximum number of conformers to search for. The default is 10.
    MaxIters : int, optional
        Maximum number of iterations in calculation. The default is 1000.

    Returns
    -------
    mol : rdkit.Chem.rdchem.Mol 
        Molecule with 3D conformers.
    int :
        ID for minimum energy conformer.

    """
    # Project to 3D
    # https://github.com/rdkit/UGM_2015/blob/master/Presentations/ETKDG.SereinaRiniker.pdf
    # http://dx.doi.org/10.1021/acs.jcim.5b00654
    # http://rdkit.org/docs/GettingStartedInPython.html#working-with-3d-molecules
    if len(mol.GetAtoms()) == 1:
        conf_id = AllChem.EmbedMolecule(mol)
        LOGGER.debug("Single atom found. ConfID:")
        LOGGER.debug(conf_id)
        return mol, conf_id
    else:
        AllChem.EmbedMultipleConfs(mol, max_confs, AllChem.ETKDG())

        # Minimise and sort
        calcEnergies = AllChem.MMFFOptimizeMoleculeConfs(mol)
        sorted_energies = sorted(
            [[calcEnergies[i], i] for i in range(len(calcEnergies))],
            key=lambda tup: tup[0][1],
        )
        min_energy_confid = sorted_energies[0][1]
        LOGGER.debug(
            "MinEnergyConfID: "
            + str(sorted_energies[0][1])
            + " Energy: "
            + str(sorted_energies[0][0][1])
            + " kcal/mol"
        )
        return mol, min_energy_confid


def generate_3d_structure_for_mol(mol, **kwargs):
    """Generate 3D structure for Mol.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol 
        RDkit molecule internal representation.
    

    Returns
    -------
    rdkit.Chem.rdchem.Mol, conf_id
        Mol and conf_id for lowest energy structure.

    """
    Chem.SanitizeMol(mol)
    mol = Chem.AddHs(mol)
    return GenerateAndMinimiseConformers(mol, **kwargs)


def generate_mol_from_smiles(smiles_string):
    """Generate rdkit Mol object from SMILES string.

    Parameters
    ----------
    smiles_string : str
        SMILES.

    Returns
    -------
    rdkit.Chem.rdchem.Mol 
        RDkit molecule internal representation.

    """
    return AllChem.MolFromSmiles(smiles_string)

"""
Parses information from Atom_df and AIP into arrays.
@author: Katarzyna Joanna Zator (kz265)
"""
import mendeleev
import numpy as np
import pandas as pd


class AtomSet():
    """Describes the atoms found in the cube and cml files and presents readable information"""

    def __init__(self):
        self.Atoms_df = None
        self.Atom = []

    def _create_from_MEPS(self, MEPS):
        self.Atoms_df = MEPS.Atoms_df
        self._get_linearised_data()
        # and initialise Atom class objects for each
        self.Atom = []
        for index, row in self.Atoms_df.iterrows():
            atom = Atom()
            atom._create_atom_from_df(index, row)
            self.Atom.append(atom)
        if any("N.pl3" in x for x in self.atom_type):
            N_atoms = np.where(["N.pl3" in x for x in self.atom_type])[0]
            for a in N_atoms:
                N_atom = self.Atom[a]
                neighbours = self.find_connecting_atom(MEPS, N_atom)
                self.check_pyrimidalisation(N_atom, neighbours)

    def _get_linearised_data(self):
        self.indices = self.Atoms_df.index.values
        self.xyz = self.Atoms_df[["x", "y", "z"]].to_numpy()
        self.atom_type = self.Atoms_df["atype"].to_numpy()
        self.atom_name = self.Atoms_df["aname"].to_numpy()

    def check_pyrimidalisation(self, N_atom, neighbours):
        """Checks whether N.pl3 atoms are indeed planar"""
        neigh1, neigh2, neigh3 = neighbours
        r1 = (N_atom.xyz - neigh1.xyz) / \
            (np.linalg.norm(N_atom.xyz - neigh1.xyz))
        r2 = (N_atom.xyz - neigh2.xyz) / \
            (np.linalg.norm(N_atom.xyz - neigh2.xyz))
        r3 = (N_atom.xyz - neigh3.xyz) / \
            (np.linalg.norm(N_atom.xyz - neigh2.xyz))
        rn = np.cross(r1, r2) / np.linalg.norm(np.cross(r1, r2))
        deviation_from_plane = abs(np.dot(rn, r3))

        # let's say beyond 20 degrees deviation is too much
        ##### FIX DEFINITION!!!
        limit = np.pi/180 * 20
        if deviation_from_plane > limit:
            if N_atom.atom_type == "N.pl3.primary":
                new_type = "N.3.primary"
            elif N_atom.atom_type == "N.pl3.secondary":
                new_type = "N.3.secondary"
            elif N_atom.atom_type == "N.pl3.tertiary":
                new_type = "N.3.tertiary"
            elif N_atom.atom_type == "N.pl3.nitro":
                new_type = "N.pl3.nitro"
            else:
                new_type = "N.3.aniline"
            
            N_atom.atom_type = new_type
            self.atom_type[N_atom.index] = new_type
            self.Atoms_df.at[N_atom.index, "atype"] = new_type

    def find_connecting_atom(self, MEPS, N_atom):
        bonds = MEPS._cml.list_bonds
        conn_atom_names = [bond.atom2 if (bond.atom1 == N_atom.atom_name) else
                           bond.atom1 if (bond.atom2 == N_atom.atom_name) else
                           np.nan for bond in bonds]
        conn_atom_names = [x for x in conn_atom_names if pd.notnull(x)]
        conn_atom = [np.where(self.atom_name == ca)[0][0]
                     for ca in conn_atom_names]
        neighbours = [self.Atom[a] for a in conn_atom]
        return neighbours


class Atom():
    '''Describes each atom with its properties as an individual object'''

    def __init__(self):
        self.index = None
        self.atomic_number = None
        self.xyz = None
        self.atom_type = None
        self.atom_name = None
        self.covalent_radius = None
        self.vdW_radius = None

    def _create_atom_from_df(self, index, atom_row):
        self.index = index
        self.atomic_number = int(atom_row[0])
        self.xyz = np.array(list(atom_row[2:5]))
        self.atom_type = atom_row[5]
        self.atom_name = atom_row[6]
        self.covalent_radius = (mendeleev.element(
            self.atomic_number).covalent_radius)/100
        self.vdW_radius = (mendeleev.element(
            self.atomic_number).vdw_radius)/100

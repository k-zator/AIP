import copy
from mendeleev import element
import pandas as pd
import numpy as np


class PdbReader:
    def __init__(self, pdb_file, skiprows=1):

        """This creates from a mol2 input file both the cml file to input to nwchem, and also a pdb file which contains the same atom names as the cml file. NB that the mol2 atom names are all overwritten, in order to ensure that the atom names do not repeat, as this would be problematic for the cml file
        """
        self.df_atom, self.df_conect = self.read_pdb(pdb_file, skiprows=skiprows)
        self.dict_count=self.get_dict_count(self.df_atom)
        self.list_bonds = self.get_bonds(self.df_conect,self.dict_count)
        self.list_atoms = self.get_atom_list(self.df_atom)
        self.xyz = np.array([i.xyz for i in self.list_atoms])
        self.colors = np.array([i.color for i in self.list_atoms])
        self.list_lines = self.get_list_lines(self.list_atoms, self.list_bonds)
        self.atoms = np.array([a.aname for a in self.list_atoms])
    @staticmethod
    def read_pdb(filename, skiprows=0):
        df = pd.read_fwf(filename, header = None, skiprows=skiprows)
        dict_count = {}
        dict_coords = {}
        count = 1
        list_lines =  []
        df_atom = copy.deepcopy(df[df[0] == "ATOM"])
        for n, i in enumerate(df_atom.iloc[0]):
            if pd.isna(i) == True:
                df_atom.drop(columns=n, inplace=True)
        #changed reset
        df_atom = df_atom.T.reset_index(drop=True).T
        df_conect = df[df[0] == "CONECT"]
        return df_atom, df_conect

    
    def get_atom_list(self, df):
        list_atoms = []
        for n, i in df.iterrows():
            if i[0]=="ATOM":
                atom = self.get_atom(i, n)
                list_atoms.append(atom)
        return list_atoms

    @staticmethod
    def get_dict_count(df):
        dict_count = {}
        count = 1
        for n, i in df.iterrows():
            if i[0]=="ATOM":
                dict_count[count] = f"{i[2]}"
                count += 1
        return dict_count

    def get_atom(atomarray, i, n):
        atom = Atom()
        if len(i[10])==1:
            element_string = i[10]
        else:
            element_string = i[10][0].upper() + i[10][1].lower()
        atom.set_element(element_string)
        atom.set_aname(f"{i[2]}")
        atom.set_color(element_string)
        atom.set_xyz(float(i[5]), float(i[6]), float(i[7]))
        return atom


    def get_bonds(self, df, dict_count):
        list_bonds = []
        for n, i in df.iterrows():
            if i[0]=="CONECT":
                dict_order = {}
                for j in i[2:]:
                    if pd.isna(j)!=True:
                        aname = dict_count[int(j)]
                        if aname not in dict_order.keys():
                            dict_order[f"{aname}"] = 1
                        else:
                            dict_order[f"{aname}"] += 1
                for j, o in dict_order.items():
                    bond = Bond()
                    bond.set_atom1(dict_count[i[1]])
                    bond.set_atom2(j)
                    bond.set_bond_order(o)
                    list_bonds.append(bond)
        return list_bonds
    @staticmethod
    def get_list_lines(list_atoms, list_bonds):
        dict_xyz = {}
        for i in list_atoms:
            dict_xyz[i.aname]=i.xyz
        list_lines = []
        for i in list_bonds:
            list_lines.append((dict_xyz[i.atom1], dict_xyz[i.atom2]))
        return list_lines 




class Bond:
    def __init__(self):
        self.atom1 = None
        self.atom2 = None
        self.bond_order = None
    def set_atom1(self, atom1):
        self.atom1 = atom1
    def set_atom2(self, atom2):
        self.atom2=atom2
    def set_bond_order (self, bond_order):
        self.bond_order = bond_order

class Atom:
    def __init__(self):
        self.aname = None
        self.element = None
        self.xyz = None
    def set_aname(self, aname):
        self.aname=aname
    def set_element(self, elem):
        self.element = elem
    def set_color(self, elem):
        self.color = str(element(elem).cpk_color)
    def set_xyz(self, x, y, z):
        self.xyz = np.array([x,y,z])

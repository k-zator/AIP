"""
Creates instance of AIP class for each individual AIP. Useful for XML file writing.
@author: Katarzyna Joanna Zator (kz265) and Maria Chiara Storer (mcs92)
"""

import numpy as np


class AIP():
    """"This class is to describe a molecule's AIPs, each as its own object."""

    def __init__(self, aip_dict_entry=None):
        if aip_dict_entry is not None:
            self.value = aip_dict_entry[0]
            self.mepsvalue = aip_dict_entry[1]
            self.xyz = aip_dict_entry[2]
            self.mepsindex = aip_dict_entry[3]
            self.type = aip_dict_entry[4]
            self.atom_owner_index = aip_dict_entry[5]
            self.atom_type = aip_dict_entry[6]
            self.atom_name = aip_dict_entry[7]
            self.dual = False
            self.define_fraction()
            if self.type in "polar":
                self.isosurface = 0.0300
            elif self.type == "hydrogen":
                self.isosurface = 0.0104
            elif self.type == "sigma":
                self.isosurface = 0.0104
            else:
                self.isosurface = 0.0020
        else:
            self.value = None
            self.mepsvalue = None
            self.xyz = None
            self.mepsindex = None
            self.type = None
            self.atom_owner_index = None
            self.atom_type = None
            self.atom_name = None
            self.fraction = None
            self.isosurface = None
            self.dual = None

    def set_value(self, value):
        self.value = value

    def set_type(self, aip_type):
        self.type = aip_type
        self.isosurface = 0.0020
    
    def set_nearest_atom(self, nearest_atom):
        self.atom_name = nearest_atom

    def set_atomindex(self, atom_index):
        self.atom_owner_index = atom_index

    def set_xyz(self, x, y, z):
        self.xyz = np.array([[x, y, z]])

    def set_atom_type(self, atom_type):
        self.atom_type = atom_type

    def set_mepsvalue(self, mepsvalue):
        self.mepsvalue = mepsvalue

    def set_mepsindex(self, mepsindex):
        self.mepsindex = mepsindex

    def set_isosurface(self, isosurface):
        self.isosurface = isosurface

    def set_AIP_area_fraction(self, aip_area_fraction):
        self.fraction = aip_area_fraction

    def define_fraction(self):
        if self.type == "sigma" and self.atom_type[:2] == "S.":
            self.fraction = 0.5
        elif self.type == "sigma":
            self.fraction = 1.0
        elif self.type == "outer-sigma":
            self.fraction = 1.0
        elif self.type == "outer-polar":
            self.fraction = 1.0
        elif self.atom_type in ["F"]:
            self.fraction = 1.0
        elif self.atom_type in ["S.2.phene"]:
            self.fraction = 1.0
        elif self.atom_type in ["Cl", "Br"] and self.type == "non-polar":
            self.fraction = 0.5
        elif self.atom_type[:2] in ["N.", "O.", "S.", "C."] and self.type == "non-polar":
            self.fraction = 0.5
        else:
            self.fraction = 1.0

    def set_dual(self, AIP_value_max, AIP_mepsvalue_max, AIP_index_max):
        self.dual = True
        self.mepsvalue_max = AIP_mepsvalue_max
        self.mepsindex_max = AIP_index_max
        self.value_max = AIP_value_max

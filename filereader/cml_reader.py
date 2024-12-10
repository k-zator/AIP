import logging
import copy
from lxml import etree
import math
from mendeleev import element
import numpy as np
import pandas as pd
import networkx as nx
from aip_atom_types.assign_aip_atom_types import assign_aip_atom_types

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.WARN)

CML_NS = "http://www.xml-cml.org/schema"

class CmlReader:
    def __init__(self, cml_file, from_file=True, scale=1, ns=CML_NS):

        """This creates from a mol2 input file both the cml file to input to nwchem, and also a pdb file which contains the same atom names as the cml file. NB that the mol2 atom names are all overwritten, in order to ensure that the atom names do not repeat, as this would be problematic for the cml file
        """
        self.ns = ns
        if from_file == True:
            parser = etree.XMLParser(remove_blank_text=True)
            self.tree = etree.parse(cml_file, parser)
        if from_file == False:
            self.tree = copy.deepcopy(etree.fromstring(cml_file))
        self.scale=scale
        self.list_bonds = self._get_bond_list_()
        self.list_atoms = self._get_atom_list_()
        self.dict_atoms = {i.aname :i for i in self.list_atoms}
        self.xyz = np.array([i.xyz for i in self.list_atoms])
        self.colors = np.array([i.color for i in self.list_atoms])
        self.atoms = np.array([a.aname for a in self.list_atoms])
        try:
            self.inchikey = self._get_inchikey_()
        except:
            pass

    def _get_inchikey_(self):
        molecule_elem = self.tree.xpath("//cml:molecule", namespaces={"cml":self.ns})
        return  molecule_elem[0].attrib['{{{}}}id'.format(self.ns)]


    def _get_atom_list_(self):
        list_atoms = []
        if self.ns != None:
            for cml_elem in self.tree.xpath("//cml:atom", namespaces={"cml":self.ns}):
                atom = self._get_atom_(cml_elem, self.scale)
                list_atoms.append(atom)
        else:
            for cml_elem in self.tree.xpath("//atom"):
                atom = self._get_atom_(cml_elem, self.scale)
                list_atoms.append(atom)
        return list_atoms
    
    def _get_atom_(self, cml_elem, scale):
        atom = Atom()
        if self.ns != None:
            name = cml_elem.attrib['{{{}}}id'.format(self.ns)]
            elem = cml_elem.attrib['{{{}}}elementType'.format(self.ns)]
            x,y,z = scale * float (cml_elem.attrib ['{{{}}}x3'.format(self.ns)]),\
                        scale * float(cml_elem.attrib ['{{{}}}y3'.format(self.ns)]),\
                        scale * float(cml_elem.attrib ['{{{}}}z3'.format(self.ns)])
            try: 
                sybyl = cml_elem.attrib['{{{}}}sybyl'.format(self.ns)]
                atom.set_sybyl(sybyl)
            except:
                pass
            try:
                aipAtomType = cml_elem.attrib['{{{}}}aipAtomType'.format(self.ns)]
                atom.set_aipAtomType(aipAtomType)
            except:
                pass

        else:
            name = cml_elem.attrib['id'.format(self.ns)]
            elem = cml_elem.attrib['elementType'.format(self.ns)]
            x,y,z = scale * float (cml_elem.attrib ['x3'.format(self.ns)]),\
                        scale * float(cml_elem.attrib ['y3'.format(self.ns)]),\
                        scale * float(cml_elem.attrib ['z3'.format(self.ns)])

        atom.set_element(elem)
        atom.set_aname(name)
        atom.set_color(elem)
        atom.set_xyz(x, y, z)
        return atom

    def _get_bond_list_(self):
        list_bonds = []
        if self.ns != None:
            for bond_elem in self.tree.xpath("//cml:bondArray/cml:bond", namespaces={"cml":self.ns}):
                list_bonds.append(self._get_bond_(bond_elem))
        else:
            for bond_elem in self.tree.xpath("//bondArray/bond"):
                list_bonds.append(self._get_bond_(bond_elem))

        return list_bonds

    def _get_bond_(self, bond_elem):
        if self.ns != None:
            bondpair = bond_elem.attrib["{{{}}}atomRefs2".format(self.ns)]
            order = bond_elem.attrib["{{{}}}order".format(self.ns)]
        else:
            bondpair = bond_elem.attrib["atomRefs2".format(self.ns)]
            order = bond_elem.attrib["order".format(self.ns)]
        b1, b2 = bondpair.split(' ')
        bond = Bond()
        bond.set_atom1(b1)
        bond.set_atom2(b2)
        bond.set_bond_order(order)
        return bond

    @staticmethod
    def _compare_atoms_(atom0, atom1):
        x0, y0, z0 = atom0.xyz
        x1, y1, z1 = atom1.xyz

        if (
                math.isclose(x0, x1, abs_tol = 0.001) and\
                math.isclose(y0, y1, abs_tol = 0.001) and\
                math.isclose(z0, z1, abs_tol = 0.001)
                ):
            return True
        else:
            return False


    def _get_sybyl_dict_(self, mol2):
        sybyl_dict = {}
        for i in self.list_atoms:
            for j in mol2.list_atoms:
                if self._compare_atoms_(i, j):
                    sybyl_dict[i.aname] = j.sybyl
        for i in [j.aname for j in self.list_atoms]:
            if i not in sybyl_dict.keys():
                LOGGER.error(f"Atom {i} has not been assigned a sybyl atom type")
        return sybyl_dict

    def get_cml_with_sybyl_info(self, mol2):
        new_tree = copy.deepcopy(self.tree)
        sybyl_dict = self._get_sybyl_dict_(mol2)
        if self.ns != None:
            for cml_elem in new_tree.xpath("//cml:atom", namespaces={"cml":self.ns}):
                new_attrib = "{{{}}}sybyl".format(self.ns)
                atom_attrib = "{{{}}}id".format(self.ns)
                cml_elem.attrib[new_attrib] = sybyl_dict[cml_elem.attrib[atom_attrib]]
        else:
            for cml_elem in new_tree.xpath("//atom"):
                new_attrib = "sybyl"
                atom_attrib = "id"
                cml_elem.attrib[new_attrib] = sybyl_dict[cml_elem.attrib[atom_attrib]]
        return new_tree
    
    def get_cml_network(self, extra_info = False):
        network = nx.Graph()
        if extra_info:
            for atom in self.list_atoms:
                network.add_node(atom.aname, elementType=atom.element, sybyl=atom.sybyl, aipAtomType = atom.aipAtomType)
            for bond in self.list_bonds:
                network.add_edge(bond.atom1, bond.atom2, bondOrder=bond.bond_order)
        else:
            for atom in self.list_atoms:
                network.add_node(atom.aname, elementType=atom.element)
            for bond in self.list_bonds:
                network.add_edge(bond.atom1, bond.atom2, bondOrder=bond.bond_order)
        return network
    
    def get_cml_network_with_sybyl(self, mol2):
        sybyl_dict = self._get_sybyl_dict_(mol2)
        network = nx.Graph()
        for atom in self.list_atoms:
            network.add_node(atom.aname, sybyl=sybyl_dict[atom.aname], elementType=atom.element)
        for bond in self.list_bonds:
            network.add_edge(bond.atom1, bond.atom2, bondOrder=bond.bond_order)
        return network
    
    def get_cml_with_sybyl_and_aip_atom_types(self, mol2, aromatic=False):
        network = self.get_cml_network_with_sybyl(mol2)
        new_tree = copy.deepcopy(self.tree)
        sybyl_dict = self._get_sybyl_dict_(mol2)
        aipatom_dict = assign_aip_atom_types(network, aromatic)
        if self.ns != None:
            for cml_elem in new_tree.xpath("//cml:atom", namespaces={"cml":self.ns}):
                new_attrib = "{{{}}}sybyl".format(self.ns)
                aipatom_attrib = "{{{}}}aipAtomType".format(self.ns)
                atom_attrib = "{{{}}}id".format(self.ns)
                cml_elem.attrib[aipatom_attrib] = aipatom_dict[cml_elem.attrib[atom_attrib]]
        else:
            for cml_elem in new_tree.xpath("//atom"):
                new_attrib = "sybyl"
                atom_attrib = "id"
                aipatom_attrib = "aipAtomType"
                cml_elem.attrib[new_attrib] = sybyl_dict[cml_elem.attrib[atom_attrib]]
                cml_elem.attrib[aipatom_attrib] = aipatom_dict[cml_elem.attrib[atom_attrib]]
        return new_tree


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
    def set_sybyl(self, sybyl):
        self.sybyl = sybyl
    def set_aipAtomType(self, aipAtomType):
        self.aipAtomType = aipAtomType


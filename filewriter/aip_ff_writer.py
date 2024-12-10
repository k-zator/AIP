import networkx as nx
from networkx.algorithms import shortest_path_length
import numpy as np
from lxml import etree
import logging
import copy
import time
import networkx.algorithms.isomorphism as iso
from mendeleev import element
from filereader.aip_reader import AipReader, SSIP_NS, CML_NS
from filewriter.compute_virtual_sites import get_weights, get_anchors


logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.INFO)


class AipFFWriter(AipReader):
    """
    A class used to represent information regarding atoms and ssips described in the ssip.xml files


    Attributes
    ----------
    atom_dict : dict
        A dictionary that describes the atoms present in the molecule. 
        The key is the cml atom name and the value is an SsipAtom object.

    list_aips : list
       A list of SsipSsip objects describing all the ssips present in the molecule.

    list_bonds : list
       A list of SsipBond objects describing the connectivity of the molecule.

    ssip_network : nx.Graph
       A nx.Graph representation of the molecule. 
       The NODES are the cml atom names of the molecule, and these have an attribute elementType which represents the element of the atom
       The EDGES represent which atoms are connected and these have an attribute bondOrder which represents bond order of the bond.

    """

    def __init__(self, ssip_file):
        super().__init__(ssip_file, scale=1)
        self.ligand_tree = None
        self.resname = "MOL"
        self.network = self.get_cml_network()
        self.atom_tree = self.get_atom_tree(ssip_file)
        self.vs_tree = self.get_vs_tree()

    def write_ligand_xml(self, out_file, vs=True):
        if vs:
            with open(out_file, "wb") as writefile:
                writefile.write(etree.tostring(
                    self.vs_tree, pretty_print=True))
                LOGGER.info("written forcefield")
            return
        else:
            with open(out_file, "wb") as writefile:
                writefile.write(etree.tostring(
                    self.atom_tree, pretty_print=True))
                LOGGER.info("written forcefield")
            return

    def write_custom_xml(self, out_file, protein_xml=None, vs=True):
        if vs:
            tree = self.vs_tree
        else:
            tree = self.atom_tree
        custom_tree = self.get_custom_xml(tree, protein_xml=protein_xml)
        with open(out_file, "wb") as writefile:
            writefile.write(etree.tostring(custom_tree, pretty_print=True))
        LOGGER.info("written custom")
        return

    def get_vs_tree(self):
        vs_tree = copy.deepcopy(self.atom_tree)
        root = vs_tree.xpath("//AtomTypes")
        new_element = etree.SubElement(root[0], "Type")
        new_element.attrib["class"] = "ssip"
        new_element.attrib["name"] = "ssip"
        new_element.attrib["mass"] = "0.0"
        root = vs_tree.xpath(
            "//Residues/Residue[@name='{}']/Atom".format(self.resname))
        root_tag = root[0].tag
        root_parent = root[0].getparent()
        self.set_anchors_weights(self.list_aips)
        for i, ssip in enumerate(self.list_aips):
            new_element = etree.SubElement(root_parent, root_tag)
            new_element.attrib["ssip_charge"] = "{}".format(ssip.value)
            name = "M{}_{}".format(i, ssip.atom_type)
            new_element.attrib["fraction"] = str(ssip.fraction)
            new_element.attrib["isosurface"] = str(ssip.isosurface)
            new_element.attrib["name"] = name
            new_element.attrib["type"] = "ssip"
            vs_element = etree.SubElement(root_parent, "VirtualSite")
            vs_element.attrib["type"] = "outOfPlane"
            vs_element.attrib["siteName"] = name
            vs_element.attrib["atomName1"] = ssip.anchor1.aname
            vs_element.attrib["atomName2"] = ssip.anchor2.aname
            vs_element.attrib["atomName3"] = ssip.anchor3.aname
            vs_element.attrib["weight12"] = str(round(ssip.weights[0], 3))
            vs_element.attrib["weight13"] = str(round(ssip.weights[1], 3))
            vs_element.attrib["weightCross"] = str(
                round(ssip.weights[2], 3)*10)  # unit conversion Angstrom to nm
        return vs_tree

    def get_atom_tree(self, ssip_file):
        atom_tree = copy.deepcopy(self.get_xml_tree())
        for k, v in self.dict_atoms.items():
            atom_name = f""
            root = atom_tree.xpath(
                '//Residues/Residue[@name="{}"]/Atom[@name="{}"]'.format(self.resname, k.upper()))
        return atom_tree

    def get_xml_tree(self):
        tree = etree.Element("ForceField")
        atype = etree.SubElement(tree, "AtomTypes")
        for k, v in self.dict_atoms.items():
            types = etree.SubElement(atype, "Type")
            types.attrib["class"] = "{}".format(k)
            types.attrib["name"] = "ligand-{}".format(k)
            types.attrib["element"] = "{}".format(v.element)
            types.attrib["mass"] = "{}".format(
                element(v.element).atomic_weight)
        residues = etree.SubElement(tree, "Residues")
        residue = etree.SubElement(residues, "Residue")
        residue.attrib["name"] = self.resname
        for k, v in self.dict_atoms.items():
            atom_elem = etree.SubElement(residue, "Atom")
            atom_elem.attrib["name"] = k
            atom_elem.attrib["type"] = "ligand-{}".format(k)
            atom_elem.attrib["ssip_charge"] = "0"
            atom_elem.attrib["fraction"] = "0"
            atom_elem.attrib["isosurface"] = "0"
        for bond in self.list_bonds:
            bond_elem = etree.SubElement(residue, "Bond")
            bond_elem.attrib["atomName1"] = bond.atom1
            bond_elem.attrib["atomName2"] = bond.atom2
        return tree

    def get_custom_xml(self, ligand_tree, protein_xml=None):
        """This function takes in the ligand force field file
        and the protein force field file and outputs a CustomNonBonded
        force field which is necessary to assign the ssip value to each atom present in the
        protein ligand system of interest."""

        tree_new = etree.Element("ForceField")
        custom_elem = etree.SubElement(tree_new, "CustomNonbondedForce")
        custom_elem.attrib["energy"] = "n_zeros*SA*SAneg*SApos*ssip_charge1*ssip_charge2"
        custom_elem.attrib["bondCutoff"] = "1"
        perpar_elem = etree.SubElement(custom_elem, "PerParticleParameter")
        perpar_elem.attrib["name"] = "ssip_charge"
        perpar_elem = etree.SubElement(custom_elem, "PerParticleParameter")
        perpar_elem.attrib["name"] = "fraction"
        perpar_elem = etree.SubElement(custom_elem, "PerParticleParameter")
        perpar_elem.attrib["name"] = "isosurface"
        global_elem = etree.SubElement(custom_elem, "GlobalParameter")
        global_elem.attrib["name"] = "SA"
        global_elem.attrib["defaultValue"] = str(self.surface.total)
        global_elem = etree.SubElement(custom_elem, "GlobalParameter")
        global_elem.attrib["name"] = "SApos"
        global_elem.attrib["defaultValue"] = str(self.surface.positive)
        global_elem = etree.SubElement(custom_elem, "GlobalParameter")
        global_elem.attrib["name"] = "SAneg"
        global_elem.attrib["defaultValue"] = str(self.surface.negative)
        global_elem = etree.SubElement(custom_elem, "GlobalParameter")
        global_elem.attrib["name"] = "n_zeros"
        global_elem.attrib["defaultValue"] = str(0)
        attfrom_elem = etree.SubElement(custom_elem, "UseAttributeFromResidue")
        attfrom_elem.attrib["name"] = "ssip_charge"
        attfrom_elem = etree.SubElement(custom_elem, "UseAttributeFromResidue")
        attfrom_elem.attrib["name"] = "fraction"
        attfrom_elem = etree.SubElement(custom_elem, "UseAttributeFromResidue")
        attfrom_elem.attrib["name"] = "isosurface"
        for i in ligand_tree.xpath("//AtomTypes/Type"):
            atom_elem = etree.SubElement(custom_elem, "Atom")
            atom_elem.attrib["type"] = i.attrib["name"]
        if protein_xml != None:
            tree_prot = copy.deepcopy(etree.parse(protein_xml))
            for i in tree_prot.xpath("//AtomTypes/Type"):
                atom_elem = etree.SubElement(custom_elem, "Atom")
                atom_elem.attrib["type"] = i.attrib["name"]
        return tree_new

    staticmethod

    def set_anchors_weights(self, list_aips):
        virtual_site_ssip_list = []
        ssip_none = []
        for aip in list_aips:
            anchor1, anchor2, anchor3 = get_anchors(
                aip.atom_neigh, self.network, self.dict_atoms)
            weights = get_weights(aip, anchor1, anchor2, anchor3)
            aip.set_anchors(anchor1, anchor2, anchor3)
            aip.set_weights(weights)

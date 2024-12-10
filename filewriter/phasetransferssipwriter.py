from lxml import etree
import copy
import os
from filereader.aip_reader import AipReader

CML_NAMESPACE_DICT = {
    "ssip": "http://www-hunter.ch.cam.ac.uk/SSIP",
    "cml": "http://www.xml-cml.org/schema",
    "xsi":"http://www.w3.org/2001/XMLSchema-instance"
}

CML_NAME = "{{{}}}".format(CML_NAMESPACE_DICT["cml"])


SSIP_NAME = "{{{}}}".format(CML_NAMESPACE_DICT["ssip"])

XSI_NAME = "{{{}}}".format(CML_NAMESPACE_DICT["xsi"])

class PhaseTransferSsipWriter():
    """The aim of this class is to write out the AIP files in the format that the solvation code can understand.
    The output of the filweriter are two SSIP files, the first one with all the AIPs, 
    so AIPs with both whole and half fractions. 
    The second file only has the SSIPs that hava  a whole fraction.
    When using these files in phase calculator you need to use both files and halve
    the molefraction of the solvent for each of the files (because the two files together 
    contribute as the whole solvet)"""
    def __init__(self, ssip_file):
        parser = etree.XMLParser(remove_blank_text=True)
        self.aip= AipReader(ssip_file)
        super_element = etree.Element(SSIP_NAME+"SSIPMolecule", nsmap=CML_NAMESPACE_DICT)
        super_element.attrib[SSIP_NAME + "ssipSoftwareVersion"] = "0.0.0-newssip"
        super_element.attrib[SSIP_NAME + "parameterVersion"] = "3.0.0"
        super_element.attrib[XSI_NAME + "schemaLocation"] = "http://www.xml-cml.org/schema http://www-hunter.ch.cam.ac.uk/schema/cmlschema.xsd http://www-hunter.ch.cam.ac.uk/SSIP http://www-hunter.ch.cam.ac.uk/schema/SSIP.xsd"
        self.tree = etree.ElementTree(super_element)
    
    def add_molecule(self, tree):
        new_tree = copy.deepcopy(tree)
        root = new_tree.getroot()
        molecule_element = etree.SubElement(root, CML_NAME+"molecule", nsmap=CML_NAMESPACE_DICT)
        molecule_element.attrib[CML_NAME + "id"] = str(self.aip.inchikey)
        molecule_element.attrib[SSIP_NAME + "stdInChIKey"] = str(self.aip.inchikey)
        self.add_atom_array(molecule_element)
        self.add_bond_array(molecule_element)
        return new_tree

    def add_atom_array(self, molecule_element):
        atom_array_element = etree.SubElement(molecule_element, CML_NAME+"atomArray", nsmap=CML_NAMESPACE_DICT)
        for atom in self.aip.list_atoms:
            self.add_atom_element(atom_array_element, atom)
        return 
    
        
    @staticmethod
    def add_atom_element(atom_array_element, atom):
        element = etree.SubElement(atom_array_element, CML_NAME+"atom", nsmap=CML_NAMESPACE_DICT)
        element.attrib[CML_NAME + "elementType"] = str(atom.element)
        element.attrib[CML_NAME + "id"] = str(atom.aname)
        element.attrib[CML_NAME + "x3"] = str(atom.xyz[0])
        element.attrib[CML_NAME + "y3"] = str(atom.xyz[1])
        element.attrib[CML_NAME + "z3"] = str(atom.xyz[2])
        return
    
    def add_bond_array(self, molecule_element):
        bond_array_element = etree.SubElement(molecule_element, CML_NAME+"bondArray", nsmap=CML_NAMESPACE_DICT)
        for bond in self.aip.list_bonds:
            self.add_bond_element(bond_array_element, bond)
        return 
    
        
    @staticmethod
    def add_bond_element(bond_array_element, bond):
        element = etree.SubElement(bond_array_element, CML_NAME+"bond", nsmap=CML_NAMESPACE_DICT)
        atom1 = bond.atom1
        atom2 = bond.atom2
        element.attrib[CML_NAME + "atomRefs2"] = f"{atom1} {atom2}"
        element.attrib[CML_NAME + "order"] = str(bond.bond_order)
        return
        
    def add_surface_information(self, tree, surface):
        new_tree = copy.deepcopy(tree)
        root = new_tree.getroot()
        aips_element = etree.SubElement(root, SSIP_NAME+"SurfaceInformation", nsmap=CML_NAMESPACE_DICT)
        surfaces_element = etree.SubElement(aips_element, SSIP_NAME+"Surfaces", nsmap=CML_NAMESPACE_DICT)
        surface_element = etree.SubElement(surfaces_element, SSIP_NAME+"Surface", nsmap=CML_NAMESPACE_DICT)
        total_element = etree.SubElement(surface_element, SSIP_NAME+"TotalSurfaceArea", nsmap=CML_NAMESPACE_DICT)
        total_element.text = str(surface.total)
        total_element.attrib[SSIP_NAME + "unit"] = "Å^2"
        negative_element = etree.SubElement(surface_element, SSIP_NAME+"NegativeSurfaceArea", nsmap=CML_NAMESPACE_DICT)
        negative_element.text = str(surface.negative)
        negative_element.attrib[SSIP_NAME + "unit"] = "Å^2"
        positive_element = etree.SubElement(surface_element, SSIP_NAME+"PositiveSurfaceArea", nsmap=CML_NAMESPACE_DICT)
        positive_element.text = str(surface.positive)
        positive_element.attrib[SSIP_NAME + "unit"] = "Å^2"
        edensity_element = etree.SubElement(surface_element, SSIP_NAME+"ElectronDensityIsosurface", nsmap=CML_NAMESPACE_DICT)
        edensity_element.text = str(0.002)
        edensity_element.attrib[SSIP_NAME + "unit"] = "e bohr^-3"
        number_meps_element = etree.SubElement(surface_element, SSIP_NAME+"NumberOFMEPSPoints", nsmap=CML_NAMESPACE_DICT)
        number_meps_element.text = str(surface.numberOFMEPSPoints)
        vdw_volume_element = etree.SubElement(surface_element, SSIP_NAME+"VdWVolume", nsmap=CML_NAMESPACE_DICT)
        vdw_volume_element.text = str(surface.VdWVolume)
        vdw_volume_element.attrib[SSIP_NAME + "unit"] = "Å^3"
        maximum_element = etree.SubElement(surface_element, SSIP_NAME+"ElectrostaticPotentialMax", nsmap=CML_NAMESPACE_DICT)
        maximum_element.text = str(surface.electrostaticPotentialMax)
        maximum_element.attrib[SSIP_NAME + "unit"] = "hartree"
        minimum_element =  etree.SubElement(surface_element, SSIP_NAME+"ElectrostaticPotentialMin", nsmap=CML_NAMESPACE_DICT)
        minimum_element.text = str(surface.electrostaticPotentialMin)
        minimum_element.attrib[SSIP_NAME + "unit"] = "hartree"

        return new_tree

    def add_aip_information(self, tree, list_aips):
        new_tree = copy.deepcopy(tree)
        root = new_tree.getroot()
        aips_element = etree.SubElement(root, SSIP_NAME+"SSIPs", nsmap=CML_NAMESPACE_DICT)
        for aip in list_aips:
            self.add_aip_element(aips_element, aip)
        return new_tree
    
    @staticmethod
    def add_aip_element(aips_elem, aip):
        element = etree.SubElement(aips_elem, SSIP_NAME+"SSIP", nsmap=CML_NAMESPACE_DICT)
        element.attrib[SSIP_NAME + "value"] = str(aip.value)
        element.attrib[SSIP_NAME + "nearestAtomID"] = str(aip.atom_neigh.aname)
        element.attrib[CML_NAME + "x3"] = str(aip.xyz[0])
        element.attrib[CML_NAME + "y3"] = str(aip.xyz[1])
        element.attrib[CML_NAME + "z3"] = str(aip.xyz[2])
        return
    
    
    def get_aip_tree(self, surface, list_aips):
        tree = copy.deepcopy(self.tree)
        molecule_tree = self.add_molecule(tree)
        surface_tree = self.add_surface_information(molecule_tree, surface)
        aip_tree = self.add_aip_information(surface_tree, list_aips)
        return aip_tree

    def get_all_aips_tree(self):
        new_tree = self.get_aip_tree(self.aip.surface, self.aip.list_aips)
        return new_tree
    
    def get_whole_aips_tree(self):
        list_whole = [i for i in self.aip.list_aips if i.aipAreaFraction == 1.0]
        new_tree = self.get_aip_tree(self.aip.surface, list_whole)
        return new_tree
    
    def get_split_half_aips_tree(self):
        new_tree0 = self.get_aip_tree(self.aip.surface, [self.aip.list_aips[0]])
        new_tree1 = self.get_aip_tree(self.aip.surface, self.aip.list_aips[1:])
        return new_tree0, new_tree1
    
    def write_all_aips_file(self, filename):
        new_tree = self.get_all_aips_tree()
        new_tree.write(filename, encoding="UTF-8", xml_declaration=True, pretty_print=True)        

    def write_whole_aips_file(self, filename):
        new_tree = self.get_whole_aips_tree()
        new_tree.write(filename, encoding="UTF-8", xml_declaration=True, pretty_print=True)
        
    def write_split_half_aips_file(self, filename0, filename1):
        new_tree0, new_tree1 = self.get_split_half_aips_tree()
        new_tree0.write(filename0, encoding="UTF-8", xml_declaration=True, pretty_print=True)
        new_tree1.write(filename1, encoding="UTF-8", xml_declaration=True, pretty_print=True)
    
    def write_files(self, directory):
        isExist = os.path.exists(directory)
        if not isExist:
            os.makedirs(directory)
        filename_a = self.aip.inchikey[:26] + "A" + "_ssip.xml"
        filename_b = self.aip.inchikey[:26] + "B" + "_ssip.xml"
        list_whole = [i for i in self.aip.list_aips if i.aipAreaFraction == 1.0]
        if len(list_whole) > 0: 
            self.write_all_aips_file(directory + "/" + filename_a)
            self.write_whole_aips_file(directory + "/" + filename_b)
        else:
            self.write_split_half_aips_file(directory + "/" + filename_a, directory + "/" + filename_b)


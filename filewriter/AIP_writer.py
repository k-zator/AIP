from lxml import etree
import copy

CML_NAMESPACE_DICT = {
    "ssip": "http://www-hunter.ch.cam.ac.uk/SSIP",
    "cml": "http://www.xml-cml.org/schema",
    "xsi":"http://www.w3.org/2001/XMLSchema-instance"
}

CML_NAME = "{{{}}}".format(CML_NAMESPACE_DICT["cml"])

SSIP_NAME = "{{{}}}".format(CML_NAMESPACE_DICT["ssip"])

XSI_NAME = "{{{}}}".format(CML_NAMESPACE_DICT["xsi"])

class AIPWriter():
    def __init__(self, cml_file):
        parser = etree.XMLParser(remove_blank_text=True)
        parsed_tree = copy.deepcopy(etree.parse(cml_file, parser))
        super_element = etree.Element(SSIP_NAME+"SSIPMolecule", nsmap=CML_NAMESPACE_DICT)
        super_element.attrib[SSIP_NAME + "ssipSoftwareVersion"] = "0.0.0-newssip"
        super_element.attrib[SSIP_NAME + "parameterVersion"] = "3.0.0"
        super_element.attrib[XSI_NAME + "schemaLocation"] = "http://www.xml-cml.org/schema http://www-hunter.ch.cam.ac.uk/schema/cmlschema_KMC.xsd http://www-hunter.ch.cam.ac.uk/SSIP http://www-hunter.ch.cam.ac.uk/schema/SSIP_KMC.xsd"
        super_element.append(parsed_tree.getroot())
        self.tree = etree.ElementTree(super_element)
        self.inchikey = self.get_inchikey()

    def get_inchikey(self):
        molecule_elem = self.tree.xpath("//cml:molecule", namespaces=CML_NAMESPACE_DICT)
        return  molecule_elem[0].attrib['{}id'.format(CML_NAME)]

    def add_surface_information(self, tree, surface):
        new_tree = copy.deepcopy(tree)
        root = new_tree.getroot()
        aips_element = etree.SubElement(root, SSIP_NAME+"SurfaceInformation", nsmap=CML_NAMESPACE_DICT)
        surfaces_element = etree.SubElement(aips_element, SSIP_NAME+"Surfaces", nsmap=CML_NAMESPACE_DICT)
        surface_element = etree.SubElement(surfaces_element, SSIP_NAME+"Surface", nsmap=CML_NAMESPACE_DICT)
        total_element = etree.SubElement(surface_element, SSIP_NAME+"TotalSurfaceArea", nsmap=CML_NAMESPACE_DICT)
        total_element.text = str(surface.total)
        total_element.attrib[SSIP_NAME + "unit"] = "Å^2"
        positive_element = etree.SubElement(surface_element, SSIP_NAME+"PositiveSurfaceArea", nsmap=CML_NAMESPACE_DICT)
        positive_element.text = str(surface.positive)
        positive_element.attrib[SSIP_NAME + "unit"] = "Å^2"
        negative_element = etree.SubElement(surface_element, SSIP_NAME+"NegativeSurfaceArea", nsmap=CML_NAMESPACE_DICT)
        negative_element.text = str(surface.negative)
        negative_element.attrib[SSIP_NAME + "unit"] = "Å^2"
        positive_polar_element = etree.SubElement(surface_element, SSIP_NAME+"PositivePolarSurfaceArea", nsmap=CML_NAMESPACE_DICT)
        positive_polar_element.text = str(surface.positive_polar)
        positive_polar_element.attrib[SSIP_NAME + "unit"] = "Å^2"
        negative_polar_element = etree.SubElement(surface_element, SSIP_NAME+"NegativePolarSurfaceArea", nsmap=CML_NAMESPACE_DICT)
        negative_polar_element.text = str(surface.negative_polar)
        negative_polar_element.attrib[SSIP_NAME + "unit"] = "Å^2"
        total_nonpolar_element = etree.SubElement(surface_element, SSIP_NAME+"TotalNonPolarSurfaceArea", nsmap=CML_NAMESPACE_DICT)
        total_nonpolar_element.text = str(surface.positive_non_polar+surface.negative_non_polar)
        total_nonpolar_element.attrib[SSIP_NAME + "unit"] = "Å^2"
        positive_nonpolar_element = etree.SubElement(surface_element, SSIP_NAME+"PositiveNonPolarSurfaceArea", nsmap=CML_NAMESPACE_DICT)
        positive_nonpolar_element.text = str(surface.positive_non_polar)
        positive_nonpolar_element.attrib[SSIP_NAME + "unit"] = "Å^2"
        negative_nonpolar_element = etree.SubElement(surface_element, SSIP_NAME+"NegativeNonPolarSurfaceArea", nsmap=CML_NAMESPACE_DICT)
        negative_nonpolar_element.text = str(surface.negative_non_polar)
        negative_nonpolar_element.attrib[SSIP_NAME + "unit"] = "Å^2"
        number_meps_element = etree.SubElement(surface_element, SSIP_NAME+"NumberOFMEPSPoints", nsmap=CML_NAMESPACE_DICT)
        number_meps_element.text = str(surface.numberOFMEPSPoints)
        maximum_element = etree.SubElement(surface_element, SSIP_NAME+"ElectrostaticPotentialMax", nsmap=CML_NAMESPACE_DICT)
        maximum_element.text = str(surface.electrostaticPotentialMax)
        maximum_element.attrib[SSIP_NAME + "unit"] = "hartree"
        minimum_element =  etree.SubElement(surface_element, SSIP_NAME+"ElectrostaticPotentialMin", nsmap=CML_NAMESPACE_DICT)
        minimum_element.text = str(surface.electrostaticPotentialMin)
        minimum_element.attrib[SSIP_NAME + "unit"] = "hartree"
        vdw_volume_element = etree.SubElement(surface_element, SSIP_NAME+"VdWVolume", nsmap=CML_NAMESPACE_DICT)
        vdw_volume_element.text = str(surface.VdWVolume)
        vdw_volume_element.attrib[SSIP_NAME + "unit"] = "Å^3"
        return new_tree

    def add_aip_information(self, tree, list_aips, dual=False):
        new_tree = copy.deepcopy(tree)
        root = new_tree.getroot()
        aips_element = etree.SubElement(root, SSIP_NAME+"SSIPs", nsmap=CML_NAMESPACE_DICT)
        for i, aip in enumerate(list_aips):
            if dual is not False:
                isdual = dual[i]
            else:
                isdual = False
            if isdual:
                self.add_aip_element(aips_element, aip, dual=True)
            else:
                self.add_aip_element(aips_element, aip)
        return new_tree

    def get_surface_aip_tree(self, surface, list_aips, dual=False):
        surface_tree = self.add_surface_information(self.tree, surface)
        aip_tree = self.add_aip_information(surface_tree, list_aips, dual)
        return aip_tree

    def write_file(self, surface, list_aips, filename, dual=False):
        if dual is not False: 
            new_tree = self.get_surface_aip_tree(surface, list_aips, dual)
        else:
            new_tree = self.get_surface_aip_tree(surface, list_aips)
        new_tree.write(filename, encoding="UTF-8", xml_declaration=True, pretty_print=True)
    
    
    @staticmethod
    def add_aip_element(aips_elem, aip, dual=False):
        element = etree.SubElement(aips_elem, SSIP_NAME+"SSIP", nsmap=CML_NAMESPACE_DICT)
        element.attrib[SSIP_NAME + "value"] = str(aip.value)
        element.attrib[SSIP_NAME + "nearestAtomID"] = str(aip.atom_name)
        element.attrib[SSIP_NAME + "aipAtomType"] = str(aip.atom_type)
        element.attrib[CML_NAME + "x3"] = str(aip.xyz[0][0])
        element.attrib[CML_NAME + "y3"] = str(aip.xyz[0][1])
        element.attrib[CML_NAME + "z3"] = str(aip.xyz[0][2])
        element.attrib[SSIP_NAME +"MEPvalue"] = str(aip.mepsvalue)
        element.attrib[SSIP_NAME + "isosurface"] = str(aip.isosurface)
        element.attrib[SSIP_NAME + "aipAreaFraction"] = str(aip.fraction)
        try:
            element.attrib[SSIP_NAME + "originInChIKey"] = str(aip.origin_inchikey)
        except:
            pass
        if dual:
            element.attrib[SSIP_NAME + "valueMax"] = str(aip.value_max)
            element.attrib[SSIP_NAME +"MEPvalueMax"] = str(aip.mepsvalue_max)
        return


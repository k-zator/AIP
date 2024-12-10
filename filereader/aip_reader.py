from filereader.cml_reader import CmlReader, CML_NS 
import numpy as np
import networkx as nx

SSIP_NS = "http://www-hunter.ch.cam.ac.uk/SSIP"

class AipReader(CmlReader):
    def __init__(self, ssip_file, scale=1, read_surface=True):
        super().__init__(ssip_file, scale)
        self.list_aips = self.get_list_aips()
        if read_surface:
            self.surface = self.get_surface()

    def get_surface(self):
        surface_path=self.tree.xpath("//ssip:SurfaceInformation/ssip:Surfaces/ssip:Surface", namespaces={"cml":CML_NS, "ssip":SSIP_NS})
        return AipSurface(surface_path[0])

    def get_list_aips(self):
        list_aips = []
        for aip_elem in self.tree.xpath("//ssip:SSIP", namespaces={"cml":CML_NS, "ssip":SSIP_NS}):
            aip = self.get_aip(aip_elem)
            list_aips.append(aip)
        return list_aips
    
    def get_aip(self, aip_elem):
        aip = Aip()
        atom_name =aip_elem.attrib[f"{{{SSIP_NS}}}nearestAtomID"]
        aip.set_atom_neigh(self.dict_atoms[atom_name])
        aip.set_value(float(aip_elem.attrib[f"{{{SSIP_NS}}}value"]))
        aip.set_atom_name(atom_name)
        x = float( aip_elem.attrib[f"{{{CML_NS}}}x3"])
        y = float( aip_elem.attrib[f"{{{CML_NS}}}y3"])
        z = float( aip_elem.attrib[f"{{{CML_NS}}}z3"])
        aip.set_xyz(x, y, z)
        aip.set_atom_type(aip_elem.attrib[f"{{{SSIP_NS}}}aipAtomType"])
        aip.set_mepsvalue(float(aip_elem.attrib[f"{{{SSIP_NS}}}MEPvalue"]))
        aip.set_isosurface(float(aip_elem.attrib[f"{{{SSIP_NS}}}isosurface"]))
        aip.set_fraction(float(aip_elem.attrib[f"{{{SSIP_NS}}}aipAreaFraction"]))
        return aip

    def get_cml_network_all(self):
        network = nx.Graph()
        network.name = self.inchikey
        for atom in self.list_atoms:
            network.add_node(atom.aname, sybyl=atom.sybyl, elementType=atom.element, aipAtomType = atom.aipAtomType)
        for bond in self.list_bonds:
            network.add_edge(bond.atom1, bond.atom2, bondOrder=bond.bond_order)
        return network

class Aip:
    #def __init__(self):

    
    def set_value(self, value):
        self.value = value
    
    def set_atom_name(self, atom_name):
        self.atom_name = atom_name

    def set_xyz(self, x, y, z):
        self.xyz = np.array([[x, y, z]])

    def set_atom_type(self, atom_type):
        self.atom_type = atom_type

    def set_mepsvalue(self, mepsvalue):
        self.mepsvalue = mepsvalue

    def set_isosurface(self, isosurface):
        self.isosurface = isosurface

    def set_fraction(self, fraction):
        self.fraction = fraction


    def set_atom_neigh(self, atom_neigh):
        self.atom_neigh = atom_neigh

    def set_anchors(self, anchor1, anchor2, anchor3):
        self.anchor1 = anchor1
        self.anchor2 = anchor2
        self.anchor3 = anchor3

    def set_weights(self, weights):
        self.weights = weights

    def set_origin_inchikey(self, origin_inchikey):
        self.origin_inchikey = origin_inchikey

class AipSurface:
        """This is a class that helps to organise information regarding the surface of the molecule as described  in the ssip.xml file
    
    
    Attributes
    ----------
    tot : float 
        Total surface area

    pos : float 
        Positive surface area
    
    neg : float 
        Negative surface area

    elec_dens_iso : float 
        The electron density of the isosurface considered find the electrostatic potential

    n_meps : int 

    """

        def __init__ (self, surface_elem):

            tot_path = surface_elem.xpath("ssip:TotalSurfaceArea", namespaces={"ssip":SSIP_NS, "cml":CML_NS})
            self.total = float(tot_path[0].text)

            pos_path = surface_elem.xpath("ssip:PositiveSurfaceArea", namespaces={"ssip":SSIP_NS, "cml":CML_NS})
            self.positive = float(pos_path[0].text)


            neg_path = surface_elem.xpath("ssip:NegativeSurfaceArea", namespaces={"ssip":SSIP_NS, "cml":CML_NS})
            self.negative = float(neg_path[0].text)

            pos_polar_path = surface_elem.xpath("ssip:PositivePolarSurfaceArea", namespaces={"ssip":SSIP_NS, "cml":CML_NS})
            self.positive_polar = float(pos_polar_path[0].text)

            neg_polar_path = surface_elem.xpath("ssip:NegativePolarSurfaceArea", namespaces={"ssip":SSIP_NS, "cml":CML_NS})
            self.negative_polar = float(neg_polar_path[0].text)
            
            tot_nonpolar_path = surface_elem.xpath("ssip:TotalNonPolarSurfaceArea", namespaces={"ssip":SSIP_NS, "cml":CML_NS})
            self.total_nonpolar = float(tot_nonpolar_path[0].text)
            
            pos_nonpolar_path = surface_elem.xpath("ssip:PositiveNonPolarSurfaceArea", namespaces={"ssip":SSIP_NS, "cml":CML_NS})
            self.positive_nonpolar = float(pos_nonpolar_path[0].text)

            neg_nonpolar_path = surface_elem.xpath("ssip:NegativeNonPolarSurfaceArea", namespaces={"ssip":SSIP_NS, "cml":CML_NS})
            self.negative_nonpolar = float(neg_nonpolar_path[0].text)
            
            mep_points_path = surface_elem.xpath("ssip:NumberOFMEPSPoints", namespaces={"ssip":SSIP_NS, "cml":CML_NS})
            self.numberOFMEPSPoints = int(mep_points_path[0].text)

            emax_path = surface_elem.xpath("ssip:ElectrostaticPotentialMax", namespaces={"ssip":SSIP_NS, "cml":CML_NS})
            self.electrostaticPotentialMax = float(emax_path[0].text)

            emin_path = surface_elem.xpath("ssip:ElectrostaticPotentialMin", namespaces={"ssip":SSIP_NS, "cml":CML_NS})
            self.electrostaticPotentialMin = float(emin_path[0].text)

            try:
                vdw_volume_path = surface_elem.xpath("ssip:VdWVolume", namespaces={"ssip":SSIP_NS, "cml":CML_NS})
                self.VdWVolume = float(vdw_volume_path[0].text)
            except:
                self.VdWVolume = None



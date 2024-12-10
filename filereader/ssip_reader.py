from filereader.cml_reader import CmlReader, CML_NS 
import numpy as np


SSIP_NS = "http://www-hunter.ch.cam.ac.uk/SSIP"

class SsipReader(CmlReader):
    def __init__(self, ssip_file, scale=1):
        super().__init__(ssip_file, scale)
        self.dict_atoms = {i.aname: i for i in self.list_atoms}
        self.list_ssips = self.get_list_ssips()
        self.surface = self.get_surface()

    def get_surface(self):
        surface_path=self.tree.xpath("//ssip:SurfaceInformation/ssip:Surfaces/ssip:Surface", namespaces={"cml":CML_NS, "ssip":SSIP_NS})
        return SsipSurface(surface_path[0])

    def get_list_ssips(self):
        list_ssips = []
        for ssip_elem in self.tree.xpath("//ssip:SSIP", namespaces={"cml":CML_NS, "ssip":SSIP_NS}):
            ssip = self.get_ssip(ssip_elem)
            list_ssips.append(ssip)
        return list_ssips
    
    def get_ssip(self, ssip_elem):
        ssip = Ssip()
        atom_neigh_aname = ssip_elem.attrib[f"{{{SSIP_NS}}}nearestAtomID"]
        ssip.set_atom_neigh(self.dict_atoms[atom_neigh_aname])
        ssip.set_value(float(ssip_elem.attrib[f"{{{SSIP_NS}}}value"]))
        x = float( ssip_elem.attrib[f"{{{CML_NS}}}x3"])
        y = float( ssip_elem.attrib[f"{{{CML_NS}}}y3"])
        z = float( ssip_elem.attrib[f"{{{CML_NS}}}z3"])
        ssip.set_xyz(x, y, z)
        return ssip

class Ssip:
    def __init__(self):
        self.atom_neigh = None
        self.value = None
    def set_atom_neigh(self, atom_neigh):
        self.atom_neigh = atom_neigh
    def set_value(self, value):
        self.value = value
    def set_xyz(self, x, y, z):
        self.xyz = np.array([x,y,z])

class SsipSurface:
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
        The electron denisty of the isosurface considered find the electrostatic potential

    n_meps : int 

    """

        def __init__ (self, surface_elem):

            tot_path = surface_elem.xpath("ssip:TotalSurfaceArea", namespaces={"ssip":SSIP_NS, "cml":CML_NS})
            self.tot = tot_path[0].text

            pos_path = surface_elem.xpath("ssip:PositiveSurfaceArea", namespaces={"ssip":SSIP_NS, "cml":CML_NS})
            self.pos = pos_path[0].text


            neg_path = surface_elem.xpath("ssip:NegativeSurfaceArea", namespaces={"ssip":SSIP_NS, "cml":CML_NS})
            self.neg = neg_path[0].text

            vol_path = surface_elem.xpath("ssip:VdWVolume", namespaces={"ssip":SSIP_NS, "cml":CML_NS})
            self.vol = vol_path[0].text

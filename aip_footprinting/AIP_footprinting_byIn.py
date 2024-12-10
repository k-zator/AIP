"""
Parses foot-printing as class.
@author: Katarzyna Joanna Zator (kz265)
"""

import logging
from aip_footprinting.read_MEPS import MEPS
from aip_footprinting.atom_class import AtomSet
from aip_footprinting.surface_class import Surface
from aip_footprinting.cross_check_dict import ccd
from aip_footprinting.footprinting import Footprinting
from filewriter.AIP_writer import AIPWriter

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.WARN)


class KMC_footprint_byIn():
    # MEPS gets read and processed, surface quantified and divided
    def __init__(self, inchikey, centre_surface_percentile=90, lp_excl_r=1.5, dualAIP=False,
                 path="/home/kate/workspace/newssip/dataset"):
        """Class for handling footprinting to find AIPs in a .ipy format. Parallel to aip_footprinting_script,
           but it handles molecule names, changed by ccd (dictionary) and reverts to default path to look for 
           the files for ease of handling."""
        self.dualAIP = dualAIP
        try:
            cube_np_file = "{}/0.0020/{}_0.0020_merged.cube".format(
                path, inchikey)
            cube_m_file = "{}/0.0104/{}_0.0104_merged.cube".format(
                path, inchikey)
            cml_file = "{}/cml/{}.cml".format(path, inchikey)
        except Exception:
            LOGGER.error("\n Could not localise cube and cml files, check the cross_check_mol/{} \
                          directory for missing files".format(inchikey))
            exit(1)
        try:
            cube_p_file = "{}/0.0300/{}_0.0300_merged.cube".format(
                path, inchikey)
        except Exception:
            LOGGER.warn("\n Could not localise the polar cube files")

        self.MEPS_np = MEPS(cube_np_file, cml_file)
        self.Atom = AtomSet()
        self.Atom._create_from_MEPS(self.MEPS_np)
        try:
            self.MEPS_m = MEPS(cube_m_file, cml_file)
            if not all(self.MEPS_np.Atoms_df == self.MEPS_m.Atoms_df):
                LOGGER.info(
                    "\n The middle cube file contains different atom coordinates and this will cause errors")
        except Exception:
            self.MEPS_m = 0
        try:
            self.MEPS_p = MEPS(cube_p_file, cml_file, own_dist_scaled=False)
            if not all(self.MEPS_np.Atoms_df == self.MEPS_p.Atoms_df):
                LOGGER.info(
                    "\n The polar cube file contains different atom coordinates and this will cause errors")
        except Exception:
            self.MEPS_p = 0

        self.AIP = Footprinting(self.MEPS_np, self.MEPS_p,
                                self.MEPS_m, self.Atom, centre_surface_percentile, lp_excl_r, dualAIP)

        self.Surface = Surface(self.MEPS_np, self.AIP, self.Atom)

    def write_xml(self, filename):
        writer = AIPWriter(self.MEPS_np._cml_file)
        if self.dualAIP:
            dualAIPbool = [a.dual for a in self.AIP.AIP]
            writer.write_file(self.Surface, self.AIP.AIP, filename, dual=dualAIPbool)
        else:
            writer.write_file(self.Surface, self.AIP.AIP, filename)

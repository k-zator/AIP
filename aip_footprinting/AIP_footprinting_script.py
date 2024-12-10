"""
Parses foot-printing as class.
@author: Katarzyna Joanna Zator (kz265)
"""

import logging
from aip_footprinting.read_MEPS import MEPS
from aip_footprinting.footprinting import Footprinting
from aip_footprinting.atom_class import AtomSet
from aip_footprinting.surface_class import Surface
from filewriter.AIP_writer import AIPWriter

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.WARNING)


class KMC_footprint_script():

    # MEPS gets read and processed, surface quantified and divided
    def __init__(self, cml_file, cube_p_file, cube_m_file, cube_np_file, centre_surface_percentile=90, lp_excl_r=1.5, dualAIP=False):
        """Class for handling footprinting to find AIPs from command line. Parallel to aip_footprinting_ipy,
           but it handles exact file paths as supplied by parser. It has a write_xml function for writing of
           the XML file in the schema similar to the original SSIP code."""
        self.dualAIP = dualAIP
        self.MEPS_np = MEPS(cube_np_file, cml_file)
        try:
            self.MEPS_p = MEPS(cube_p_file, cml_file, own_dist_scaled=False)
            if not all(self.MEPS_np.Atoms_df == self.MEPS_p.Atoms_df):
                LOGGER.warn(
                    "\n The cube files contain different atom coordinates and this will cause errors")
        except Exception:
            self.MEPS_p = 0
            LOGGER.warn(
                "\n Could not localise the 0.0300 cube file, proceeding without this class object")
        try:
            self.MEPS_m = MEPS(cube_m_file, cml_file)
            if not all(self.MEPS_np.Atoms_df == self.MEPS_p.Atoms_df):
                LOGGER.warn(
                    "\n The cube files contain different atom coordinates and this will cause errors")
        except Exception:
            self.MEPS_m = 0
            LOGGER.warn(
                "\n Could not localise the 0.0104 cube file, proceeding without this class object")

        self.Atom = AtomSet()
        self.Atom._create_from_MEPS(self.MEPS_np)
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

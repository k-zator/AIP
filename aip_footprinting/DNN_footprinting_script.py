import logging
from aip_footprinting.read_MEPS import MEPS
from aip_footprinting.footprinting import Footprinting
from aip_footprinting.atom_class import AtomSet
from aip_footprinting.surface_class import Surface
from filewriter.AIP_writer import AIPWriter

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.WARNING)


class DNN_footprint_script():

    # MEPS gets read and processed, surface quantified and divided
    def __init__(self, cml_file, cube_file, centre_surface_percentile=90, lp_excl_r=2.0, dualAIP=False):
        """Class for handling footprinting to find AIPs from command line. Parallel to aip_footprinting_ipy,
           but it handles exact file paths as supplied by parser. It has a write_xml function for writing of
           the XML file in the schema similar to the original SSIP code."""
        self.MEPS_np = MEPS(cube_file, cml_file)
        self.MEPS_p = self.MEPS_np
        self.MEPS_m = self.MEPS_np
        self.Atom = AtomSet()
        self.Atom._create_from_MEPS(self.MEPS_np)
        self.AIP = Footprinting(self.MEPS_np, self.MEPS_p,
                                self.MEPS_m, self.Atom, centre_surface_percentile, lp_excl_r, dualAIP)
        self.Surface = Surface(self.MEPS_np, self.AIP, self.Atom)

    def write_xml(self, filename=False):
        writer = AIPWriter(self.MEPS_np._cml_file)
        if filename != False:
            writer.write_file(self.Surface, self.AIP.AIP, filename)
        else:
            writer.write_file(self.Surface, self.AIP.AIP,
                              writer.inchikey + "_dnn_aip.xml")

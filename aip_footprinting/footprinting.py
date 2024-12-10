"""
Class that carries out foot-printing calculation using all the MEP surfaces
reates instance of AIP class for each individual AIP. Useful for XML file writing.
@author: Katarzyna Joanna Zator (kz265) and Maria Chiara Storer (mcs92)
"""

from aip_footprinting.polar_search import polar_AIP_search
from aip_footprinting.non_polar_search import non_polar_AIP_search
from aip_footprinting.constants import subset_r_nn, Atom_based_exclusion
import logging
import numpy as np
import pandas as pd
from scipy import spatial
from scipy.linalg import expm, norm

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.INFO)


class Footprinting():
    """"This class footprints all the MEP surfaces of the ligand to calculate AIPs output
        as a list of AIP objects. Defines how to determine them for different atom types,
        which surface to use and which linear fit to use. """

    def __init__(self, MEPS_np, MEPS_p, MEPS_m, Atom, centre_surface_percentile=80, lp_excl_r=1.74, dualAIP=False):
        """
        Parameters
        ----------
        MEPS_np : class instance
            describes the non-polar (0.0020 eBohr-3) MEPS surface for non-polar areas

        MEPS_p : class instance
            describes the polar (0.0300 eBohr-3) MEPS surface for N and O lone pairs

        MEPS_m : class instance
            describes the middle (0.0104 eBohr-3) MEPS surface for H atoms

        centre_surface_percentile : float
            default 90, determines which fractional density of the atom-owned surface is considered.
            Created so that atom-owned surface edges are excluded from consideration

        dualAIP: bool
            default False, determines if the algorithm uses the cluster centre MEPS value, or the extreme one
        """
        self.csp = centre_surface_percentile
        self.lp_excl_r = lp_excl_r
        self._dualAIP = dualAIP
        self.MEPS = MEPS_np
        self.MEPS_p = MEPS_p
        self.MEPS_m = MEPS_m
        self._Atom = Atom

        self.AIP = []  # list of AIPs defined by their class

        self.surface_areas_initial = []
        self.surface_areas_remain = []
        self.sai_03 = []
        self.sai_01 = []

        # A loop running over all atoms of the ligand, adding AIPs for polar, hydrogens,
        # and non-polar areas from relevant MEPS
        for atom in range(self.MEPS.no_atoms):
            atomclass = self._Atom.Atom[atom]
            min_area = Atom_based_exclusion[atomclass.atom_type]
            self.FootprintByAtom(atomclass, min_area)

    def FootprintByAtom(self, atomclass, min_area):
        """An instance of the loop over atoms to assign relevant AIPs to the MEPS associated with the atom.
           A subsject of a set of definitions for variuos atom types, however, each is softened to also consider
           the effect of molecular geometry, and therefore MEPS availability"""

        # initial areas at different isosurfaces
        ia_002 = sum(self.MEPS.MEPS_owner == atomclass.index) * \
            self.MEPS.total_area/self.MEPS.no_MEPS
        if self.MEPS_p is not 0:
            ia_03 = sum(self.MEPS_p.MEPS_owner == atomclass.index) * \
                self.MEPS_p.total_area/self.MEPS_p.no_MEPS
        else:
            ia_03 = None
        if self.MEPS_m is not 0:
            ia_01 = sum(self.MEPS_m.MEPS_owner == atomclass.index) * \
                self.MEPS_m.total_area/self.MEPS_m.no_MEPS
        else:
            ia_01 = None
        self.surface_areas_initial.append(ia_002)
        self.sai_03.append(ia_03)
        self.sai_01.append(ia_01)

        if ia_002 > min_area:
            # define surface to search
            if atomclass.atom_type[0] == "H":
                MEPS = self.MEPS_m
            elif atomclass.atom_type == "F":
                MEPS = self.MEPS_p
            else:
                MEPS = self.MEPS

            MEPS_df_sub = self.EdgeDetection(MEPS, atomclass, self.csp)
            if atomclass.atom_type in ["S.3", "S.2.phene", "S.2.ps", "S.O", "Cl", "Br", "I"]:
                MEPS_df_sub = self.EdgeDetection(
                    self.MEPS_m, atomclass, self.csp)

            # only select atom types can be assigned lone pairs on the inner surface
            if atomclass.atom_type[:2] in ["N.", "O."] \
                    and "no_lp" not in atomclass.atom_type \
                    and "pl3" not in atomclass.atom_type \
                    or atomclass.atom_type == "S.2" \
                    or atomclass.atom_type == "S.2.ps":
                # locate the polar AIPs on inner surface
                MEPS_df_sub = polar_AIP_search(self.MEPS_p, self.AIP, self._Atom,
                                               atomclass, MEPS_df_sub, self.lp_excl_r)

                # if it failed to add AIPs due to the MEPS being too positive
                no_polar_AIPs = [s.atom_name for s in self.AIP].count(
                    atomclass.atom_name)
                if no_polar_AIPs == 0:
                    MEPS_df_sub = polar_AIP_search(self.MEPS, self.AIP, self._Atom,
                                                   atomclass, MEPS_df_sub, self.lp_excl_r, outer_meps=True)

            # update and check the surface area requirement again
            remainder_area = len(MEPS_df_sub) * MEPS.total_area/MEPS.no_MEPS
            if remainder_area > min_area * self.csp/100:
                non_polar_AIP_search(self.AIP, self._Atom, MEPS, self.MEPS_p,
                                     MEPS_df_sub, atomclass, self._dualAIP)

    @staticmethod
    def EdgeDetection(MEPS, atom, csp):
        """Detects MEPS points near atom-owned surface fragment edge based on local density of points.
        Takes centre_surface_percentile of densest points of the surface as atom-owned surface"""
        # the percentile cannot be 0, it ends up crashing the code
        if csp == 0:
            csp = 0.1
            LOGGER.warn("\n Centre surface percentile cannot be zero as it \
                            negates the point of the calculation, i.e. \
                            no MEPS is available to locate the AIPs. \
                            Using 0.1 instead")

        # matrix of distances; ndarray
        MEPS_df_atom = MEPS.MEPS_df[MEPS.MEPS_owner == atom.index]
        subset_mat_dis = spatial.distance.cdist(MEPS_df_atom[["x", "y", "z"]],
                                                MEPS_df_atom[["x", "y", "z"]])

        # non-averaged per MEPS point surface area; array
        subset_nonav_pMA = 1 / \
            (sum(subset_mat_dis < subset_r_nn)/(subset_r_nn**2*np.pi))
        edge_mask = subset_nonav_pMA <= np.percentile(subset_nonav_pMA, csp)
        MEPS_df = MEPS_df_atom.loc[edge_mask]

        return MEPS_df

"""
Class to describe 0.0020 e.Bohr-3  surface for the XML file output.
@author: Katarzyna Joanna Zator (kz265)
"""

from re import A
from aip_footprinting.linear_fit_aip import alpha_linear_002, beta_linear_all_002


class Surface():
    """"This class is to describe a molecule's 0.0020 surface."""

    def __init__(self, MEPS_np, AIP, Atom):
        self.total = MEPS_np.total_area
        self.numberOFMEPSPoints = len(MEPS_np.MEPS_df)
        self.positive = self.total * \
            sum(MEPS_np.MEPS_df["charge"] > 0)/self.numberOFMEPSPoints
        self.negative = self.total * \
            sum(MEPS_np.MEPS_df["charge"] < 0)/self.numberOFMEPSPoints
        self.electrostaticPotentialMax = max(MEPS_np.MEPS_df["charge"])
        self.electrostaticPotentialMin = min(MEPS_np.MEPS_df["charge"])
        self.isosurface = 0.0020
        self.VdWVolume = MEPS_np.vdw_volume
        self.get_surface_types(MEPS_np, AIP, Atom)

    def get_surface_types(self, MEPS_np, AIP, Atom):
        """determines the contributions to the total surface area from positive, negative, polar 
           and non-polar areas of the surface"""

        sa = AIP.surface_areas_initial  # for a list of areas total
        pp = []
        pn = []
        np = []
        nn = []
        for atom in Atom.Atom:
            MEP_points = MEPS_np.MEPS_df["charge"][MEPS_np.MEPS_owner ==
                                                   atom.index].values
            if len(MEP_points) > 0:
                fraction = sa[atom.index] / len(MEP_points)
                MEP_points_pos = MEP_points[MEP_points > 0]
                MEP_points_neg = MEP_points[MEP_points < 0]

                if len(MEP_points_pos) > 0:
                    if atom.atom_type in alpha_linear_002.keys():
                        a0, a1 = alpha_linear_002[atom.atom_type]
                    else:
                        a0, a1 = alpha_linear_002['H.soft']
                    AIP_values_pos = a0 + a1 * MEP_points_pos
                    pp.append(sum(AIP_values_pos > 1.5) * fraction)
                    pn.append(sum(AIP_values_pos <= 1.5) * fraction)

                if len(MEP_points_neg) > 0:
                    if atom.atom_type in beta_linear_all_002.keys():
                        b0, b1 = beta_linear_all_002[atom.atom_type]
                    else:
                        b0, b1 = beta_linear_all_002['C.ar']
                    AIP_values_neg = b0 + b1 * MEP_points_neg
                    # beta is a positive value
                    np.append(sum(AIP_values_neg > 2.5) * fraction)
                    nn.append(sum(AIP_values_neg <= 2.5) * fraction)

        self.positive_polar = sum(pp)
        self.negative_polar = sum(np)
        self.positive_non_polar = sum(pn)
        self.negative_non_polar = sum(nn)

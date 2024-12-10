from tables import Atom
from aip_footprinting.polar_search import polar_AIP_search
import pathlib
import unittest
import numpy as np
import logging
import pandas as pd
from aip_footprinting.AIP_class import AIP as AIPclass
from aip_footprinting.read_MEPS import MEPS
from aip_footprinting.atom_class import AtomSet, Atom
from aip_footprinting.define_AIP import define_extreme_geometric_AIP

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.INFO)


class TestDefineAIPPolar(unittest.TestCase):

    def setUp(self):
        parent_directory = pathlib.Path(__file__).resolve().parents[0]
        path = "test_files/select_surfaces/"

        self.N1_03 = pd.read_csv(
            "{}/{}benzonitrile_N_03.csv".format(parent_directory, path), index_col=0)
        #self.N2_03 = pd.read_csv("{}/{}_N_03.csv".format(parent_directory, path), index_col=0)
        self.Nar_03 = pd.read_csv(
            "{}/{}pyridine_N_03.csv".format(parent_directory, path), index_col=0)
        self.O2aldehyde_03 = pd.read_csv(
            "{}/{}acetal_O_03.csv".format(parent_directory, path), index_col=0)
        self.O2am_03 = pd.read_csv(
            "{}/{}acetamide_O_03.csv".format(parent_directory, path), index_col=0)
        self.O2carbonyl_03 = pd.read_csv(
            "{}/{}acetone_O_03.csv".format(parent_directory, path), index_col=0)
        self.O2nitro_03 = pd.read_csv(
            "{}/{}nitro_O_03.csv".format(parent_directory, path), index_col=0)
        self.O2noxide_03 = pd.read_csv(
            "{}/{}n_oxide_O_03.csv".format(parent_directory, path), index_col=0)
        #self.O2other_03 = pd.read_csv("{}/{}_O_03.csv".format(parent_directory, path), index_col=0)
        self.O2po_03 = pd.read_csv(
            "{}/{}po_O_03.csv".format(parent_directory, path), index_col=0)
        self.O2sulfone_03 = pd.read_csv(
            "{}/{}dms_O_03.csv".format(parent_directory, path), index_col=0)
        self.O2sulfoxide_03 = pd.read_csv(
            "{}/{}dmso_O_03.csv".format(parent_directory, path), index_col=0)
        self.O3alcohol_03 = pd.read_csv(
            "{}/{}ethanol_O_03.csv".format(parent_directory, path), index_col=0)
        self.O3any_03 = pd.read_csv(
            "{}/{}diethyl_ether_O_03.csv".format(parent_directory, path), index_col=0)
        self.O3water_03 = pd.read_csv(
            "{}/{}water_O_03.csv".format(parent_directory, path), index_col=0)
        self.S2_03 = pd.read_csv(
            "{}/{}thiocamphor_S_03.csv".format(parent_directory, path), index_col=0)
        self.S2ps_03 = pd.read_csv(
            "{}/{}ps_S_03.csv".format(parent_directory, path), index_col=0)
        self.S3_03 = pd.read_csv(
            "{}/{}thioether_S_03.csv".format(parent_directory, path), index_col=0)
        self.Spos_03 = pd.read_csv(
            "{}/{}positive_S_03.csv".format(parent_directory, path), index_col=0)
        self.N1_002 = pd.read_csv(
            "{}/{}benzonitrile_N_002.csv".format(parent_directory, path), index_col=0)
        self.Nar_002 = pd.read_csv(
            "{}/{}pyridine_N_002.csv".format(parent_directory, path), index_col=0)
        self.O2aldehyde_002 = pd.read_csv(
            "{}/{}acetal_O_002.csv".format(parent_directory, path), index_col=0)
        self.O2am_002 = pd.read_csv(
            "{}/{}acetamide_O_002.csv".format(parent_directory, path), index_col=0)
        self.O2carbonyl_002 = pd.read_csv(
            "{}/{}acetone_O_002.csv".format(parent_directory, path), index_col=0)
        self.O2nitro_002 = pd.read_csv(
            "{}/{}nitro_O_002.csv".format(parent_directory, path), index_col=0)
        self.O2sulfoxide_002 = pd.read_csv(
            "{}/{}dmso_O_002.csv".format(parent_directory, path), index_col=0)
        self.S2_002 = pd.read_csv(
            "{}/{}thiocamphor_S_002.csv".format(parent_directory, path), index_col=0)
        self.Spos_002 = pd.read_csv(
            "{}/{}positive_S_002.csv".format(parent_directory, path), index_col=0)
        """
        missing 002 files for polar-nonpolar handover:
        "N2": 0,
        "O2other":0,
        """

        self.MEPS_list = {"N1": self.N1_03,
                          "Nar": self.Nar_03,
                          "O2aldehyde": self.O2aldehyde_03,
                          "O2am": self.O2am_03,
                          "O2carbonyl": self.O2carbonyl_03,
                          "O2nitro": self.O2nitro_03,
                          "O2noxide": self.O2noxide_03,
                          # "O2other": self.O2other_03,
                          "O2po": self.O2po_03,
                          "O2sulfone": self.O2sulfone_03,
                          "O2sulfoxide": self.O2sulfoxide_03,
                          "O3alcohol": self.O3alcohol_03,
                          "O3any": self.O3any_03,
                          "O3water": self.O3water_03,
                          "S2": self.S2_03,
                          "S2ps": self.S2ps_03,
                          "S3": self.S3_03,
                          "N1_002": self.N1_002,
                          "Nar_002": self.Nar_002,
                          "O2aldehyde_002": self.O2aldehyde_002,
                          "O2am_002": self.O2am_002,
                          "O2carbonyl_002": self.O2carbonyl_002,
                          "O2nitro_002": self.O2nitro_002,
                          "O2sulfoxide_002": self.O2sulfoxide_002,
                          "S2_002": self.S2_002}

        self.twolps_atom_types = ["O.2.aldehyde", "O.2.am", "O.2.carbonyl", "O.2.nitro",
                                  "O.2.sulfoxide", "O.3.alcohol", "O.3.any", "O.3.water",
                                  "S.2", "S.3"]
        self.lpandpi_atom_types = ["N.1", "N.2", "N.ar","O.2.aldehyde", "O.2.am", "O.2.carbonyl",
                                    "O.2.nitro", "O.2.sulfoxide", "S.2"]
        self.threelps_atom_types = ["O.2.noxide", "O.2.po","O.2.sulfone", "S.2.ps"]

    def test_finding_both_lone_pairs(self):
        LOGGER.info(
            "Testing two correct MEP minima are designated as lone pairs")
        atom_indices_dict = {"O2aldehyde": [602, 1973],
                             "O2am": [958, 3141],
                             "O2carbonyl": [2473, 4588],
                             "O2nitro": [3953, 15673],
                             # "O2other": [],
                             "O2sulfoxide": [3284, 5515],
                             "O3alcohol": [937, 2330],
                             "O3any": [952, 2378],
                             "O3water": [468, 1751],
                             "S2": [21359, 22411],
                             "S3": [3823, 7689]}

        for i, (atype, indices) in enumerate(atom_indices_dict.items()):
            MEPS_class = MEPS
            MEPS_class.MEPS_df = self.MEPS_list[atype]
            MEPS_class.MEPS_owner = np.array([0] * len(MEPS_class.MEPS_df))
            AIP = []  # list of AIPs
            _Atom = AtomSet()
            atom = Atom()
            atom.atom_type = self.twolps_atom_types[i]
            atom.index = 0
            MEPS_np = pd.DataFrame(np.array([[0, 0, 0, 0]]), columns=[
                                   "x", "y", "z", "charge"])

            _ = polar_AIP_search(MEPS_class, AIP, _Atom, atom, MEPS_np)
            self.assertEqual(AIP[0].mepsindex, indices[0])
            self.assertEqual(AIP[1].mepsindex, indices[1])

    def test_three_lone_pairs(self):
        LOGGER.info(
            "Testing correct MEP points are designated as three lone pairs")

        atom_indices_xyz_dict = {"O2noxide": [[12299, 11699, 5072], [np.array([0.000326, -2.635925, 0.012248]), np.array([0.000115, -1.361599, 0.006343])]],
                                 "O2po": [[4183, 623, 4016], [np.array([-0.003018, 0.846169, 1.954379]), np.array([-0.001123, 0.249667, 0.576620])]],
                                 "O2sulfone": [[3363, 5435, 11298], [np.array([-0.345054, 1.561353, 1.173517]), np.array([-0.044259, 0.787797, -0.039875])]],
                                 "S2ps": [[26570, 19093, 12811], [np.array([-1.002335, 0.603370, 3.545830]), np.array([-0.511889, 0.010156, 1.720519])]]}

        for atype, values in atom_indices_xyz_dict.items():
            indices = values[0]
            xyz = values[1]
            MEPS_p = self.MEPS_list[atype]
            atom = Atom()
            atom.atom_type = atype
            atom.xyz = xyz[0]
            PS = Atom()
            PS.xyz = xyz[1]

            AIP = []  # list of AIPs
            aip0_xyz = MEPS_p.loc[indices[0]][[
                "x", "y", "z"]].to_numpy().reshape(1, 3)
            AIP += [AIPclass([-2.0, -0.1, aip0_xyz, indices[0],
                              "polar", 0, "O.2.po", 0])]

            AIP += define_extreme_geometric_AIP(MEPS_p, AIP[0], atom, PS, 3)
            centre_xyz = self.get_centre_of_circle(atom, PS, AIP[0])
            self.assertEqual(AIP[0].mepsindex, indices[0])
            self.assertEqual(AIP[1].mepsindex, indices[1])
            self.assertEqual(AIP[2].mepsindex, indices[2])

            self.assertAlmostEqual(np.degrees(self.get_angle(
                AIP[0].xyz, centre_xyz, AIP[1].xyz)), 120, places=-2)
            self.assertAlmostEqual(np.degrees(self.get_angle(
                AIP[0].xyz, centre_xyz, AIP[2].xyz)), 120, places=-1)

    def test_polar_to_nonpolar_handover(self):
        LOGGER.info(
            "Testing correct non-polar MEPS fraction was excised by lone pairs")
        atom_frac_dict = {"N1": 2096,
                          # "N2": 0,
                          "Nar": 737,
                          "O2aldehyde": 1902,
                          "O2am": 997,
                          "O2carbonyl": 568,
                          "O2nitro": 754,
                          "O2sulfoxide": 681,
                          "S2": 3711}

        for i, (atype, exp_fraction) in enumerate(atom_frac_dict.items()):
            MEPS_class = MEPS
            MEPS_class.MEPS_df = self.MEPS_list[atype]
            MEPS_class.MEPS_owner = np.array([0] * len(MEPS_class.MEPS_df))
            AIP = []  # list of AIPs
            _Atom = AtomSet()
            atom = Atom()
            atom.atom_type = self.lpandpi_atom_types[i]
            atom.index = 0
            MEPS_np = self.MEPS_list[atype+"_002"]

            MEPS_df_after = polar_AIP_search(
                MEPS_class, AIP, _Atom, atom, MEPS_np)
            self.assertEqual(len(MEPS_df_after), exp_fraction)

    @staticmethod
    def get_angle(a, b, c):
        ba = a - b
        bc = c - b
        cosine_angle = np.dot(ba, bc.T) / \
            (np.linalg.norm(ba) * np.linalg.norm(bc))
        angle = np.arccos(cosine_angle)
        return float(angle)

    def get_centre_of_circle(self, O, PS, first_AIP):
        norm_OPS = np.linalg.norm(O.xyz - PS.xyz)
        norm_OAIP = np.linalg.norm(O.xyz - first_AIP.xyz)
        PS_O_lp_angle = self.get_angle(PS.xyz, O.xyz, first_AIP.xyz)
        centre_xyz = O.xyz + \
            np.cos(np.pi - PS_O_lp_angle) * \
            (O.xyz - PS.xyz) * norm_OAIP / norm_OPS
        return centre_xyz


if __name__ == '__main__':
    unittest.main()

from unicodedata import decimal
from tables import Atom
from aip_footprinting.polar_search import polar_AIP_search
import pathlib
import unittest
import numpy as np
import logging
import pandas as pd
from aip_footprinting.AIP_class import AIP as AIPclass
from aip_footprinting.atom_class import AtomSet, Atom
from aip_footprinting.read_MEPS import MEPS
from aip_footprinting.polar_search import polar_AIP_search as sulfur_AIP_search
from aip_footprinting.define_AIP import get_AIP_value, define_extreme_AIP, \
    define_extreme_geometric_AIP, define_kmedoid_AIP_value_critical_point, \
    define_kmedoid_AIP_value_polar, define_single_cluster

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.INFO)


class TestDefineAIPNonPolar(unittest.TestCase):

    def setUp(self):
        parent_directory = pathlib.Path(__file__).resolve().parents[0]
        path = "test_files/select_surfaces/"
        self.F_03 = pd.read_csv(
            "{}/{}fluorobenzene_F_03.csv".format(parent_directory, path), index_col=0)
        self.S3_002 = pd.read_csv(
            "{}/{}mes2s_S3_002.csv".format(parent_directory, path), index_col=0)
        self.Br_002 = pd.read_csv(
            "{}/{}bromobenzene_Br_002.csv".format(parent_directory, path), index_col=0)
        self.Cl_002 = pd.read_csv(
            "{}/{}chlorobenzene_Cl_002.csv".format(parent_directory, path), index_col=0)
        self.C1_002 = pd.read_csv(
            "{}/{}benzonitrile_C1_002.csv".format(parent_directory, path), index_col=0)
        self.Car_002 = pd.read_csv(
            "{}/{}benzene_C_002.csv".format(parent_directory, path), index_col=0)
        self.Car_002_small = pd.read_csv(
            "{}/{}benzene_C_002_small.csv".format(parent_directory, path), index_col=0)
        self.Car_polarised_002 = pd.read_csv(
            "{}/{}pyridine_C_002.csv".format(parent_directory, path), index_col=0)
        self.C2_002 = pd.read_csv(
            "{}/{}/acetone_C2_002.csv".format(parent_directory, path), index_col=0)
        self.SO_002 = pd.read_csv(
            "{}/{}/dmso_SO_002.csv".format(parent_directory, path), index_col=0)
        self.I_002 = pd.read_csv(
            "{}/{}/iodobenzene_I_002.csv".format(parent_directory, path), index_col=0)
        self.N1_002 = pd.read_csv(
            "{}/{}/benzonitrile_N_002.csv".format(parent_directory, path), index_col=0)
        self.Nar_002 = pd.read_csv(
            "{}/{}/pyridine_N_002.csv".format(parent_directory, path), index_col=0)
        self.O2_pi_only_002 = pd.read_csv(
            "{}/{}O2_pi_only.csv".format(parent_directory, path), index_col=0)
        self.O2_pi_only_002_small = pd.read_csv(
            "{}/{}O2_pi_only_small.csv".format(parent_directory, path), index_col=0)
        self.paracyclophane_Car_002 = pd.read_csv(
            "{}/{}paracyclophane_C_002.csv".format(parent_directory, path), index_col=0)
        self.MEPS_list = {"Br": self.Br_002,
                          "Cl": self.Cl_002,
                          "I": self.I_002,
                          "C1": self.C1_002,
                          "Car": self.Car_002,
                          "Car_polarised": self.Car_polarised_002,
                          "C2": self.C2_002,
                          "S.O": self.SO_002,
                          "I": self.I_002,
                          "N1": self.N1_002,
                          "Nar": self.Nar_002,
                          "O2_pi_only": self.O2_pi_only_002}

    def test_halide_sp1_AIP_extreme(self):
        LOGGER.info(
            "Testing correct extremum is chosen for geometric halide/sp1 AIPs")

        atom_indices_value_dict = {"Cl": [12485, -0.024471],
                                   "Br": [21929, -0.023515],
                                   "I": [8439, -0.020589],
                                   "C1": [24194, -0.028194],
                                   "N1": [17237, -0.022652]}

        for atype, values in atom_indices_value_dict.items():
            index = values[0]
            mepsvalue = values[1]
            atom = Atom()
            atom.atom_type = atype
            AIP = []  # list of AIPs
            if atype == "N1":
                AIP += [define_extreme_AIP(self.MEPS_list[atype],
                                           atom, minimum=False, polar=False)]
            else:
                AIP += [define_extreme_AIP(self.MEPS_list[atype],
                                           atom, minimum=True, polar=False)]

            self.assertEqual(AIP[0].mepsindex, index)
            self.assertAlmostEqual(AIP[0].mepsvalue, mepsvalue)

    def test_halide_sp1_AIP_distribution(self):
        LOGGER.info(
            "Testing correct MEP points are designated as geometric halide/sp1 AIPs")

        atom_indices_xyz_dict = {"Cl": [[12485, 9623, 17000, 8208], [np.array([0.072846, 3.084815, 0.028936]), np.array([0.031117, 1.324146, 0.012199])]],
                                 "Br": [[21929, 23333, 8221, 6463], [np.array([0.076438, 3.232995, 0.029999]), np.array([0.030990, 1.309352, 0.011969])]],
                                 "I": [[8439, 4313, 15741, 28498], [np.array([3.423660, -0.281233, -0.117434]), np.array([1.291790, -0.106294, -0.044508])]],
                                 "C1": [[24194, 7441, 8827, 23898], [np.array([1.720099, 1.244253, 0.095085]), np.array([0.432303, 0.613914, 0.044704])]],
                                 "N1": [[17237, 13582, 24513, 13571], [np.array([2.764005, 1.755980, 0.135970]), np.array([1.720099, 1.244253, 0.095085])]]}

        for atype, values in atom_indices_xyz_dict.items():
            indices = values[0]
            xyz = values[1]
            atom = Atom()
            atom.atom_type = atype
            atom.xyz = xyz[0]
            PS = Atom()
            PS.xyz = xyz[1]

            AIP = []  # list of AIPs
            MEPS_np = self.MEPS_list[atype]

            aip0_xyz = MEPS_np.loc[indices[0]][[
                "x", "y", "z"]].to_numpy().reshape(1, 3)
            AIP += [AIPclass(
                [-2.0, -0.1, aip0_xyz, indices[0], "non_polar", 0, "Cl", 0])]
            AIP += define_extreme_geometric_AIP(
                MEPS_np, AIP[0], atom, PS, 4, polar=False)
            centre_xyz = self.get_centre_of_circle(
                atom, PS, AIP[0])
            self.assertAlmostEqual(np.degrees(
                self.get_angle(AIP[0].xyz, centre_xyz, AIP[1].xyz)), 90, places=-1)
            self.assertAlmostEqual(np.degrees(
                self.get_angle(AIP[0].xyz, centre_xyz, AIP[2].xyz)), 180, places=-1)
            self.assertAlmostEqual(np.degrees(
                self.get_angle(AIP[0].xyz, centre_xyz, AIP[3].xyz)), 90, places=-1)
            
            self.assertEqual(AIP[1].mepsindex, indices[1])
            self.assertEqual(AIP[2].mepsindex, indices[2])
            self.assertEqual(AIP[3].mepsindex, indices[3])

    def test_kmedoid_AIPs(self):
        LOGGER.info("Testing correct choice of medoid points")
        atom_indices_dict = {"Car": [4707, 13529], "Car_polarised": [17847, 17874],
                             "C2": [6294, 5194], "O2_pi_only": [6805, 13008]}

        for atype, indices in atom_indices_dict.items():
            atom = Atom()
            atom.atom_type = "C.ar"
            atom.vdW_radius = 1.7

            MEPS_np = self.MEPS_list[atype]
            AIP = []
            if atype == "O2_pi_only":
                AIP += define_kmedoid_AIP_value_polar(MEPS_np, atom, sigma=False)
            else:
                AIP += define_kmedoid_AIP_value_critical_point(
                    MEPS_np, atom, minimum=True, sigma=False)

            self.assertEqual(AIP[0].mepsindex, indices[0])
            self.assertEqual(AIP[1].mepsindex, indices[1])

    def test_sulfoxide_S_AIP(self):
        LOGGER.info(
            "Testing assignment of AIP of sulfoxide sulfur on outer MEPS")
        atom = Atom()
        atom.atom_type = "S.O"
        AIP = []
        AIP += [define_single_cluster(self.SO_002, atom)]
        self.assertEqual(AIP[0].mepsindex, 4639)

    def test_S3_AIP(self): 
            atom = Atom()
            atom.atom_type = "S.3"
            atom.index = 0
            MEPS_m = pd.DataFrame(np.array([[0, 0, 0, 0]]), columns=["x", "y", "z", "charge"])
            MEPS_class = MEPS
            MEPS_class.MEPS_df = self.S3_002
            MEPS_class.MEPS_owner = np.array([0] * len(MEPS_class.MEPS_df))
            AIP = []
            _ = sulfur_AIP_search(MEPS_class, AIP, Atom, atom, MEPS_m, outer_meps=True)
            self.assertEqual(AIP[0].mepsindex, 13971)
            self.assertEqual(AIP[1].mepsindex, 15239)

    def test_sa_edge_cases_polar(self):
        LOGGER.info(
            "Testing medoid adjusts based on cluster size")
        atom = Atom()
        atom.atom_type = "O.2.carbonyl"
        atom.vdW_radius = 1.7
        AIP = []
        AIP += define_kmedoid_AIP_value_polar(self.O2_pi_only_002, atom, sigma=False)
        self.assertEqual(len(AIP), 2)

        AIP = []
        AIP += define_kmedoid_AIP_value_polar(self.O2_pi_only_002_small, atom, sigma=False)
        self.assertEqual(len(AIP), 1)

    def test_sa_edge_cases_nonpolar(self):
        atom = Atom()
        atom.atom_type = "C.ar"
        atom.vdW_radius = 1.7
        AIP = []
        AIP += define_kmedoid_AIP_value_critical_point(self.Car_002, atom, sigma=False)
        self.assertEqual(len(AIP), 2)

        AIP = []
        AIP += define_kmedoid_AIP_value_critical_point(
            self.Car_002_small, atom, sigma=False)
        self.assertEqual(len(AIP), 1)

    def test_sa_edge_unequal(self):
        LOGGER.info(
            "Testing medoid adjusts based on uneven cluster sizes")

        atom = Atom()
        atom.atom_type = "C.ar"
        atom.vdW_radius = 1.7
        AIP = []
        AIP += define_kmedoid_AIP_value_critical_point(
            self.paracyclophane_Car_002, atom, sigma=False)
        self.assertEqual(len(AIP), 1)

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

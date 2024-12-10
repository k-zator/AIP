import pathlib
import unittest
import numpy as np
import logging
from scipy.linalg import expm, norm
from aip_footprinting.read_MEPS import MEPS
from aip_footprinting.atom_class import AtomSet
from aip_footprinting.footprinting import Footprinting
from aip_footprinting.AIP_class import AIP as AIPclass
from aip_footprinting.surface_class import Surface

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.INFO)


class TestAIPWhole(unittest.TestCase):
    """Test AIP Parse class' handling of AIP identification"""

    def setUp(self):
        """Set up instances before testing"""
        centre_surface_percentile = 70
        parent_directory = pathlib.Path(__file__).resolve().parents[0]
        path = parent_directory / "test_files"
        water = "XLYOFNOQVPJJNP-UHFFFAOYSA-N"
        pyridine = "JUJWROOIHBZHMG-UHFFFAOYSA-N"
        mes2s = "QMMFVYPAHWMCMS-UHFFFAOYSA-N"
        cube_inner_file = (
            (path / water / "XLYOFNOQVPJJNP-UHFFFAOYSA-N_0.0300_merged.cube").absolute().as_posix())
        cube_middle_file = (
            (path / water / "XLYOFNOQVPJJNP-UHFFFAOYSA-N_0.0104_merged.cube").absolute().as_posix())
        cube_outer_file = (
            (path / water / "XLYOFNOQVPJJNP-UHFFFAOYSA-N_0.0020_merged.cube").absolute().as_posix())
        cml_file = (
            (path / water / "XLYOFNOQVPJJNP-UHFFFAOYSA-N.cml").absolute().as_posix())
        self.MEPS_p = MEPS(cube_inner_file, cml_file)
        self.MEPS_np = MEPS(cube_outer_file, cml_file)
        self.MEPS_m = MEPS(cube_middle_file, cml_file)
        self.Atom = AtomSet()
        self.Atom._create_from_MEPS(self.MEPS_np)
        self.AIP = Footprinting(self.MEPS_np, self.MEPS_p,
                                self.MEPS_m, self.Atom, 80)

        cube_inner_file_2 = (
            (path / pyridine / "JUJWROOIHBZHMG-UHFFFAOYSA-N_0.0300_merged.cube").absolute().as_posix())
        cube_middle_file_2 = (
            (path / pyridine / "JUJWROOIHBZHMG-UHFFFAOYSA-N_0.0104_merged.cube").absolute().as_posix())
        cube_outer_file_2 = (
            (path / pyridine / "JUJWROOIHBZHMG-UHFFFAOYSA-N_0.0020_merged.cube").absolute().as_posix())
        cml_file_2 = (
            (path / pyridine / "JUJWROOIHBZHMG-UHFFFAOYSA-N.cml").absolute().as_posix())
        self.MEPS_p_pyr = MEPS(cube_inner_file_2, cml_file_2)
        self.MEPS_np_pyr = MEPS(cube_outer_file_2, cml_file_2)
        self.MEPS_m_pyr = MEPS(cube_middle_file_2, cml_file_2)
        self.Atom_pyr = AtomSet()
        self.Atom_pyr._create_from_MEPS(self.MEPS_np_pyr)
        self.AIP_pyr = Footprinting(
            self.MEPS_np_pyr, self.MEPS_p_pyr, self.MEPS_m_pyr, self.Atom_pyr, 80)
        
        cube_inner_file_3 = (
            (path / mes2s / "QMMFVYPAHWMCMS-UHFFFAOYSA-N_0.0300_merged.cube").absolute().as_posix())
        cube_middle_file_3 = (
            (path / mes2s / "QMMFVYPAHWMCMS-UHFFFAOYSA-N_0.0104_merged.cube").absolute().as_posix())
        cube_outer_file_3 = (
            (path / mes2s / "QMMFVYPAHWMCMS-UHFFFAOYSA-N_0.0020_merged.cube").absolute().as_posix())
        cml_file_3 = (
            (path / mes2s / "QMMFVYPAHWMCMS-UHFFFAOYSA-N.cml").absolute().as_posix())
        self.MEPS_p_mes = MEPS(cube_inner_file_3, cml_file_3)
        self.MEPS_np_mes = MEPS(cube_outer_file_3, cml_file_3)
        self.MEPS_m_mes = MEPS(cube_middle_file_3, cml_file_3)
        self.Atom_mes = AtomSet()
        self.Atom_mes._create_from_MEPS(self.MEPS_np_mes)
        self.AIP_mes = Footprinting(
            self.MEPS_np_mes, self.MEPS_p_mes, self.MEPS_m_mes, self.Atom_mes, 80)

    def tearDown(self):
        """Clean up after tests"""
        del self.MEPS_p
        del self.MEPS_m
        del self.MEPS_np
        del self.Atom
        del self.AIP

        del self.MEPS_p_pyr
        del self.MEPS_m_pyr
        del self.MEPS_np_pyr
        del self.Atom_pyr
        del self.AIP_pyr

        del self.MEPS_p_mes
        del self.MEPS_m_mes
        del self.MEPS_np_mes
        del self.Atom_mes
        del self.AIP_mes

    def test_hydrogen_surface_assignment(self):
        LOGGER.info("Testing no negative MEP polar-hydrogens exist")
        structure_Hs = self.MEPS_m.Atoms_df[
            self.MEPS_m.Atoms_df["amass"] == 1].index.tolist()
        MEPS_Hs_mask = [
            True if i in structure_Hs else False for i in self.MEPS_m.MEPS_owner]
        MEPS_Hs_min_value = min(self.MEPS_m.MEPS_df[MEPS_Hs_mask]["charge"])
        self.assertTrue(MEPS_Hs_min_value > 0)

    def test_hydrogen_surface_pyr_assignment(self):
        LOGGER.info("Testing no negative MEP non-polar-hydrogens exist")
        structure_Hs = self.MEPS_m_pyr.Atoms_df[self.MEPS_m_pyr.Atoms_df["amass"] == 1].index.tolist()
        MEPS_Hs_mask = [
            True if i in structure_Hs else False for i in self.MEPS_m_pyr.MEPS_owner]
        MEPS_Hs_min_value = min(
            self.MEPS_m_pyr.MEPS_df[MEPS_Hs_mask]["charge"])
        self.assertTrue(MEPS_Hs_min_value > 0)

    def test_AIP_positions(self):
        LOGGER.info("Testing calculation of AIPs with all its attributes")
        AIP_values = [-4.5, -4.5, 2.8, 2.8]
        AIP_MEPS_values = [-0.088262, -0.088128, 0.164896, 0.164782]
        AIP_positions = [[-0.483867277729, -0.866746945955, 0.915166112278],
                         [-0.521142505609,  0.7897527508090001, 0.994967062232],
                         [1.504737554111, 0, 0.088252964998],
                         [-0.856759259257,  0.015526582357, -1.245078867043]]
        AIP_index = [468, 1751, 9381, 5112]
        AIP_owner = [0, 0, 1, 2]
        frac_values = [1, 1, 1, 1]

        self.assertSequenceEqual([a.value for a in self.AIP.AIP], AIP_values)
        self.assertSequenceEqual([a.mepsvalue for a in self.AIP.AIP], AIP_MEPS_values)
        self.assertSequenceEqual([list(a.xyz[0]) for a in self.AIP.AIP], AIP_positions)
        self.assertSequenceEqual([a.mepsindex for a in self.AIP.AIP], AIP_index)
        self.assertSequenceEqual([a.atom_owner_index for a in self.AIP.AIP], AIP_owner)
        self.assertSequenceEqual([a.fraction for a in self.AIP.AIP], frac_values)


    def test_AIP_positions_pyr(self):
        AIP_values = [-0.83, -0.82, -0.85, -0.85, -0.85, -0.85, -1.88, -1.92, -1.94, -1.8, -7.23, -3.01, -2.85, 1.24, 1.19, 1.18, 0.94, 0.95]
        AIP_MEPS_values = [-0.01281, -0.012646, -0.013105, -0.013113, -0.013108, -0.013112, -0.028879, -0.029536, -0.029827, 
                           -0.027704, -0.094217, -0.046377, -0.04387, 0.091682, 0.089079, 0.088842, 0.07695, 0.077224]
        AIP_positions = [[-0.117902223131, 1.7491786981900002, -1.669902162643],
                         [0.053016126922, 1.847858566796, 1.607251839259],
                         [-1.81553749399, 0.723015593454, -1.5981351779030002],
                         [-1.68385179754, 0.792587022821, 1.711041440977],
                         [1.620257951565, 0.829111877715, -1.7323027144829999],
                         [1.738961348674, 0.7697408642, 1.607251839259],
                         [1.6269954331290002, -1.076385177098, 1.630850487574],
                         [1.4766028005520002, -1.190027525379, -1.714326571793],
                         [-1.634996060192, -1.3392125760420002, -1.551289783978],
                         [-1.486323252865, -1.333451426043, 1.6999477743490001],
                         [-0.066148712531, -2.725554184881, 0.101359620934],
                         [-0.00729735083, -2.193401093433, -1.59528873482],
                         [0.07074778983800001, -1.996280544225, 1.692370488886],
                         [-0.004860490745, 3.106286452841, -0.061059617322],
                         [-2.888835981647, 1.3650824520409999, 0.101359620934],
                         [2.873411001274, 1.402833939221, -0.0946697653],
                         [2.759348486455, -1.938511878666, -0.06914861694400001],
                         [-2.7632035409, -1.924882925031, 0.107733028722]]
        AIP_index = [1240, 4467, 18022, 6422, 13119, 1778, 15560, 3161, 17874, 14794, 1815, 12036, 10064, 38769, 11673, 23808, 35221, 35860]
        AIP_owner = [0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 5, 6, 7, 8, 9, 10]
        frac_values = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1. , 0.5, 0.5, 1. , 1. , 1. , 1. , 1. ]
        self.assertSequenceEqual([a.value for a in self.AIP_pyr.AIP], AIP_values)
        self.assertSequenceEqual([a.mepsvalue for a in self.AIP_pyr.AIP], AIP_MEPS_values)
        self.assertSequenceEqual([list(a.xyz[0]) for a in self.AIP_pyr.AIP], AIP_positions)
        self.assertSequenceEqual([a.mepsindex for a in self.AIP_pyr.AIP], AIP_index)
        self.assertSequenceEqual([a.atom_owner_index for a in self.AIP_pyr.AIP], AIP_owner)
        self.assertSequenceEqual([a.fraction for a in self.AIP_pyr.AIP], frac_values)

    def test_Edge_Detection(self):
        LOGGER.info("Testing EdgeDetection removes right fraction")
        ox_edge_MEPS = self.AIP.EdgeDetection(
            self.AIP.MEPS, self.Atom.Atom[0], 70)
        self.assertEqual(len(ox_edge_MEPS), 2174)

    def test_Edge_Detection_pyr(self):
        ox_edge_MEPS = self.AIP_pyr.EdgeDetection(
            self.AIP_pyr.MEPS,  self.Atom.Atom[0], 70)
        self.assertEqual(len(ox_edge_MEPS), 1647)

    
    def test_vector_AIP_alignment(self):
        h_3_xyz = np.array([-2.29399605,  0.67886265, -0.05993935]) #connected to c_1
        c_1_xyz = np.array([-1.39287271,  0.06327846,  0.00451758])
        c_2_xyz = np.array([ 1.3940078,  -0.02719599, -0.00971992])
        s_0_xyz = np.array([ 0.03782663,  1.18773672, -0.14851141])

        h3_df = self.MEPS_m_mes.MEPS_df[self.MEPS_m_mes.MEPS_owner == 3]
        s0_df = self.MEPS_m_mes.MEPS_df[self.MEPS_m_mes.MEPS_owner == 0]

        h3_AIP = 21868
        s0_AIP = [30029, 14152] #sigma holes not lone pairs

        #for hydrogen
        dist = [norm(np.cross(h_3_xyz-p, c_1_xyz-h_3_xyz))/norm(c_1_xyz-h_3_xyz) for p in h3_df[["x", "y", "z"]].to_numpy()]
        AIP_index = h3_df.index[np.argmin(dist)]
        self.assertEqual(AIP_index, h3_AIP)

        #for sulfur
        dist1 = [norm(np.cross(s_0_xyz-p, c_1_xyz-s_0_xyz))/norm(c_1_xyz-s_0_xyz) for p in s0_df[["x", "y", "z"]].to_numpy()]
        AIP_index1= s0_df.index[np.argmin(dist1)]
        self.assertEqual(AIP_index1, s0_AIP[0])
        dist2 = [norm(np.cross(s_0_xyz-p, c_2_xyz-s_0_xyz))/norm(c_2_xyz-s_0_xyz) for p in s0_df[["x", "y", "z"]].to_numpy()]
        AIP_index2= s0_df.index[np.argmin(dist2)]
        self.assertEqual(AIP_index2, s0_AIP[1])

if __name__ == '__main__':
    unittest.main()

import logging
import pathlib
import unittest
import numpy as np
import pandas as pd
from aip_footprinting.constants import bohr_to_Angstrom
from aip_footprinting.read_MEPS import MEPS

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.WARN)

class TestCubeParse(unittest.TestCase):
    """Test .cube file parsing and handling"""

    def setUp(self):
        """Set up instances before testing"""
        parent_directory = pathlib.Path(__file__).resolve().parents[0]
        cube_inner_file = ((parent_directory / "test_files/test_inner.cube").absolute().as_posix())
        cube_outer_file = ((parent_directory / "test_files/test_outer.cube").absolute().as_posix())        
        cml_file = ((parent_directory / "test_files/test.cml").absolute().as_posix())
        self.MEPS_p = MEPS(cube_inner_file, cml_file, own_dist_scaled=False)
        self.MEPS_np = MEPS(cube_outer_file, cml_file)

    def tearDown(self):
        """Clean up after tests"""
        del self.MEPS_p
        del self.MEPS_np

    def test_area(self):
        total_area_np = 37.249475717728636
        total_area_p = 16.594608547199424
        LOGGER.info("Testing area calcualtion")
        self.assertAlmostEqual(self.MEPS_np.total_area, total_area_np)
        self.assertAlmostEqual(self.MEPS_p.total_area, total_area_p)

    def test_atoms_df(self):
        LOGGER.info("Testing parsing of Atoms from .cube files into Atoms_df")
        no_atoms = 3
        self.assertEqual(self.MEPS_np.no_atoms, no_atoms)
        self.assertEqual(self.MEPS_p.no_atoms, no_atoms)

        atoms_df = pd.DataFrame([[8.0,0.0,-0.197756,0.0,0.346603,"O.3.water","a1"],
                                 [1.0,0.0,0.760547,0.0,0.204588,"H.O","a2"],
                                 [1.0,0.0,-0.561792,0.0,-0.551191,"H.O","a3"]],
                                 columns=['amass','charge','x','y','z','atype','aname'])

        pd.testing.assert_frame_equal(self.MEPS_np.Atoms_df, atoms_df)
        self.assertDictEqual(self.MEPS_np.Atoms_df.round(6).to_dict(), atoms_df.to_dict())
        pd.testing.assert_frame_equal(self.MEPS_p.Atoms_df, atoms_df)
        self.assertDictEqual(self.MEPS_p.Atoms_df.round(6).to_dict(), atoms_df.to_dict())

    
    def test_atom_type(self):
        LOGGER.info("Testing atom types")
        atom_types = ["O.3.water", "H.O", "H.O"]
        self.assertSequenceEqual(list(self.MEPS_np.Atoms_df.atype), atom_types)

    def test_meps_df(self):
        LOGGER.info("Testing parsing of MEPS.cube information in class")
        no_meps_np = 5128
        no_meps_p = 2592
        self.assertEqual(self.MEPS_np.no_MEPS, no_meps_np)
        self.assertEqual(self.MEPS_p.no_MEPS, no_meps_p)

        first_series_np = [-1.738829, -0.327869, 0.128389, -0.010121]
        first_series_p = [-1.175471, -0.240000, 0.436361,  0.004092]
        self.assertSequenceEqual(self.MEPS_np.MEPS_df.iloc[0,:4].round(6).tolist(), first_series_np)
        self.assertSequenceEqual(self.MEPS_p.MEPS_df.iloc[0,:4].round(6).tolist(), first_series_p)

    def test_indices_reset(self):
        LOGGER.info("Testing indices of MEPS_df")
        true_index = self.MEPS_np.MEPS_df.index.to_list()
        sorted_index = true_index.copy()
        sorted_index.sort()
        self.assertEqual(true_index, sorted_index)


if __name__ == '__main__':
    unittest.main()




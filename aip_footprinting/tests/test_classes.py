import pathlib
import unittest
import numpy as np
import logging
from aip_footprinting.constants import bohr_to_Angstrom
from aip_footprinting.read_MEPS import MEPS
from aip_footprinting.atom_class import AtomSet, Atom
from aip_footprinting.footprinting import Footprinting
from aip_footprinting.AIP_class import AIP as AIPclass
from aip_footprinting.surface_class import Surface

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.INFO)


class TestClasses(unittest.TestCase):
    """Test AIP Parse class' handling of AIP identification"""

    def setUp(self):
        """Set up instances before testing"""
        parent_directory = pathlib.Path(__file__).resolve().parents[0]
        path = parent_directory / "test_files"
        water = "XLYOFNOQVPJJNP-UHFFFAOYSA-N"
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
                                self.MEPS_m, self.Atom, centre_surface_percentile=70)
        self.Surface = Surface(self.MEPS_np, self.AIP, self.Atom)

    def tearDown(self) -> None:
        del self.MEPS_p
        del self.MEPS_np
        del self.MEPS_m
        del self.Atom
        del self.AIP
        del self.Surface

    def test_instances(self):
        LOGGER.info("Testing class instances")
        self.assertIsInstance(self.Atom, AtomSet)
        self.assertIsInstance(self.Surface, Surface)
        self.assertIsInstance(self.MEPS_np, MEPS)
        self.assertIsInstance(self.MEPS_m, MEPS)
        self.assertIsInstance(self.MEPS_p, MEPS)
        self.assertIsInstance(self.AIP, Footprinting)

    def test_AIP_class(self):
        LOGGER.info("Testing AIP class")
        aip = AIPclass()
        aip.set_value(-5)
        aip.set_nearest_atom("a0")
        aip.set_atomindex(1)
        aip.set_xyz(0, 1, 2)
        aip.set_atom_type("C.ar")
        aip.set_mepsvalue(-0.1)
        aip.set_mepsindex(100)
        aip.set_isosurface(0.0020)
        aip.set_AIP_area_fraction(0.5)

        aip2 = AIPclass([-5, -0.1, np.array([[0, 1, 2]]), 100, "non-polar",
                         1, "C.ar", "a0"])

        for a in [aip, aip2]:
            self.assertEqual(a.value, -5)
            self.assertEqual(a.atom_owner_index, 1)
            self.assertSequenceEqual(a.xyz.tolist(), [[0, 1, 2]])
            self.assertEqual(a.atom_type, "C.ar")
            self.assertEqual(a.atom_name, "a0")
            self.assertEqual(a.mepsvalue, -0.1)
            self.assertEqual(a.mepsindex, 100)
            self.assertEqual(a.isosurface, 0.0020)
            self.assertEqual(a.fraction, 0.5)

    def test_AIP_class_water(self):
        self.assertEqual(len(self.AIP.AIP), 4)
        self.assertIsInstance(self.AIP.AIP[0], AIPclass)
        self.assertEqual(self.AIP.AIP[0].value, -4.5)
        self.assertEqual(self.AIP.AIP[0].mepsvalue, -0.088262)
        self.assertSequenceEqual(self.AIP.AIP[0].xyz.tolist(), [
                                 [-0.483867277729, -0.866746945955, 0.915166112278]])
        self.assertEqual(self.AIP.AIP[0].mepsindex, 468)
        self.assertEqual(self.AIP.AIP[0].type, "polar")
        self.assertEqual(self.AIP.AIP[0].atom_owner_index, 0)
        self.assertEqual(self.AIP.AIP[0].atom_type, "O.3.water")
        self.assertEqual(self.AIP.AIP[0].atom_name, "a1")
        self.assertEqual(self.AIP.AIP[0].fraction, 1)
        self.assertEqual(self.AIP.AIP[0].isosurface, 0.0300)

    def test_Atom_class(self):
        LOGGER.info("Testing Atom class")
        atom = Atom()
        atom_row = self.MEPS_np.Atoms_df.loc[0]
        atom._create_atom_from_df(0, atom_row)
        self.assertEqual(atom.index, 0)
        self.assertEqual(atom.atomic_number, 8)
        self.assertSequenceEqual(
            atom.xyz.tolist(), [-0.19775556160799998, -0.0, 0.34660299734500005])
        self.assertEqual(atom.atom_type, "O.3.water")
        self.assertEqual(atom.atom_name, "a1")
        self.assertEqual(atom.covalent_radius, 0.63)
        self.assertEqual(atom.vdW_radius, 1.52)

    def test_AtomSet_class(self):
        LOGGER.info("Testing AtomSet class")
        xyz_coordinates = np.array([[-0.373704, 0.000000, 0.654985],
                                    [1.437226, 0.000000, 0.386615],
                                    [-1.061633, 0.000000, -1.041601]]) * bohr_to_Angstrom
        atom_types = ["O.3.water", "H.O.water", "H.O.water"]
        self.assertSequenceEqual(
            self.Atom.xyz.tolist(), xyz_coordinates.tolist())
        self.assertSequenceEqual(
            self.Atom.atom_type.tolist(), atom_types)

    def test_Surface_class(self):
        LOGGER.info("Testing Surface class")
        self.assertEqual(self.Surface.total, 37.249475717728636)
        self.assertEqual(self.Surface.numberOFMEPSPoints, 5128)
        self.assertEqual(self.Surface.positive, 19.19858898634103)
        self.assertEqual(self.Surface.negative, 18.05088673138761)
        self.assertEqual(self.Surface.electrostaticPotentialMax, 0.084445)
        self.assertEqual(self.Surface.electrostaticPotentialMin, -0.070811)
        self.assertEqual(self.Surface.isosurface, 0.0020)

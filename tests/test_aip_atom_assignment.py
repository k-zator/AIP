import unittest
from filereader.cml_reader import CmlReader
from filereader.mol2reader import Mol2Reader
from aip_atom_types.assign_aip_atom_types import assign_aip_atom_types
import os
from lxml import etree
import logging
import networkx as nx

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.INFO)

FIXTURE_DIR = os.path.join(os.path.dirname(
    os.path.realpath(__file__)), 'test_files', '')

cml_water = f"{FIXTURE_DIR}/water.cml"
mol2_water = f"{FIXTURE_DIR}/water.mol2"

cml_phenol = f"{FIXTURE_DIR}/phenol.cml"
mol2_phenol = f"{FIXTURE_DIR}/phenol.mol2"

cml_methyl_pyrazole = f"{FIXTURE_DIR}/methyl_pyrazole.cml"
mol2_methyl_pyrazole = f"{FIXTURE_DIR}/methyl_pyrazole.mol2"

cml_thiazolone = f"{FIXTURE_DIR}/thiazolone.cml"
mol2_thiazolone = f"{FIXTURE_DIR}/thiazolone.mol2"

cml_n_nitromethanamine = f"{FIXTURE_DIR}/n-nitromethanamine.cml"
mol2_n_nitromethanamine = f"{FIXTURE_DIR}/n-nitromethanamine.mol2"

cml_urea = f"{FIXTURE_DIR}/urea.cml"
mol2_urea = f"{FIXTURE_DIR}/urea.mol2"

cml_methylthiourea = f"{FIXTURE_DIR}/methylthiourea.cml"
mol2_methylthiourea = f"{FIXTURE_DIR}/methylthiourea.mol2"

cml_aniline_like = f"{FIXTURE_DIR}/aniline_like.cml"
mol2_aniline_like = f"{FIXTURE_DIR}/aniline_like.mol2"

cml_dmso = f"{FIXTURE_DIR}/dmso.cml"
mol2_dmso = f"{FIXTURE_DIR}/dmso.mol2"

cml_phosphonate = f"{FIXTURE_DIR}/phosphonate.cml"
mol2_phosphonate = f"{FIXTURE_DIR}/phosphonate.mol2"

cml_sulfonate = f"{FIXTURE_DIR}/sulfonate.cml"
mol2_sulfonate = f"{FIXTURE_DIR}/sulfonate.mol2"

cml_acetic_acid = f"{FIXTURE_DIR}/acetic_acid.cml"
mol2_acetic_acid = f"{FIXTURE_DIR}/acetic_acid.mol2"

cml_tenoxicam = f"{FIXTURE_DIR}/tenoxicam.cml"
mol2_tenoxicam = f"{FIXTURE_DIR}/tenoxicam.mol2"

cml_thiadiazole = f"{FIXTURE_DIR}/thiadiazole.cml"
mol2_thiadiazole = f"{FIXTURE_DIR}/thiadiazole.mol2"

cml_haloethane = f"{FIXTURE_DIR}/haloethane.cml"
mol2_haloethane = f"{FIXTURE_DIR}/haloethane.mol2"


class TestCmlMol2Comparison(unittest.TestCase):
    def setUp(self):
        self.cml_water = CmlReader(cml_water, ns=None)
        self.mol2_water = Mol2Reader(mol2_water, from_file=True)

        self.cml_phenol = CmlReader(cml_phenol, ns=None)
        self.mol2_phenol = Mol2Reader(mol2_phenol, from_file=True)

        self.cml_methyl_pyrazole = CmlReader(cml_methyl_pyrazole, ns=None)
        self.mol2_methyl_pyrazole = Mol2Reader(
            mol2_methyl_pyrazole, from_file=True)

        self.cml_thiazolone = CmlReader(cml_thiazolone, ns=None)
        self.mol2_thiazolone = Mol2Reader(mol2_thiazolone, from_file=True)

        self.cml_n_nitromethanamine = CmlReader(
            cml_n_nitromethanamine, ns=None)
        self.mol2_n_nitromethanamine = Mol2Reader(
            mol2_n_nitromethanamine, from_file=True)

        self.cml_urea = CmlReader(cml_urea, ns=None)
        self.mol2_urea = Mol2Reader(mol2_urea, from_file=True)

        self.cml_methylthiourea = CmlReader(cml_methylthiourea, ns=None)
        self.mol2_methylthiourea = Mol2Reader(
            mol2_methylthiourea, from_file=True)

        self.cml_aniline_like = CmlReader(cml_aniline_like, ns=None)
        self.mol2_aniline_like = Mol2Reader(mol2_aniline_like, from_file=True)

        self.cml_dmso = CmlReader(cml_dmso, ns=None)
        self.mol2_dmso = Mol2Reader(mol2_dmso, from_file=True)

        self.cml_phosphonate = CmlReader(cml_phosphonate, ns=None)
        self.mol2_phosphonate = Mol2Reader(mol2_phosphonate, from_file=True)

        self.cml_sulfonate = CmlReader(cml_sulfonate, ns=None)
        self.mol2_sulfonate = Mol2Reader(mol2_sulfonate, from_file=True)

        self.cml_acetic_acid = CmlReader(cml_acetic_acid, ns=None)
        self.mol2_acetic_acid = Mol2Reader(mol2_acetic_acid, from_file=True)

        self.cml_tenoxicam = CmlReader(cml_tenoxicam, ns=None)
        self.mol2_tenoxicam = Mol2Reader(mol2_tenoxicam, from_file=True)

        self.cml_thiadiazole = CmlReader(cml_thiadiazole, ns=None)
        self.mol2_thiadiazole = Mol2Reader(mol2_thiadiazole, from_file=True)

        self.cml_haloethane = CmlReader(cml_haloethane, ns=None)
        self.mol2_haloethane = Mol2Reader(mol2_haloethane, from_file=True)

    def tearDown(self):
        del self.cml_water
        del self.mol2_water

        del self.cml_phenol
        del self.mol2_phenol

        del self.cml_methyl_pyrazole
        del self.mol2_methyl_pyrazole

        del self.cml_thiazolone
        del self.mol2_thiazolone

        del self.cml_urea
        del self.mol2_urea

        del self.cml_methylthiourea
        del self.mol2_methylthiourea

        del self.cml_n_nitromethanamine
        del self.mol2_n_nitromethanamine

        del self.cml_aniline_like
        del self.mol2_aniline_like

        del self.cml_dmso
        del self.mol2_dmso

        del self.cml_phosphonate
        del self.mol2_phosphonate

        del self.cml_sulfonate
        del self.mol2_sulfonate

        del self.cml_acetic_acid
        del self.mol2_acetic_acid

        del self.cml_tenoxicam
        del self.mol2_tenoxicam

        del self.cml_thiadiazole
        del self.mol2_thiadiazole

        del self.cml_haloethane
        del self.mol2_haloethane

    def test_water_assignment(self):
        LOGGER.info("TESTING WATER ASSIGNMENT")
        network = self.cml_water.get_cml_network_with_sybyl(self.mol2_water)
        match = assign_aip_atom_types(network)
        self.assertEqual(match["a1"], "O.3.water")

    def test_phenol_assignment(self):
        LOGGER.info("TESTING PHENOL ASSIGNMENT")
        network = self.cml_phenol.get_cml_network_with_sybyl(self.mol2_phenol)
        match = assign_aip_atom_types(network)
        self.assertEqual(match["a7"], "O.3.any")

    def test_methyl_pyrazaole_assignment(self):
        LOGGER.info("TESTING METHYL-PYRAZOLE ASSIGNMENT")
        network = self.cml_methyl_pyrazole.get_cml_network_with_sybyl(
            self.mol2_methyl_pyrazole)
        match = assign_aip_atom_types(network)
        self.assertEqual(match["a4"], "N.ar.no_lp")
        self.assertEqual(match["a5"], "N.ar")

    def test_thiazolone_assignment(self):
        LOGGER.info("TESTING THIAZOLONE ASSIGNMENT")
        network = self.cml_thiazolone.get_cml_network_with_sybyl(
            self.mol2_thiazolone)
        match = assign_aip_atom_types(network)
        self.assertEqual(match["a4"], "N.ar.no_lp")
        self.assertEqual(match["a6"], "O.2.am")

    def test_n_nitromethanamine_assignment(self):
        LOGGER.info("TESTING N-NITROMETHANAMINE ASSIGNMENT")
        network = self.cml_n_nitromethanamine.get_cml_network_with_sybyl(
            self.mol2_n_nitromethanamine)
        match = assign_aip_atom_types(network)
        self.assertEqual(match["a2"], "N.pl3.secondary")
        self.assertEqual(match["a3"], "N.pl3.nitro")
        self.assertEqual(match["a4"], "O.2.nitro")
        self.assertEqual(match["a5"], "O.2.nitro")

    def test_urea_assignment(self):
        LOGGER.info("TESTING UREA ASSIGNMENT")
        network = self.cml_urea.get_cml_network_with_sybyl(self.mol2_urea)
        match = assign_aip_atom_types(network)
        self.assertEqual(match["a1"], "N.pl3.am")
        self.assertEqual(match["a2"], "C.2")
        self.assertEqual(match["a3"], "O.2.am")
        self.assertEqual(match["a4"], "N.pl3.am")
        self.assertEqual(match["a5"], "H.N")
        self.assertEqual(match["a6"], "H.N")
        self.assertEqual(match["a7"], "H.N")
        self.assertEqual(match["a8"], "H.N")

    def test_methylthiourea_assignment(self):
        LOGGER.info("TESTING METHYL THIOUREA ASSIGNMENT")
        network = self.cml_methylthiourea.get_cml_network_with_sybyl(
            self.mol2_methylthiourea)
        match = assign_aip_atom_types(network)
        self.assertEqual(match["a1"], "N.pl3.am")
        self.assertEqual(match["a2"], "C.2")
        self.assertEqual(match["a3"], "S.2")
        self.assertEqual(match["a4"], "N.pl3.am")
        self.assertEqual(match["a5"], "C.3")
        self.assertEqual(match["a6"], "H.N")
        self.assertEqual(match["a7"], "H.N")
        self.assertEqual(match["a8"], "H.N")
        self.assertEqual(match["a9"], "H.soft")
        self.assertEqual(match["a10"], "H.soft")
        self.assertEqual(match["a11"], "H.soft")

    def test_aniline_like_assignment(self):
        LOGGER.info("TESTING ANILINE-LIKE ASSIGNMENT")
        network = self.cml_aniline_like.get_cml_network_with_sybyl(
            self.mol2_aniline_like)
        match = assign_aip_atom_types(network)
        self.assertEqual(match["a1"], "N.pl3.aniline")
        self.assertEqual(match["a2"], "C.ar")
        self.assertEqual(match["a3"], "N.ar")
        self.assertEqual(match["a4"], "C.ar")
        self.assertEqual(match["a5"], "C.ar")
        self.assertEqual(match["a6"], "C.ar")
        self.assertEqual(match["a7"], "N.ar")

    def test_dmso_assignment(self):
        LOGGER.info("TESTING DMSO ASSIGNMENT")
        network = self.cml_dmso.get_cml_network_with_sybyl(self.mol2_dmso)
        match = assign_aip_atom_types(network)
        self.assertEqual(match["a1"], "C.3")
        self.assertEqual(match["a2"], "S.O")
        self.assertEqual(match["a3"], "C.3")
        self.assertEqual(match["a4"], "O.2.sulfoxide")

    def test_phosphonate_assignment(self):
        LOGGER.info("TESTING PHOSPHONATE ASSIGNMENT")
        network = self.cml_phosphonate.get_cml_network_with_sybyl(
            self.mol2_phosphonate)
        match = assign_aip_atom_types(network)
        self.assertEqual(match["a1"], "C.3")
        self.assertEqual(match["a2"], "O.2.one_lp")
        self.assertEqual(match["a3"], "P.3")
        self.assertEqual(match["a4"], "C.3")
        self.assertEqual(match["a5"], "O.2.po")
        self.assertEqual(match["a6"], "O.2.one_lp")
        self.assertEqual(match["a7"], "C.3")

    def test_sulfonate_assignment(self):
        LOGGER.info("TESTING SULFONE ASSIGNMENT")
        network = self.cml_sulfonate.get_cml_network_with_sybyl(
            self.mol2_sulfonate)
        match = assign_aip_atom_types(network)
        self.assertEqual(match["a1"], "C.3")
        self.assertEqual(match["a2"], "O.2.one_lp")
        self.assertEqual(match["a3"], "S.O2")
        self.assertEqual(match["a4"], "C.3")
        self.assertEqual(match["a5"], "O.2.sulfone")
        self.assertEqual(match["a6"], "O.2.sulfone")

    def test_acetic_acid_assignment(self):
        LOGGER.info("TESTING ACETIC ACID ASSIGNMENT")
        network = self.cml_acetic_acid.get_cml_network_with_sybyl(
            self.mol2_acetic_acid)
        match = assign_aip_atom_types(network)
        self.assertEqual(match["a1"], "C.3")
        self.assertEqual(match["a2"], "C.2")
        self.assertEqual(match["a3"], "O.2.carbonyl")
        self.assertEqual(match["a4"], "O.2.one_lp")

    def test_tenoxicam_assignment(self):
        LOGGER.info("TESTING TENOXICAM ASSIGNMENT")
        network = self.cml_tenoxicam.get_cml_network_with_sybyl(
            self.mol2_tenoxicam)
        match = assign_aip_atom_types(network)
        self.assertEqual(match["a1"], "C.3")
        self.assertEqual(match["a2"], "N.pl3.tertiary")
        self.assertEqual(match["a3"], "C.2")
        self.assertEqual(match["a4"], "C.2")
        self.assertEqual(match["a5"], "O.2.one_lp")
        self.assertEqual(match["a6"], "C.ar")
        self.assertEqual(match["a7"], "S.2.phene")
        self.assertEqual(match["a8"], "C.ar")
        self.assertEqual(match["a9"], "C.ar")
        self.assertEqual(match["a10"], "C.ar")
        self.assertEqual(match["a11"], "S.O2")
        self.assertEqual(match["a12"], "O.2.sulfone")
        self.assertEqual(match["a13"], "O.2.sulfone")
        self.assertEqual(match["a14"], "C.2")
        self.assertEqual(match["a15"], "O.2.am")
        self.assertEqual(match["a16"], "N.pl3.am")
        self.assertEqual(match["a17"], "C.ar")
        self.assertEqual(match["a18"], "C.ar")
        self.assertEqual(match["a19"], "C.ar")
        self.assertEqual(match["a20"], "C.ar")
        self.assertEqual(match["a21"], "C.ar")
        self.assertEqual(match["a22"], "N.ar")
        self.assertEqual(match["a23"], "H.soft")
        self.assertEqual(match["a24"], "H.soft")
        self.assertEqual(match["a25"], "H.soft")
        self.assertEqual(match["a26"], "H.O")
        self.assertEqual(match["a27"], "H.soft")
        self.assertEqual(match["a28"], "H.soft")
        self.assertEqual(match["a29"], "H.N")
        self.assertEqual(match["a30"], "H.soft")
        self.assertEqual(match["a31"], "H.soft")
        self.assertEqual(match["a32"], "H.soft")
        self.assertEqual(match["a33"], "H.soft")

    def test_thiadiazole_assignment(self):
        LOGGER.info("TESTING THIADIAZOLE ASSIGNMENT")
        network = self.cml_thiadiazole.get_cml_network_with_sybyl(
            self.mol2_thiadiazole)
        match = assign_aip_atom_types(network)
        self.assertEqual(match["a1"], "C.3")
        self.assertEqual(match["a2"], "N.3.tertiary")
        self.assertEqual(match["a3"], "N.3.secondary")
        self.assertEqual(match["a4"], "C.3")
        self.assertEqual(match["a5"], "C.2")
        self.assertEqual(match["a6"], "S.2")
        self.assertEqual(match["a7"], "S.3")

    def test_haloethane_assignment(self):
        LOGGER.info("TESTING HALOETHANE ASSIGNMENT")
        network = self.cml_haloethane.get_cml_network_with_sybyl(
            self.mol2_haloethane)
        match = assign_aip_atom_types(network)
        self.assertEqual(match["a1"], "F")
        self.assertEqual(match["a2"], "C.3")
        self.assertEqual(match["a3"], "Cl")
        self.assertEqual(match["a4"], "I")
        self.assertEqual(match["a5"], "C.3")
        self.assertEqual(match["a6"], "Br")
        self.assertEqual(match["a7"], "H.soft")


if __name__ == '__main__':
    unittest.main()

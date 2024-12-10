import pathlib
import unittest
import numpy as np
import logging
import pandas as pd
from aip_footprinting.atom_class import Atom
from aip_footprinting.define_AIP import define_extreme_AIP, define_hydrogen, get_AIP_value
from aip_footprinting.linear_fit_aip import alpha_linear_104, alpha_linear_002, \
    beta_linear_002, beta_linear_03
logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.INFO)


class TestDefineAIP(unittest.TestCase):

    def setUp(self):
        parent_directory = pathlib.Path(__file__).resolve().parents[0]
        path = "test_files/select_surfaces/"
        self.H_01 = pd.read_csv(
            "{}/{}water_H_01.csv".format(parent_directory, path), index_col=0)
        self.N1_03 = pd.read_csv(
            "{}/{}benzonitrile_N_03.csv".format(parent_directory, path), index_col=0)
        self.C1_002 = pd.read_csv(
            "{}/{}benzonitrile_C1_002.csv".format(parent_directory, path), index_col=0)


        self.atom_types_polar_neg = ['N.1', 'N.3.ammonia', 'N.3.aniline', 'N.3.primary',
                                     'N.3.secondary', 'N.3.tertiary',
                                     'N.2', 'N.ar', 'O.2.po', 'O.2.sulfone', 'O.2.carbonyl',
                                     'O.2.am', 'O.2.other', 'O.2.aldehyde', 'O.2.sulfoxide',
                                     'O.2.nitro', 'O.2.noxide', 'O.3.water', 'O.3.alcohol',
                                     'O.3.any', 'S.2', 'S.2.ps', 'F']
        self.atom_types_all_neg = ['C.ar', 'Cl', 'Br', 'I', 'F',
                                   'N.1', 'N.3.ammonia', 'N.3.aniline', 'N.3.primary',
                                   'N.3.secondary', 'N.3.tertiary', 'N.sp2.donor',
                                   'N.2', 'N.ar', 'O.2.po', 'O.2.sulfone', 'O.2.carbonyl',
                                   'O.2.am', 'O.2.other', 'O.2.aldehyde', 'O.2.sulfoxide',
                                   'O.2.nitro', 'O.2.noxide', 'O.3.water', 'O.3.alcohol',
                                   'O.3.any', 'S.2', 'S.2.ps', 'S.3', 'Se']
        self.atom_types_pos = ["H.soft", "H.O", "H.O.water", "H.N", "H.N.group2"]

    def test_finding_lone_pair_as_minimum(self):
        LOGGER.info("Testing MEP minimum is designated as a lone pair")
        minimum_index = self.N1_03[self.N1_03.charge ==
                                   self.N1_03.charge.min()].index[0]
        N1 = Atom()
        N1.atom_type = "N.1"
        N1_AIP = define_extreme_AIP(self.N1_03, N1, polar=True)
        self.assertEqual(N1_AIP.mepsindex, minimum_index)

    def test_finding_nonpolar_minimum(self):
        LOGGER.info("Testing MEP minimum is found")
        minimum_index = self.C1_002[self.C1_002.charge ==
                                    self.C1_002.charge.min()].index[0]
        C1 = Atom()
        C1.atom_type = "C.1"
        C1_AIP = define_extreme_AIP(self.C1_002, C1, minimum=True, polar=False)
        self.assertEqual(C1_AIP.mepsindex, minimum_index)

    def test_finding_nonpolar_maximum(self):
        LOGGER.info("Testing MEP maximum is found")
        maximum_index = self.C1_002[self.C1_002.charge ==
                                    self.C1_002.charge.max()].index[0]
        C1 = Atom()
        C1.atom_type = "C.1"
        C1_AIP = define_extreme_AIP(
            self.C1_002, C1, minimum=False, polar=False)
        self.assertEqual(C1_AIP.mepsindex, maximum_index)

    def test_AIP_H_as_maximum(self):
        LOGGER.info("Testing MEP maximum is as hydrogen")
        maximum_index = self.H_01[self.H_01.charge ==
                                  self.H_01.charge.max()].index[0]
        H = Atom()
        H.atom_type = "H.O"
        H_AIP = define_extreme_AIP(self.H_01, H, minimum=False, polar=True)
        self.assertEqual(H_AIP.mepsindex, maximum_index)

    def test_get_AIP_value_positive(self):
        LOGGER.info(
            "Testing correct positive MEP value conversion to AIP value")
        posmepsvalue = 0.1

        for atype in self.atom_types_pos:
            atom = Atom()
            atom.atom_type = atype
            aipvalue_nonpolar = get_AIP_value(atom, posmepsvalue, polar=False)
            aipvalue_polar = get_AIP_value(atom, posmepsvalue, polar=True)
            fit_polar = alpha_linear_104[atype]
            fit_nonpolar = alpha_linear_104[atype]
            
            ref_calc_polar = fit_polar[0] + posmepsvalue * fit_polar[1]
            self.assertAlmostEqual(aipvalue_polar, round(ref_calc_polar, 2))
            ref_calc_nonpolar = fit_nonpolar[0] + posmepsvalue * fit_nonpolar[1]
            self.assertAlmostEqual(aipvalue_nonpolar, round(ref_calc_nonpolar, 2))

    def test_get_AIP_value_negative(self):
        LOGGER.info(
            "Testing correct negative MEP value conversion to AIP value")
        negmepsvalue = -0.1

        for atype in self.atom_types_polar_neg:
            atom = Atom()
            atom.atom_type = atype
            aipvalue_polar = get_AIP_value(atom, negmepsvalue, polar=True)
            fit_polar = beta_linear_03[atype]
            ref_calc_polar = - (fit_polar[0] + negmepsvalue * fit_polar[1])
            self.assertAlmostEqual(aipvalue_polar, round(ref_calc_polar, 2))

        for atype in self.atom_types_all_neg:
            atom = Atom()
            atom.atom_type = atype
            aipvalue_nonpolar = get_AIP_value(atom, negmepsvalue, polar=False)
            if atype in beta_linear_002.keys():
                fit_nonpolar = beta_linear_002[atype]
            else:
                fit_nonpolar = beta_linear_002["C.ar"]
            ref_calc_nonpolar = - (fit_nonpolar[0] +
                                   negmepsvalue * fit_nonpolar[1])
            self.assertAlmostEqual(aipvalue_nonpolar, round(ref_calc_nonpolar, 2))


if __name__ == '__main__':
    unittest.main()

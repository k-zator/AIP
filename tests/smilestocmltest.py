# -*- coding: utf-8 -*-
#    cmlgenerator creates fully namespaced CML for molecules from input structures.
#    Copyright (C) 2019  Mark D. Driver
#
#    cmlgenerator is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Affero General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Affero General Public License for more details.
#
#    You should have received a copy of the GNU Affero General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""
Script to test smilestocml script.

@author: mark
"""

import logging
import unittest
import re
from rdkit.Chem import AllChem
import cmlgenerator.smilestocml as smicml

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.INFO)


class SmilesToCMLTestCase(unittest.TestCase):
    """Test case for smilestocml methods
    """

    def setUp(self):
        """Set up before tests.
        """
        self.maxDiff = None
        self.input_smiles = "CCO"
        self.expected_molblock = r"""
     RDKit          

  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
M  END
"""
        self.expected_molblock1 = """
        RDKit          3D

  9  8  0  0  0  0  0  0  0  0999 V2000
   -0.8883    0.1670   -0.0273 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4658   -0.5116   -0.0368 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4311    0.3229    0.5867 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8487    1.1175   -0.5695 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6471   -0.4704   -0.4896 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1964    0.3978    0.9977 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.7920   -0.7224   -1.0597 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.4246   -1.4559    0.5138 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.4671    1.1550    0.0848 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  1  4  1  0
  1  5  1  0
  1  6  1  0
  2  7  1  0
  2  8  1  0
  3  9  1  0
M  END
"""
        self.expected_molblock2 = r"""
     RDKit          3D

  9  8  0  0  0  0  0  0  0  0999 V2000
\s+(-)?([0-9][.])([0-9]){4}\s+(-)?([0-9][.])([0-9]){4}\s+(-)?([0-9][.])([0-9]){4} C   0  0  0  0  0  0  0  0  0  0  0  0
\s+(-)?([0-9][.])([0-9]){4}\s+(-)?([0-9][.])([0-9]){4}\s+(-)?([0-9][.])([0-9]){4} C   0  0  0  0  0  0  0  0  0  0  0  0
\s+(-)?([0-9][.])([0-9]){4}\s+(-)?([0-9][.])([0-9]){4}\s+(-)?([0-9][.])([0-9]){4} O   0  0  0  0  0  0  0  0  0  0  0  0
\s+(-)?([0-9][.])([0-9]){4}\s+(-)?([0-9][.])([0-9]){4}\s+(-)?([0-9][.])([0-9]){4} H   0  0  0  0  0  0  0  0  0  0  0  0
\s+(-)?([0-9][.])([0-9]){4}\s+(-)?([0-9][.])([0-9]){4}\s+(-)?([0-9][.])([0-9]){4} H   0  0  0  0  0  0  0  0  0  0  0  0
\s+(-)?([0-9][.])([0-9]){4}\s+(-)?([0-9][.])([0-9]){4}\s+(-)?([0-9][.])([0-9]){4} H   0  0  0  0  0  0  0  0  0  0  0  0
\s+(-)?([0-9][.])([0-9]){4}\s+(-)?([0-9][.])([0-9]){4}\s+(-)?([0-9][.])([0-9]){4} H   0  0  0  0  0  0  0  0  0  0  0  0
\s+(-)?([0-9][.])([0-9]){4}\s+(-)?([0-9][.])([0-9]){4}\s+(-)?([0-9][.])([0-9]){4} H   0  0  0  0  0  0  0  0  0  0  0  0
\s+(-)?([0-9][.])([0-9]){4}\s+(-)?([0-9][.])([0-9]){4}\s+(-)?([0-9][.])([0-9]){4} H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  1  4  1  0
  1  5  1  0
  1  6  1  0
  2  7  1  0
  2  8  1  0
  3  9  1  0
M  END
"""
        self.expected_conversion = """<molecule>
 <atomArray>
  <atom id="a1" elementType="C" hydrogenCount="3" x3="-0.888300" y3="0.167000" z3="-0.027300"/>
  <atom id="a2" elementType="C" hydrogenCount="2" x3="0.465800" y3="-0.511600" z3="-0.036800"/>
  <atom id="a3" elementType="O" hydrogenCount="1" x3="1.431100" y3="0.322900" z3="0.586700"/>
 </atomArray>
 <bondArray>
  <bond atomRefs2="a1 a2" order="1"/>
  <bond atomRefs2="a2 a3" order="1"/>
 </bondArray>
</molecule>
"""
        self.expected_cml = """<molecule>
 <atomArray>
  <atom id="a1" elementType="C" x3="(-)?([0-9][.])([0-9]){6}" y3="(-)?([0-9][.])([0-9]){6}" z3="(-)?([0-9][.])([0-9]){6}"/>
  <atom id="a2" elementType="C" x3="(-)?([0-9][.])([0-9]){6}" y3="(-)?([0-9][.])([0-9]){6}" z3="(-)?([0-9][.])([0-9]){6}"/>
  <atom id="a3" elementType="O" x3="(-)?([0-9][.])([0-9]){6}" y3="(-)?([0-9][.])([0-9]){6}" z3="(-)?([0-9][.])([0-9]){6}"/>
  <atom id="a4" elementType="H" x3="(-)?([0-9][.])([0-9]){6}" y3="(-)?([0-9][.])([0-9]){6}" z3="(-)?([0-9][.])([0-9]){6}"/>
  <atom id="a5" elementType="H" x3="(-)?([0-9][.])([0-9]){6}" y3="(-)?([0-9][.])([0-9]){6}" z3="(-)?([0-9][.])([0-9]){6}"/>
  <atom id="a6" elementType="H" x3="(-)?([0-9][.])([0-9]){6}" y3="(-)?([0-9][.])([0-9]){6}" z3="(-)?([0-9][.])([0-9]){6}"/>
  <atom id="a7" elementType="H" x3="(-)?([0-9][.])([0-9]){6}" y3="(-)?([0-9][.])([0-9]){6}" z3="(-)?([0-9][.])([0-9]){6}"/>
  <atom id="a8" elementType="H" x3="(-)?([0-9][.])([0-9]){6}" y3="(-)?([0-9][.])([0-9]){6}" z3="(-)?([0-9][.])([0-9]){6}"/>
  <atom id="a9" elementType="H" x3="(-)?([0-9][.])([0-9]){6}" y3="(-)?([0-9][.])([0-9]){6}" z3="(-)?([0-9][.])([0-9]){6}"/>
 </atomArray>
 <bondArray>
  <bond atomRefs2="a1 a2" order="1"/>
  <bond atomRefs2="a2 a3" order="1"/>
  <bond atomRefs2="a1 a4" order="1"/>
  <bond atomRefs2="a1 a5" order="1"/>
  <bond atomRefs2="a1 a6" order="1"/>
  <bond atomRefs2="a2 a7" order="1"/>
  <bond atomRefs2="a2 a8" order="1"/>
  <bond atomRefs2="a3 a9" order="1"/>
 </bondArray>
</molecule>
"""

    def tearDown(self):
        """Clean up after tests.
        """
        del self.input_smiles
        del self.expected_molblock
        del self.expected_molblock2
        del self.expected_cml

    def test_generate_cml_and_inchikey_from_smiles(self):
        """Test to see if expected InChIKey and cml blocks are returned.
        """
        actual_inchikey, actual_cml = smicml.generate_cml_and_inchikey_from_smiles(
            self.input_smiles
        )
        self.assertEqual("LFQSCWFLJHTTHZ-UHFFFAOYSA-N", actual_inchikey)
        self.assertTrue(type(re.match(self.expected_cml, actual_cml)) is not None)

    def test_generate_cml_from_smiles(self):
        """Test to see if expected cml is procuded.
        """
        actual_cml = smicml.generate_cml_from_smiles(self.input_smiles)
        self.assertTrue(type(re.match(self.expected_cml, actual_cml)) is not None)

    def test_convert_structure_to_cml(self):
        """Test to see if expected cml is produced.
        """
        mol = AllChem.MolFromMolBlock(self.expected_molblock1)
        actual_cml = smicml.convert_structure_to_cml(mol)
        self.assertMultiLineEqual(self.expected_conversion, actual_cml)

    def test_convert_structure_to_molblock(self):
        """Test to see if expected molblock is produced.
        """
        mol, conf_id = smicml.generate_3d_structure_for_mol(
            smicml.generate_mol_from_smiles(self.input_smiles)
        )
        actual_molblock = smicml.convert_structure_to_molblock(mol, confId=conf_id)
        self.assertTrue(
            type(re.match(self.expected_molblock2, actual_molblock)) is not None
        )

    def test_generate_3d_structure_for_mol(self):
        """Test to see if expected mol is produced
        """
        mol, conf_id = smicml.generate_3d_structure_for_mol(
            smicml.generate_mol_from_smiles(self.input_smiles)
        )
        LOGGER.info("mol type: %s", str(type(mol)))
        LOGGER.info("conf_id type: %s", str(type(conf_id)))
        LOGGER.info("conf_id : %s", str(conf_id))
        actual_molblock = AllChem.MolToMolBlock(mol, confId=conf_id)
        self.assertTrue(
            type(re.match(self.expected_molblock2, actual_molblock)) is not None
        )

    def test_generate_mol_from_smiles(self):
        """Test to see if expected molecule is generated.
        """
        mol = smicml.generate_mol_from_smiles(self.input_smiles)
        actual_molblock = AllChem.MolToMolBlock(mol)
        self.assertTrue(
            type(re.match(self.expected_molblock2, actual_molblock)) is not None
        )

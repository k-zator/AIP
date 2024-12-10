import unittest
import pathlib
import logging
import aip_footprinting.AIP_parser as fscript

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.INFO)


class Test_footprinting_script(unittest.TestCase):
    """Test footprinting script's parsing ability"""

    def test_parser(self):
        parent_directory = pathlib.Path(__file__).resolve().parents[0]
        path = "test_files/XLYOFNOQVPJJNP-UHFFFAOYSA-N/XLYOFNOQVPJJNP-UHFFFAOYSA-N"
        cube_polar = "{}/{}_0.0300_merged.cube".format(parent_directory, path)
        cube_middle = "{}/{}_0.0104_merged.cube".format(parent_directory, path)
        cube_nonpolar = "{}/{}_0.0020_merged.cube".format(
            parent_directory, path)
        cml = "{}/{}.cml".format(parent_directory, path)

        p = fscript.create_parser()
        args = p.parse_args(["--cml_file", cml,
                             "--cube_polar", cube_polar,
                             "--cube_nonpolar", cube_nonpolar,
                             "--cube_middle", cube_middle,
                             "--centre_surface_percentile", "70",
                             ])
        LOGGER.info("Testing parser of cube and cml files")
        self.assertEqual(args.cube_polar, cube_polar)
        self.assertEqual(args.cube_nonpolar, cube_nonpolar)
        self.assertEqual(args.cube_middle, cube_middle)
        self.assertEqual(args.cml, cml)
        self.assertEqual(args.centre_surface_percentile, 70)
        LOGGER.info("Testing parser parses all arguments")
        q = fscript.footprint(args)
        self.assertIsNone(q)


if __name__ == '__main__':
    unittest.main()

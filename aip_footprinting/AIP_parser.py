"""linear_fit_aip.use_dnn
Parses input for foot-printing from command line.
@author: Katarzyna Joanna Zator (kz265)
"""

import logging
import argparse
import aip_footprinting.linear_fit_aip as linear_fit_aip


logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.INFO)


def create_parser():
    help_text = 'Footprinting algorithm to find atom interactions points (AIPs) for a small moleucle. \
                 Requires input of three .cube files (at 0.0020, 0.0104, and 0.0300 e.Bohr-3) and .cml file \
                 from which it extracts atom-owned MEP surfaces which it searches for extreme points to be \
                 designated as an AIP. CML file contains atom types for specific calculation.'
    sign_off = 'Author: Katarzyna Joanna Zator <kz265>'

    parser = argparse.ArgumentParser(description=help_text, epilog=sign_off)

    parser.add_argument(
        '--cml_file',
        '-c',
        dest='cml',
        type=str,
        default=".",
        action='store',
        help='cml file containing the molecule along with aipAtomType information',
        metavar='cml'
    )
    parser.add_argument(
        '--cube_polar',
        '-m',
        dest='cube_polar',
        type=str,
        default=".",
        action='store',
        help='cube file containing the molecule''s MEPS information at the 0.0300 e.Bohr-3',
        metavar='cube_polar'
    )
    parser.add_argument(
        '--cube_nonpolar',
        '-n',
        dest='cube_nonpolar',
        type=str,
        default=".",
        action='store',
        help='cube file containing the molecule''s MEPS information at the 0.0020 e.Bohr-3',
        metavar='cube_nonpolar'
    )
    parser.add_argument(
        '--cube_middle',
        '-b',
        dest='cube_middle',
        type=str,
        default=".",
        action='store',
        help='cube file containing the molecule''s MEPS information at the 0.0104 e.Bohr-3',
        metavar='cube_middle'
    )
    parser.add_argument(
        '--write',
        '-w',
        dest='write',
        type=str,
        default=False,
        action='store',
        help='if present, writes out xml file with footprinting information'
    )
    parser.add_argument(
        '--dualAIP',
        '-dA',
        dest='dualAIP',

        default=False,
        action='store_true',
        help='if present, dualAIPs are determined for the pi-systems'
    )
    parser.add_argument(
        '--centre_surface_percentile',
        '-csp',
        dest='centre_surface_percentile',
        type=float,
        default=90,
        action='store',
        help='percentile of surface MEP density below which points are excluded for the identification as AIPs'
    )
    parser.add_argument(
        '--lp_excl_r',
        '-lpr',
        dest='lp_excl_r',
        type=float,
        default=1.5,
        action='store',
        help='radius of lone pair on the outer 0.0020 MEPS, exluding this area from pi system consideration'
    )

    parser.add_argument(
        '--dnn',
        '-d',
        dest='dnn',
        default=False,
        action='store_true',
        help='radius of lone pair on the outer 0.0020 MEPS, exluding this area from pi system consideration'
    )
    parser.set_defaults(func=footprint)
    return parser


def footprint(args):
    """Actually runs foot-printing."""
    if args.dnn:
        linear_fit_aip.use_dnn = True
        from aip_footprinting.DNN_footprinting_script import DNN_footprint_script
        k = DNN_footprint_script(args.cml, args.cube_nonpolar, args.centre_surface_percentile, args.lp_excl_r, args.dualAIP)
        k.write_xml(args.write)

    else:
        linear_fit_aip.use_dnn = False
        from aip_footprinting.AIP_footprinting_script import KMC_footprint_script
        k = KMC_footprint_script(args.cml, args.cube_polar, args.cube_middle,
                             args.cube_nonpolar, args.centre_surface_percentile, args.lp_excl_r, args.dualAIP)
        if args.write !=False:
            k.write_xml(args.write)
def main():
    """Main function. Returns logger outputs if encounters failures."""

    parser = create_parser()
    args = parser.parse_args()
    LOGGER.info("parsed args:")
    LOGGER.info(args)
    processed_args = args.func(args)
    LOGGER.info("processed all args except:")
    LOGGER.info(processed_args)


if __name__ == "__main__":
    main()

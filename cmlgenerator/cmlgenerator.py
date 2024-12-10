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
Script for parsing command line input and running conversion.

:Authors:
    Mark Driver <mdd31>

Attributes
----------

GENERATE_DESCRIPTION : str
    Generate subparser CLI description.

GENERATE_EPILOG : str
    Generate subparser CLI epilog containing examples.
    
CONVERT_DESCRIPTION :str
    Convert subparser CLI description.
"""

import logging
import argparse
import glob
import cmlgenerator.cmlmaker as cmlmaker
import cmlgenerator.cmlnamespacing as cmlname

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.INFO)

GENERATE_DESCRIPTION = """Script to carry out transformation to generate CML with full namespacing.

This can operate in two modes:
    1) single use mode:
        takes either a single SMILES or 3D structure file at the command line
        and outputs 1 cml file.
    2) batch mode:
        takes list of 3D structure files or a file containing a list of SMILES
        (one perline) and performs conversion
"""
GENERATE_EPILOG = """Example usages:

        # For a single SMILES, example of benzene
        python -m cmlgenerator generate -s 'c1ccccc1'
        # For multiple SMILES strings, each on a new line in a file called smiles.txt
        python -m cmlgenerator generate -s -b smiles.txt
"""


CONVERT_DESCRIPTION = """Converison of CML with unqualified attributes to qualified attributes."""

def process_args(args):
    """Process CLI Arguments.

    Parameters
    ----------
    args : argparse.Namespace
        Arguments from reading command line.

    Returns
    -------
    output_values : 
        Result of argument processing- collection of failed conversions.

    """
    return args.func(args)


def process_generator_args(args):
    """Process Generator CLI subparser options.
    
    Faciltates conversion of molecule to CML from specified file types.

    Parameters
    ----------
    args : argparse.ArgumentParser
        cmlgenerator CLI argumnet parser.

    Returns
    -------
    output_values : 
        Result of argument processing- collection of failed conversions.

    """
    output_values = None
    if args.smiles and args.batch:
        LOGGER.info("number of confs per molecule: %i", args.confs)
        smiles_list = parse_smiles_from_file(args.molecule[0])
        for smiles in smiles_list:
            output_values = []
            try:
                cmlmaker.generate_cml_file_from_smiles(
                    smiles_string=smiles, max_confs=args.confs,
                    directory=args.out_dir,
                    add_sybyl = args.sybyl,
                    add_aip_atom_types = args.aip_atom_types,
                    aromatic = args.aromatic
                )
            except Exception as e:
                output_values.append((smiles, e))
    elif args.smiles:
        LOGGER.info("number of confs per molecule: %i", args.confs)
        try:
            cmlmaker.generate_cml_file_from_smiles(
                smiles_string=args.molecule[0], max_confs=args.confs,
                directory=args.out_dir,
                add_sybyl=args.sybyl,
                add_aip_atom_types = args.aip_atom_types,
                aromatic = args.aromatic
            )
        except Exception as e:
            output_values = (args.molecule[0], e)
    else:
        LOGGER.info("molecules: %s", args.molecule)
        output_values = []
        for molecule in args.molecule:
            LOGGER.info("molecule: %s", molecule)
            try:
                cmlmaker.generate_cml_file_from_file(molecule, args.f,
                                                     directory=args.out_dir, add_sybyl=args.sybyl, add_aip_atom_types = args.aip_atom_types, aromatic = args.aromatic)
            except Exception as e:
                output_values.append((molecule, e))
    return output_values

def process_convert_args(args):
    """Process CLI arguments of conversion subparser.

    Parameters
    ----------
    args : argparse.ArgumentParser
        cmlgenerator CLI argumnet parser.

    Returns
    -------
    output_values : 
        Result of argument processing- collection of failed conversions.

    """
    output_values = []
    cml_regex = args.cmlregex
    unnamespaced_files = glob.glob(cml_regex)
    
    for unnamespaced_file in unnamespaced_files:
        try:
            cmlname.convert_unnamespaced_file(unnamespaced_file, args.out_dir)
        except Exception as e:
            output_values.append((unnamespaced_file, e))
    
    
    return output_values

def create_argparser():
    """Create Argument parser for CMLGenerator CLI.

    Returns
    -------
    parser : argparse.ArgumentParser
        parser for cmlgenerator CLI.

    """
    description = """CLI for cmlgenerator module. See submodule help for more information."""
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawDescriptionHelpFormatter)
    subparsers = parser.add_subparsers(title="Subparsers")
    subparsers.add_parser("generate", parents=[create_generate_argparser()],
                          description=GENERATE_DESCRIPTION, epilog=GENERATE_EPILOG,
                          formatter_class=argparse.RawDescriptionHelpFormatter,
                          add_help=False,
                          help="Generate attribute qualified CML from many file formats")
    subparsers.add_parser("convert", parents=[create_convert_parser()],
                          description=CONVERT_DESCRIPTION,
                          formatter_class=argparse.RawDescriptionHelpFormatter,
                          add_help=False,
                          help="Convert attribute unqualified CML to attribute qualified CML.")
    return parser

def create_generate_argparser():
    """Create Generate options parser for the command line interface.

    Returns
    -------
    parser : argparse.ArgumentParser
        parser for cmlgenerator generate CLI subparser.

    """
    parser = argparse.ArgumentParser(description=GENERATE_DESCRIPTION, epilog=GENERATE_EPILOG, formatter_class=argparse.RawDescriptionHelpFormatter)
    smiles_or_file_mode = parser.add_mutually_exclusive_group()
    smiles_or_file_mode.add_argument(
        "-f", type=str, help="Specify files to be read, of the given type."
    )
    smiles_or_file_mode.add_argument(
        "-s", "--smiles", action="store_true", help="Specify input is SMILES"
    )
    parser.add_argument(
        "-b",
        "--batch",
        action="store_true",
        help="Specify if multiple SMILES in file rather than single on command line",
    )
    parser.add_argument("--confs",type=int, help="number of conformers to generate from SMILES", default=1000)
    parser.add_argument(
        "molecule", type=str, nargs="+", help="This is the SMILES or filename"
    )
    parser.add_argument("-o", "--out_dir", type=str, default=None,
                        help="Output directory location")
    
    parser.add_argument(
        "-y",
        "--sybyl",
        action="store_true",
        default=False,
        help="Specify if you want to add sybyl atom type to the cml")

    parser.add_argument(
        "-t",
        "--aip_atom_types",
        action="store_true",
        default=False,
        help="Specify if you want to add aip_atom_types and sybyl to the cml")
    
    parser.add_argument(
        "-a",
        "--aromatic",
        action="store_true",
        default=False,
        help="Specify if you want to add aip_atom_types and sybyl to the cml")

    parser.set_defaults(func=process_generator_args)
    return parser

def create_convert_parser():
    """Create convert options parser for the command line interface.

    Returns
    -------
    parser : argparse.ArgumentParser
        parser for cmlgenerator convert CLI subparser.

    """
    parser = argparse.ArgumentParser(description=CONVERT_DESCRIPTION,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-o", "--out_dir", type=str, default=None,
                        help="Output directory location")
    parser.add_argument("cmlregex", type=str, 
                        help="REGEX expression for CML file location.")
    parser.set_defaults(func=process_convert_args)
    return parser

def parse_smiles_from_file(filename):
    """Parse text file containing a single SMILES per line, returning SMILES.
    
    This strips all whitespace surrounding SMILES so they can be resolved.

    Parameters
    ----------
    filename : str
        SMILES filename.

    Returns
    -------
    list
        list of SMILES strings.

    """
    with open(filename, "r") as smiles_file:
        lines = smiles_file.readlines()
        return [line.strip() for line in lines]


def main():
    """Main function to run CLI.

    Returns
    -------
    None.
        Logger outputs information about failures.

    """
    # Create parser
    parser = create_argparser()
    # parse args
    args = parser.parse_args()
    LOGGER.info("parsed args:")
    LOGGER.info(args)
    processed_args = process_args(args)
    LOGGER.info(processed_args)


if __name__ == "__main__":
    main()

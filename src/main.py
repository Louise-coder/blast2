"""
This script executes a custom BLAST 2 program using the `GappedBlast` class defined in the `gappedblast` module.

Command-line Arguments:
- -d, --db : (required) Path to the FASTA file containing the database.
- -q, --query : (required) Path to the FASTA file containing the query sequence.
- -o, --output : (optional) Path to the output file. The default value is "out".
- -k, --k : (optional) Word length. The default value is 3.
- -m, --matrix : (optional) Substitution matrix to use for scoring alignments. The default value is "blosum62".
- -e, --evalue : (optional) E-value threshold. The default value is 0.001.

Example Usage:
    python main.py -d path/to/db.fasta -q path/to/query.fasta -o output_file -k 4 -m pam250 -e 0.01

This will create an instance of `GappedBlast` with the provided parameters.
"""

import argparse
from argparse import Namespace
import logging

from gappedblast import GappedBlast

# Set up logging configuration
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

# Metadata
__author__ = "Louise LAM"
__email__ = "louise.lam@etu.u-paris.fr"
__year__ = "2024"
__license__ = "Université Paris Cité"


def get_parameters() -> Namespace:
    """Parse command-line arguments and return parameters.

    Returns
    -------
    Namespace
        Object containing the command-line parameters.
    """
    parser = argparse.ArgumentParser(
        description="Your custom BLAST 2 program."
    )
    parser.add_argument(
        "-d",
        "--db",
        type=str,
        required=True,
        help="Path to the FASTA file of the database.",
    )
    parser.add_argument(
        "-q",
        "--query",
        type=str,
        required=True,
        help="Path to the FASTA file of the query.",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="out",
        help="Path to the output file.",
    )
    parser.add_argument(
        "-e",
        "--evalue",
        type=float,
        default=0.001,
        help="E-value threshold.",
    )
    parser.add_argument(
        "-k", "--k", type=int, default=3, help="Word length."
    )
    parser.add_argument(
        "-m",
        "--matrix",
        type=str,
        choices=[
            "blosum45",
            "blosum62",
            "blosum80",
            "pam30",
            "pam70",
            "pam250",
        ],
        default="blosum62",
        help="Choose the substitution matrix: blosum45, blosum62, blosum80, pam30, pam70 or pam250 (default: blosum62)",
    )
    return parser.parse_args()


def display_title():
    """Display a welcome message."""
    msg = "WELCOME to your custom BLAST 2 !"
    border_char = "*"
    border = f"\033[35m{border_char * (len(msg) + 4)}*\033[0m"
    print(border)
    print(f"\033[35m* {msg} *\033[0m")
    print(border)


def display_metadata():
    """Display metadata information."""
    print(f"\033[90mAuthor: {__author__}\033[0m")
    print(f"\033[90mEmail: {__email__}\033[0m")
    print(f"\033[90mYear: {__year__}\033[0m")
    print(f"\033[90mLicense: {__license__}\033[0m")


if __name__ == "__main__":
    display_title()
    display_metadata()
    params = get_parameters()
    gblast = GappedBlast(params)
    gblast.run()

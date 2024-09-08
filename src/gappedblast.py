"""This module defines the `GappedBlast` class."""

from Bio.Align import substitution_matrices
import logging

from database import Database
from sequence import Sequence

# Configure logger
logger = logging.getLogger(__name__)


class GappedBlast:
    """A class to represent the custom BLAST 2 process.

    Attributes
    ----------
    db : Database
        An instance of the `Database` class containing sequence records.
    query : Sequence
        An instance of the `Sequence` class representing the query sequence.
    output : str
        Path to the output file.
    evalue : float
        E-value threshold for BLAST.
    k : int
        Length of the word for BLAST.
    matrix : substitution_matrices.MatrixInfo
        Chosen substitution matrix for scoring alignments.
    """

    def __init__(self, params):
        """Initialize the GappedBlast instance.

        Parameters
        ----------
        params : Namespace
            Object containing the parameters of a blast run.
        """
        self.db_fasta = params.db
        self.query_fasta = params.query
        self.output = params.output
        self.evalue = params.evalue
        self.matrix = substitution_matrices.load(params.matrix.upper())
        self.k = params.k
        Sequence.set_word_length(self.k)
        self.db = None
        self.query = None

    def load_data(self):
        """Load data from the provided FASTA files.

        Notes
        -----
        Attempts to load the database and query sequences.
        If loading fails, prompts the user for new file paths.
        """
        logger.info("Gapped-BLAST: Loading data...")
        ok = False
        while not ok:
            try:
                self.db = Database(self.db_fasta)
                self.query = Sequence.from_fasta(self.query_fasta)
                ok = True
            except FileNotFoundError as e:
                logger.error(f"\033[31m{e}\033[0m")
            except ValueError as e:
                logger.error(f"\033[31m{e}\033[0m")
            finally:
                if not ok:
                    self.db_fasta = input(
                        "\033[33mPlease re-enter the path to the database file: \033[0m"
                    )
                    self.query_fasta = input(
                        "\033[33mPlease re-enter the path to the query file: \033[0m"
                    )

    def run(self):
        """Execute the BLAST process."""
        self.load_data()

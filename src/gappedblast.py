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
    """

    def __init__(self, params):
        """Initialize the GappedBlast instance.

        Parameters
        ----------
        params : Namespace
            Object containing the parameters of a blast run.
        """
        self.output = params.output
        self.evalue = params.evalue
        self.matrix = substitution_matrices.load(params.matrix.upper())
        self.k = params.k
        Sequence.set_word_length(self.k)
        self.db = None
        self.query = None
        self._load_data(params.db, params.query)

    def _load_data(self, db_fasta: str, query_fasta: str):
        """Load data from the provided FASTA files.

        Parameters
        ----------
        db_fasta : str
            Path to the FASTA file of the database.
        query_fasta : str
            Path to the FASTA file of the query.

        Notes
        -----
        Attempts to load the database and query sequences.
        If loading fails, prompts the user for new file paths.
        """
        ok = False
        while not ok:
            try:
                self.db = Database(db_fasta)
                self.query = Sequence.from_fasta(query_fasta)
                ok = True
            except FileNotFoundError as e:
                logger.error(f"\033[31m{e}\033[0m")
            except ValueError as e:
                logger.error(f"\033[31m{e}\033[0m")
            finally:
                if not ok:
                    db_fasta = input(
                        "\033[33mPlease re-enter the path to the database file: \033[0m"
                    )
                    query_fasta = input(
                        "\033[33mPlease re-enter the path to the query file: \033[0m"
                    )

    def run(self):
        """Execute the BLAST process."""
        logger.info("Blast 2 running...")
        # TODO: Implement the BLAST process

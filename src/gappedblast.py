"""This module defines the `GappedBlast` class."""

from collections import defaultdict
import logging
from typing import List, Tuple

from config import Config
from database import Database
from sequence import Sequence

# Configure logger
logger = logging.getLogger(__name__)


class GappedBlast:
    """A class to represent the custom BLAST 2 process.

    Attributes
    ----------
    query_fasta : str
        Path to the FASTA file of the query.
    db_fasta : str
        Path to the FASTA file of the database.
    db : Database
        An instance of the `Database` class containing sequence records.
    query : Sequence
        An instance of the `Sequence` class representing the query sequence.
    output : str
        Path to the output file.
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
        if params.k != 3:
            Config.update_param("K", params.k)
        if params.matrix != "blosum62":
            Config.update_param("MATRIX", params.matrix.upper())
        if params.evalue != 0.001:
            Config.update_param("EVALUE", params.evalue)

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

    def compute_alignment_score(self, word_a: str, word_b: str) -> float:
        """Compute the alignment score between two words.

        Parameters
        ----------
        word_a : str
            The first word.
        word_b : str
            The second word.

        Returns
        -------
        float
            The alignment score between the two words.
        """
        score = 0
        for aa_a, aa_b in zip(word_a, word_b):
            score += self.matrix[aa_a, aa_b]
        return score

    def hits_detection(self) -> Dict[str, List[str]]:
        """Detects hits between the query sequence and the database.

        Returns
        -------
        Dict[str, List[str]]
            Query k-mers are keys and matching database k-mers are values.

        Notes
        -----
        A hit is detected if the alignment score is greater than `T`.
        """
        logger.info("Gapped-BLAST: Indexation...")
        self.db.load_index()
        logger.info("Gapped-BLAST: Searching hits...")
        hits = defaultdict(list)
        for q_word in self.query.words:
            for db_word, db_position in self.db.index.items():
                score = self.compute_alignment_score(q_word, db_word)
                if score <= Config.T:
                    continue
                hits[q_word].append(db_word)
        return hits

    def run(self):
        """Execute the BLAST process."""
        self.load_data()
        self.hits_detection()
        # TODO: Implement hit extension and output generation

"""This module defines the `Database` class."""

from Bio import SeqIO
from collections import defaultdict
import logging
from typing import Dict, List

from sequence import Sequence

# Configure logger
logger = logging.getLogger(__name__)


class Database:
    """A class to represent a database of biological sequences.

    Attributes
    ----------
    records : List[Sequence]
        A list of `Sequence` objects loaded from a FASTA file.
    """

    def __init__(self, fasta_file: str):
        """Initialize the Database instance from a FASTA file.

        Parameters
        ----------
        fasta_file : str
            Path to the FASTA file containing the sequences.
        """
        self.records: List[Sequence] = []
        self._load_records(fasta_file)

    def _load_records(self, fasta_file: str):
        """Store the sequences from the FASTA file in the `records` attribute.

        Parameters
        ----------
        fasta_file : str
            Path to the FASTA file containing the sequences.
        """
        for record in SeqIO.parse(fasta_file, "fasta"):
            self.records.append(Sequence.from_seqrecord(record))
        logger.info(
            f"{len(self.records)} sequences have been loaded into the database."
        )

    def __iter__(self):
        """Provide an iterator over the sequences in the database.

        Returns
        -------
        iterator
            An iterator over the list of sequences.
        """
        return iter(self.records)

    def get_index(self) -> Dict[str, Dict[int, List[int]]]:
        """Create an index of k-mers from the sequences in the database.

        Returns
        -------
        Dict[str, Dict[int, List[int]]]
            A dictionnary mapping k-mers to their positions in the sequences.

        Notes
        -----
        This index maps each k-mer to a dictionary where:
        - The key is the sequence index in the database.
        - The value is a list of positions where the k-mer appears in that sequence.
        """
        logger.info("Creating the index...")
        index = defaultdict(lambda: defaultdict(list))
        for i, record in enumerate(self.records):
            for word, positions in record.words.items():
                index[word][i].extend(positions)
        return index

"""This module defines the `Sequence` class."""

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
import logging
from typing import List

from config import Config

# Configure logger
logger = logging.getLogger(__name__)


class Sequence(SeqRecord):
    """A class to represent a biological sequence.

    The `Sequence` class extends the `SeqRecord` class from Biopython.
    It provides additional methods for handling sequences.

    Attributes
    ----------
    words : dict[str, List[int]]
        A dictionary where keys are k-mers and values are lists of positions.
    """

    def __init__(self, seq: str, id=None, name="", description=""):
        """Initialize a `Sequence` object.

        Parameters
        ----------
        seq : str
            The sequence data as a string.
        id : Optional
            Identifier for the sequence.
        name : str, optional
            Name of the sequence.
        description : str, optional
            Description of the sequence.
        """
        super().__init__(
            seq=seq, id=id, name=name, description=description
        )
        self.words = self._get_words()

    @classmethod
    def from_seqrecord(cls, seq_record: SeqRecord):
        """Create a `Sequence` instance from a `SeqRecord` object.

        Parameters
        ----------
        seq_record : SeqRecord
            The `SeqRecord` object to convert.

        Returns
        -------
        Sequence
            A `Sequence` object initialized from the provided `SeqRecord`.
        """
        return cls(
            seq=seq_record.seq,
            id=seq_record.id,
            name=seq_record.name,
            description=seq_record.description,
        )

    @classmethod
    def from_fasta(cls, fasta_file: str):
        """Create a `Sequence` instance from a FASTA file.

        Parameters
        ----------
        fasta_file : str
            Path to the FASTA file.

        Returns
        -------
        Sequence
            A `Sequence` object initialized from the FASTA file.
        """
        record = next(SeqIO.parse(fasta_file, "fasta"))
        sequence = cls.from_seqrecord(record)
        logger.info(f"{sequence.name} has been loaded as the query.")
        return sequence

    def _get_words(self) -> dict[str, List[int]]:
        """Extract all possible k-mers from the sequence.

        Returns
        -------
        dict[str, List[int]]
            A dictionary where keys are k-mers and values are positions.

        Raises
        ------
        ValueError
            If `k` is greater than the length of the sequence.
        """
        k = Config.K
        if k > len(self.seq):
            raise ValueError(
                "k cannot be greater than the length of the sequence."
            )
        all_words = defaultdict(list)
        for i in range(len(self.seq) - k + 1):
            word = str(self.seq[i : i + k])
            word = word.replace("U", "C")
            all_words[word].append(i)
        return all_words

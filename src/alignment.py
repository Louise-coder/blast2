"""This module defines the `Alignment` class."""

from typing import Dict, List, Self, Tuple

from config import Config


class Alignment:
    """A class to represent an alignment between two sequences.

    Attributes
    ----------
    seq_a : str
        The first sequence.
    seq_b : str
        The second sequence.
    start_a : int
        The starting position of the alignment in the first sequence.
    start_b : int
        The starting position of the alignment in the second sequence.
    len : int
        The length of the alignment.
    """

    def __init__(
        self,
        seq_a: str,
        seq_b: str,
        start_a: int,
        start_b: int,
        length: int,
    ):
        """Initialize an Alignment instance."""
        self.seq_a = seq_a
        self.seq_b = seq_b
        self.start_a = start_a
        self.start_b = start_b
        self.len = length
        self.score = -1

    def __str__(self):
        """Redefining the print behavior of an Alignment.

        Returns
        -------
        str
            A string representation of the alignment.
        """
        fragment_a = self.seq_a[self.start_a : (self.start_a + self.len)]
        fragment_b = self.seq_b[self.start_b : (self.start_b + self.len)]
        msg = f"{fragment_a}\n{fragment_b}\n"
        if self.score != -1:
            msg += f"score = {self.score}\n"
        else:
            msg += "Not computed yet.\n"
        return msg

    @staticmethod
    def compute_ungapped_score(
        peptide_a: str,
        peptide_b: str,
    ) -> float:
        """Compute the alignment score between two peptides.

        Parameters
        ----------
        peptide_a : str
            The first peptide.
        peptide_b : str
            The second peptide.

        Returns
        -------
        float
            The alignment score between the two peptides.
        """
        score = 0
        for aa_a, aa_b in zip(peptide_a, peptide_b):
            score += Config.get_matrix_score(aa_a, aa_b)
        return score


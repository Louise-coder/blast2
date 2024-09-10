"""This module defines the `Alignment` class."""

from typing import Dict, List, Self, Tuple

from config import Config
from sequence import Sequence
from utils import normalize_ungapped_score


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

    def copy(self) -> Self:
        """Create a copy of the current alignment.

        Returns
        -------
        Alignment
            A copy of the current alignment.
        """
        copy = Alignment(
            self.seq_a, self.seq_b, self.start_a, self.start_b, self.len
        )
        if self.score != -1:
            copy.score = self.score
        return copy

    @staticmethod
    def _is_two_hits_valid(
        current_hit: Tuple[int, int], recent_hits: Dict
    ) -> bool:
        """Verify if the current and recent hits are valid for extension.

        Parameters
        ----------
        current_hit : Tuple[int, int]
            The current hit to be validated.
        recent_hits : Dict[int, Tuple[int, int]]
            A dictionary mapping each diagonal to recent hits.

        Returns
        -------
        bool
            Whether the two hits are valid for extension.

        Notes
        -----
        The two hits are valid for extension if they are:
        - on the same diagonal
        - non overlapping
        - within `Config.A` bases of each other
        """
        diagonal = current_hit[1] - current_hit[0]
        if diagonal in recent_hits:
            recent_hit = recent_hits[diagonal]
            q_distance = current_hit[0] - recent_hit[0]
            db_distance = current_hit[1] - recent_hit[1]
            return ((Config.K - 1) < db_distance <= Config.A) and (
                q_distance == db_distance
            )
        return False

    def _one_extension(self):
        """Extend the alignment by one residue (forward and backward).

        Notes
        -----
        The extension cannot exceed the boundaries of the sequences.
        """
        a_start, b_start = self.start_a, self.start_b
        a_end, b_end = (
            self.start_a + self.len - 1,
            self.start_b + self.len - 1,
        )
        a_len, b_len = len(self.seq_a), len(self.seq_b)
        b_step = 0 if ((a_start - 1 < 0) or (b_start - 1 < 0)) else 1
        f_step = 0 if ((a_end + 1 < a_len) or (b_end + 1 < b_len)) else 1
        self.start_a = a_start - b_step
        self.start_b = b_start - b_step
        self.len += b_step + f_step
        self.score += Config.get_matrix_score(
            self.seq_a[self.start_a],
            self.seq_b[self.start_b],
        ) + Config.get_matrix_score(
            self.seq_a[a_end + f_step],
            self.seq_b[b_end + f_step],
        )

    def _ungapped_extension(self) -> Self:
        """Extend the alignment without allowing for gaps.

        Returns
        -------
        Alignment
            The extended alignment with the highest score.
        """
        q_start, db_start = self.start_a, self.start_b
        q_word = str(self.seq_a[q_start : (q_start + self.len)])
        db_word = str(self.seq_b[db_start : (db_start + self.len)])
        self.score = Alignment.compute_ungapped_score(q_word, db_word)
        current, top = self, self.copy()
        Su = 0.1 * top.score
        i = 0
        while i < 20 and (top.score - current.score) <= Su:
            current._one_extension()
            if current.score > top.score:
                top = current.copy()
                Su = 0.1 * top.score
            i += 1
        return top

    @classmethod
    def extend_to_hsp(
        cls,
        q_record: Sequence,
        db_record: Sequence,
        hits: List[Tuple[str, str]],
    ) -> List[Self]:
        """Extend hits between the query and a database sequence without gaps.

        Parameters
        ----------
        q_record : Sequence
            The query `Sequence` being compared to the database.
        db_record : Sequence
            The database `Sequence` being compared to the query.
        hits : List[Tuple[str, str]]
            A list of tuples where each tuple contains a query word and a matching database word.

        Notes
        -----
        This method attempts to extend hits without allowing for gaps by
        comparing each `current_hit` to recently encountered hits stored in
        `recent_hits` which contains for each diagonal, the most recent hit.

        - **Comparison of Hits**:
        - For each pair of query and database word positions,
        - The couple (recent_hit, current_hit) is evaluated to determine if they are valid for extension.
        - If valid, an ungapped extension between the two hits is performed to generate an HSP.
        - If the normalized score is greater than the threshold (`SG`), the alignment is added to the results.

        Returns
        -------
        List[Alignment]
            A list of valid HSP to extend with gaps in the next step.
        """
        recent_hits = {}
        all_hsp = []
        for q_word, db_word in hits:
            for q_pos in q_record.words[q_word]:
                for db_pos in db_record.words[db_word]:
                    diagonal = db_pos - q_pos
                    current_hit = (q_pos, db_pos)
                    if cls._is_two_hits_valid(current_hit, recent_hits):
                        recent_hit = recent_hits[diagonal]
                        hsp = cls(
                            q_record.seq,
                            db_record.seq,
                            recent_hit[0],
                            recent_hit[1],
                            current_hit[0] - recent_hit[0] + Config.K,
                        )
                        extended_hsp = hsp._ungapped_extension()
                        if (
                            normalize_ungapped_score(extended_hsp.score)
                            > Config.SG
                        ):
                            all_hsp.append(extended_hsp)
                    recent_hits[diagonal] = (q_pos, db_pos)
        return all_hsp

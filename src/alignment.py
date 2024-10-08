"""This module defines the `Alignment` class."""

import numpy as np
from numpy import ndarray
from typing import Dict, List, Self, Tuple

from config import Config
from sequence import Sequence
from utils import (
    compute_evalue,
    normalize_gapped_score,
    normalize_ungapped_score,
)


def worker_gapped_extension(hsp: Self) -> Self:
    """Extend a hit with gaps.

    Parameters
    ----------
    hsp : Alignment

    Returns
    -------
    Alignment
        The extended alignment with gaps.
    """
    seed = hsp.find_best_seed()
    forward = Alignment(
        hsp.seq_a[seed[0] :],
        hsp.seq_b[seed[1] :],
        0,
        0,
        1,
    ).needleman_wunsch_local_alignment()
    sub_a = hsp.seq_a[: seed[0]]
    sub_b = hsp.seq_b[: seed[1]]
    backward = Alignment(
        sub_a[::-1],
        sub_b[::-1],
        0,
        0,
        1,
    ).needleman_wunsch_local_alignment()
    backward.seq_a = backward.seq_a[::-1]
    backward.seq_b = backward.seq_b[::-1]
    gapped_alignment = forward.merge(backward)
    gapped_alignment.seq_id = hsp.seq_id
    return gapped_alignment


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
    seq_id : int
        The id of the matching sequence in the database.
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
        self.seq_id = -1

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
            q_distance = abs(current_hit[0] - recent_hit[0])
            db_distance = abs(current_hit[1] - recent_hit[1])
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
        f_step = 0 if ((a_end + 1 >= a_len) or (b_end + 1 >= b_len)) else 1
        self.start_a = a_start - b_step
        self.start_b = b_start - b_step
        self.len += b_step + f_step
        if f_step == 1:
            self.score += Config.get_matrix_score(
                self.seq_a[a_end + f_step],
                self.seq_b[b_end + f_step],
            )
        if b_step == 1:
            self.score += Config.get_matrix_score(
                self.seq_a[a_start - b_step],
                self.seq_b[b_start - b_step],
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
        Su = Config.SU_PERC * top.score
        while (top.score - current.score) <= Su:
            current._one_extension()
            if current.score > top.score:
                top = current.copy()
                Su = Config.SU_PERC * top.score
            if (current.start_a == 0 or current.start_b == 0) and (
                current.start_a + current.len == len(self.seq_a)
                or current.start_b + current.len == len(self.seq_b)
            ):
                break
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
                            extended_hsp.seq_id = db_record.id
                            all_hsp.append(extended_hsp)
                    recent_hits[diagonal] = (q_pos, db_pos)
        return all_hsp

    def find_best_seed(self) -> Tuple[int, int]:
        """Find the best seed for gapped extension.

        Returns
        -------
        Tuple[int, int]
            A tuple containing the seed positions for the query and database.

        Notes
        -----
        The seed is defined as the center of the best 11-ungapped alignment.
        If the HSP is shorter than 11, the seed is the center of the HSP.
        """
        size = 11
        q_seed = self.start_a + (self.len // 2)
        db_seed = self.start_b + (self.len // 2)
        if self.len < size:
            return (q_seed, db_seed)
        top_score = Alignment.compute_ungapped_score(
            self.seq_a[self.start_a : self.start_a + size],
            self.seq_b[self.start_b : self.start_b + size],
        )
        for i in range(1, self.len - size + 1):
            a_start = self.start_a + i
            b_start = self.start_b + i
            current_score = Alignment.compute_ungapped_score(
                self.seq_a[a_start : a_start + size],
                self.seq_b[b_start : b_start + size],
            )
            if current_score > top_score:
                top_score = current_score
                q_seed = a_start + (size // 2)
                db_seed = b_start + (size // 2)
        return (q_seed, db_seed)

    def find_best_seed_2(self) -> Tuple[int, int]:
        """Another version of the best seed finder.

        Notes
        -----
        The seed is defined as the center of hsp.
        """
        q_seed = self.start_a + (self.len // 2)
        db_seed = self.start_b + (self.len // 2)
        return (q_seed, db_seed)

    def needleman_wunsch_local_alignment(self) -> Self:
        """Perform a local alignment using the Needleman-Wunsch algorithm.

        Returns
        -------
        Alignment
            The extended local alignment with gaps between the two sequences.

        Notes
        -----
        F is the matrix of scores for the alignment.
        M is the matrix of scores for the alignment ending in a match.
        Ix is the matrix of scores for the alignment ending in a gap in the query sequence.
        Iy is the matrix of scores for the alignment ending in a gap in the database sequence.
        """
        a_seq, b_seq = self.seq_a, self.seq_b
        a_len, b_len = len(a_seq), len(b_seq)
        gap_opening = Config.GAP_OPENING_PENALTY
        gap_extension = Config.GAP_EXTENSION_PENALTY
        F = np.zeros((b_len + 1, a_len + 1))
        M = np.zeros((b_len + 1, a_len + 1))  # top left
        Ix = np.zeros((b_len + 1, a_len + 1))  # left
        Iy = np.zeros((b_len + 1, a_len + 1))  # top
        Ix[0, 0], Iy[0, 0] = gap_opening, gap_opening
        for i in range(1, b_len + 1):
            Ix[i, 0] = gap_opening + i * gap_extension
            Iy[i, 0] = -np.inf
            M[i, 0] = Iy[i, 0]
        for j in range(1, a_len + 1):
            Iy[0, j] = gap_opening + j * gap_extension
            Ix[0, j] = -np.inf
            M[0, j] = Ix[0, j]
        for i in range(1, b_len + 1):
            for j in range(1, a_len + 1):
                matrix_score = Config.get_matrix_score(
                    a_seq[j - 1], b_seq[i - 1]
                )
                M[i, j] = max(
                    M[i - 1, j - 1] + matrix_score,
                    Ix[i - 1, j - 1] + matrix_score,
                    Iy[i - 1, j - 1] + matrix_score,
                )
                Ix[i, j] = max(
                    M[i - 1, j] + gap_opening,
                    Ix[i - 1, j] + gap_extension,
                )
                Iy[i, j] = max(
                    M[i, j - 1] + gap_opening,
                    Iy[i, j - 1] + gap_extension,
                )
                F[i, j] = max(M[i, j], Ix[i, j], Iy[i, j])
        return self._backtrack(F, M, Ix, Iy)

    def _backtrack(
        self, F: ndarray, M: ndarray, Ix: ndarray, Iy: ndarray
    ) -> Self:
        """Backtrack to find the local alignment.

        Parameters
        ----------
        F : ndarray
            The matrix of scores for the alignment.
        M : ndarray
            The matrix of scores for the alignment ending in a match.
        Ix : ndarray
            The matrix of scores for the alignment ending in a gap in the query sequence.
        Iy : ndarray
            The matrix of scores for the alignment ending in a gap in the database sequence.

        Returns
        -------
        Alignment
            The extended local alignment with gaps between the two sequences.
        """
        i, j = np.unravel_index(np.argmax(F), F.shape)
        aligned_a, aligned_b = [], []
        while F[i, j] > 0 and i > 0 and j > 0:
            current_score = F[i, j]
            if current_score == M[i, j]:
                aligned_a.append(self.seq_a[j - 1])
                aligned_b.append(self.seq_b[i - 1])
                i -= 1
                j -= 1
            elif current_score == Ix[i, j]:
                aligned_a.append(self.seq_a[j - 1])
                aligned_b.append("-")
                j -= 1
            elif current_score == Iy[i, j]:
                aligned_a.append("-")
                aligned_b.append(self.seq_b[i - 1])
                i -= 1
        while j > 0:
            aligned_a.append(self.seq_a[j - 1])
            aligned_b.append("-")
            j -= 1
        while i > 0:
            aligned_a.append("-")
            aligned_b.append(self.seq_b[i - 1])
            i -= 1
        aligned_a.reverse()
        aligned_b.reverse()
        local = Alignment(
            "".join(aligned_a), "".join(aligned_b), 0, 0, len(aligned_a)
        )
        local.score = np.max(F)
        return local

    def merge(self, other: Self) -> Self:
        """Merge two alignments into a single one.

        Parameters
        ----------
        other : Alignment
            The other alignment to merge with the current one.

        Returns
        -------
        Alignment
            The merged alignment.
        """
        seq_a = other.seq_a + self.seq_a
        seq_b = other.seq_b + self.seq_b
        length = other.len + self.len - 1
        res = Alignment(seq_a, seq_b, 0, 0, length)
        res.score = other.score + self.score
        return res

    def compute_statistics(self, q_len: int, db_len: int):
        """Compute the statistics of the alignment.

        Parameters
        ----------
        q_len : int
            The length of the query sequence.
        db_len : int
            The total number of residues in the database.

        Notes
        -----
        The statistics include:
        - The normalized score (for gapped alignments).
        - The e-value.
        - The number of matches, mismatches, and gaps.
        """
        nb_matches, nb_mismatches, nb_gaps = 0, 0, 0
        normalized_score = normalize_gapped_score(self.score)
        evalue = compute_evalue(normalized_score, q_len, db_len)
        for i in range(self.len):
            if self.seq_a[i] == self.seq_b[i]:
                nb_matches += 1
            elif self.seq_a[i] == "-" or self.seq_b[i] == "-":
                nb_gaps += 1
            else:
                nb_mismatches += 1
        self.normalized_score = normalized_score
        self.evalue = evalue
        self.n_matches = nb_matches
        self.n_mismatches = nb_mismatches
        self.n_gaps = nb_gaps

    def keep_best_alignments(
        gapped_alignments: List[Self],
    ) -> List[Self]:
        """Keep only the best alignment for each sequence.

        Parameters
        ----------
        gapped_alignments : List[Self]
            A list of gapped alignments

        Returns
        -------
        List[Self]
            A list of the best alignment for each sequence.
        """
        unique_alignments = {}
        for alignment in gapped_alignments:
            if alignment.seq_id in unique_alignments:
                if (
                    alignment.score
                    > unique_alignments[alignment.seq_id].score
                ):
                    unique_alignments[alignment.seq_id] = alignment
            else:
                unique_alignments[alignment.seq_id] = alignment
        return list(unique_alignments.values())

    def display_results(self, q_record: Sequence, db_record: Sequence):
        """Display the results of an alignment.

        Parameters
        ----------
        db_record : Sequence
            The database `Sequence` being compared to the query.
        """
        print(f"\033[33m>{db_record.name} {db_record.description}\033[0m")
        print(f"Length={len(db_record.seq)}\n")

        print(
            f"Score={self.normalized_score:.1f} bits ({self.score}), \033[36mExpect={self.evalue:.1e}\033[0m"
        )
        print(
            f"Identities={self.n_matches}/{self.len} ({int(self.n_matches*100/self.len)}%), Positives={self.n_mismatches}/{self.len} ({int(self.n_mismatches*100/self.len)}%), Gaps={self.n_gaps}/{self.len} ({int(self.n_gaps*100/self.len)}%)\n"
        )
        q_start = str(q_record.seq.strip()).find(
            self.seq_a.replace("-", "").strip()
        )
        db_start = str(db_record.seq.strip()).find(
            self.seq_b.replace("-", "").strip()
        )
        segment_width = 50
        if self.len < 50:
            print(f"Query:  {q_start:5d}\t{self.seq_a:<50}\t{q_start+50}")
            print(
                f"Sbjct:  {db_start:5d}\t{self.seq_b:<50}\t{db_start+50}\n"
            )
        else:
            total_length = len(self.seq_a)
            for start in range(0, total_length, segment_width):
                end = min(start + segment_width, total_length)
                print(
                    f"Query:  {q_start + start + 1:5d} {self.seq_a[start:end]:<50} {q_start + end:5d}"
                )
                print(
                    f"Sbjct:  {db_start + start + 1:5d} {self.seq_b[start:end]:<50} {db_start + end:5d}\n"
                )

    def get_results(self, q_record: Sequence, db_record: Sequence) -> str:
        """Display the results of an alignment.

        Parameters
        ----------
        db_record : Sequence
            The database `Sequence` being compared to the query.

        Returns
        -------
        str
            A string containing the results of the alignment.
        """
        content = f"{db_record.name} {db_record.description}\n"
        content += f"Length={len(db_record.seq)}\n\n"
        content += f"Score={self.normalized_score:.1f} bits ({self.score}), Expect={self.evalue:.1e}\n"
        content += f"Identities={self.n_matches}/{self.len} ({int(self.n_matches*100/self.len)}%), Positives={self.n_mismatches}/{self.len} ({int(self.n_mismatches*100/self.len)}%), Gaps={self.n_gaps}/{self.len} ({int(self.n_gaps*100/self.len)}%)\n\n"
        q_start = str(q_record.seq.strip()).find(
            self.seq_a.replace("-", "").strip()
        )
        db_start = str(db_record.seq.strip()).find(
            self.seq_b.replace("-", "").strip()
        )
        segment_width = 50
        if self.len <= segment_width:
            content += f"Query:  {q_start + 1:5d} {self.seq_a} {q_start + self.len:5d}\n"
            content += f"Sbjct:  {db_start + 1:5d} {self.seq_b} {db_start + self.len:5d}\n\n"
        else:
            total_length = len(self.seq_a)

            for start in range(0, total_length, segment_width):
                end = min(start + segment_width, total_length)
                content += f"Query:  {q_start + start + 1:5d} {self.seq_a[start:end]:<50} {q_start + end:5d}\n"
                content += f"Sbjct:  {db_start + start + 1:5d} {self.seq_b[start:end]:<50} {db_start + end:5d}\n\n"
        return content + "\n"

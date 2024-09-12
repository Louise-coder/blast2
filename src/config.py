"""This module defines the `Config` class."""

from Bio.Align import substitution_matrices


class Config:
    """Configuration class for Gapped-BLAST parameters.

    Attributes
    ----------
    MATRIX : substitution_matrices
        The substitution matrix used for scoring alignments.
    K : int
        Length of the k-mer used for indexing the database and query sequence.
    T : int
        Minimum score threshold for a single k-mer hit to be considered significant.
    A : int
        Maximum allowable distance between two hits to trigger extension.
    SG : float
        Minimum score threshold for a HSP to be considered significant.
    LU : float
        Lambda value used for normalizing scores in ungapped alignments.
    KU : float
        Kappa value used for normalizing scores in ungapped alignments.
    LG : float
        Lambda value used for normalizing scores in gapped alignments.
    KG : float
        Kappa value used for normalizing scores in gapped alignments.
    GAP_OPENING_PENALTY : int
        Penalty for opening a gap in the alignment.
    GAP_EXTENSION_PENALTY : int
        Penalty for extending a gap in the alignment.
    EVALUE : float
        E-value threshold for reporting significant alignments.
    NB_RESULTS : int
        Number of results to display in the output.
    """

    # Default configuration
    MATRIX = substitution_matrices.load("BLOSUM62")
    K = 3
    T = 11
    A = 40
    SG = 22.0
    LU = 0.314
    KU = 0.131
    LG = 0.267
    KG = 0.0410
    GAP_OPENING_PENALTY = -11
    GAP_EXTENSION_PENALTY = -1
    EVALUE = 0.1
    NB_RESULTS = 5

    def __str__(cls):
        """Redefine __str__ to display the current configuration."""
        result = ["GAPPED-BLAST CONFIGURATION"]
        for attr, value in cls.__class__.__dict__.items():
            if not attr.startswith("__"):
                result.append(f"{attr}: {value}")
        return "\n".join(result)

    @classmethod
    def update_param(cls, param_name: str, new_value: any):
        """Update a class attribute."""
        setattr(cls, param_name, new_value)

    @classmethod
    def get_matrix_score(cls, aa_a: str, aa_b: str) -> float:
        """Get the score of aligning two amino acids.

        Parameters
        ----------
        aa_a : str
            The first amino acid to consider.
        aa_b : str
            The second amino acid to consider.

        Returns
        -------
        float
            The score of aligning the two amino acids.

        Notes
        -----
        The Selenocysteine (U) and Cysteine (C) amino acids are considered equivalent.
        """
        return cls.MATRIX[aa_a.replace("U", "C"), aa_b.replace("U", "C")]

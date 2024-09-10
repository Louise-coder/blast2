"""This module defines the `Config` class."""


class Config:
    """Configuration class for Gapped-BLAST parameters.

    Attributes
    ----------
    MATRIX : str
        Name of the substitution matrix used for scoring alignments.
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
    EVALUE : float
        E-value threshold for reporting significant alignments.
    """

    # Default configuration
    MATRIX = "BLOSUM62"
    K = 3
    T = 11
    A = 10
    SG = 22.0
    LU = 0.3176
    KU = 0.134
    EVALUE = 0.001

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

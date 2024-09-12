from numpy import log

from config import Config


def compute_evalue(
    normalized_score: float, q_len: int, db_len: int
) -> float:
    """Compute the e-value from the normalized score.

    Parameters
    ----------
    normalized_score : float
        The normalized score.

    Returns
    -------
    float
        The resulting e-value.
    """
    evalue = (q_len * db_len) / 2**normalized_score
    return evalue


def normalize_gapped_score(score: float) -> float:
    """Normalize the gapped score using the formula.

    Parameters
    ----------
    score : float
        The gapped score to normalize.

    Returns
    -------
    float
        The normalized score.
    """
    return (Config.LG * score - log(Config.KG)) / log(2)


def normalize_ungapped_score(score: float) -> float:
    """Normalize the ungapped score using the formula.

    Parameters
    ----------
    score : float
        The ungapped score to normalize.

    Returns
    -------
    float
        The normalized score.
    """
    return (Config.LU * score - log(Config.KU)) / log(2)

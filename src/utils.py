from numpy import log

from config import Config


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

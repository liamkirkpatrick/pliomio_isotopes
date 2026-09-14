"""Stable-isotope notation helpers used by the legacy SWIM model.

All delta values are expressed in per mil (‰).  The constants and excess
definitions intentionally match ``legacy_matlab/distillation_2020.m`` during
the parity phase.
"""

from typing import TypeAlias

import numpy as np
from numpy.typing import ArrayLike, NDArray

FloatResult: TypeAlias = np.float64 | NDArray[np.float64]

R18O_VSMOW = 0.00200520
RD_VSMOW = 0.00015576
R17O_VSMOW = 0.0003799


def delta_to_ratio(
    delta_permil: ArrayLike, standard_ratio: float
) -> FloatResult:
    """Convert delta notation (‰) to an absolute isotope ratio."""
    delta = np.asarray(delta_permil, dtype=np.float64)
    return (1.0 + delta / 1000.0) * standard_ratio


def ratio_to_delta(
    isotope_ratio: ArrayLike, standard_ratio: float
) -> FloatResult:
    """Convert an absolute isotope ratio to delta notation (‰)."""
    ratio = np.asarray(isotope_ratio, dtype=np.float64)
    return (ratio / standard_ratio - 1.0) * 1000.0


def log_delta(delta_permil: ArrayLike) -> FloatResult:
    """Return the logarithmic delta value used by SWIM, in per mil.

    This is ``1000 * ln(1 + delta / 1000)``.  Keeping the factor of 1000 is
    important because the legacy logarithmic d-excess polynomial assumes it.
    """
    delta = np.asarray(delta_permil, dtype=np.float64)
    return 1000.0 * np.log1p(delta / 1000.0)


def linear_deuterium_excess(
    delta_d_permil: ArrayLike, delta_18o_permil: ArrayLike
) -> FloatResult:
    """Return conventional deuterium excess, ``delta-D - 8*delta-18O`` (‰)."""
    delta_d = np.asarray(delta_d_permil, dtype=np.float64)
    delta_18o = np.asarray(delta_18o_permil, dtype=np.float64)
    return delta_d - 8.0 * delta_18o


def logarithmic_deuterium_excess(
    delta_d_permil: ArrayLike, delta_18o_permil: ArrayLike
) -> FloatResult:
    """Return the logarithmic deuterium excess used by legacy SWIM (‰)."""
    delta_d_ln = log_delta(delta_d_permil)
    delta_18o_ln = log_delta(delta_18o_permil)
    return delta_d_ln - (-0.0285 * delta_18o_ln**2 + 8.47 * delta_18o_ln)


def oxygen_17_excess(
    delta_17o_permil: ArrayLike, delta_18o_permil: ArrayLike
) -> FloatResult:
    """Return 17O excess using the legacy exponent 0.528, in per meg."""
    delta_17o = np.asarray(delta_17o_permil, dtype=np.float64)
    delta_18o = np.asarray(delta_18o_permil, dtype=np.float64)
    return 1.0e6 * (
        np.log1p(delta_17o / 1000.0)
        - 0.528 * np.log1p(delta_18o / 1000.0)
    )


__all__ = [
    "R17O_VSMOW",
    "R18O_VSMOW",
    "RD_VSMOW",
    "delta_to_ratio",
    "linear_deuterium_excess",
    "log_delta",
    "logarithmic_deuterium_excess",
    "oxygen_17_excess",
    "ratio_to_delta",
]

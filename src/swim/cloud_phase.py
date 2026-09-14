"""Temperature-dependent cloud liquid and ice fractions used by SWIM."""

from typing import Literal, TypeAlias

import numpy as np
from numpy.typing import ArrayLike, NDArray

CloudPhaseMethod: TypeAlias = Literal["mid", "iir", "adj"]
FloatResult: TypeAlias = np.float64 | NDArray[np.float64]


def cloud_phase_fractions(
    temperature_c: ArrayLike, method: CloudPhaseMethod = "mid"
) -> tuple[FloatResult, FloatResult]:
    """Return cloud ice and liquid fractions in the legacy MATLAB order.

    The three polynomial/logistic fits correspond to the options in
    ``fraction_il_brm_H10.m``.  The active ``distillation_2020.m`` pathway
    explicitly selects ``"adj"``; ``"mid"`` remains the helper's default to
    preserve its standalone MATLAB behavior.
    """
    temperature = np.asarray(temperature_c, dtype=np.float64)

    if method == "mid":
        coefficients = (5.3608, 0.4025, 0.08387, 0.007182, 2.39e-4, 2.87e-6)
    elif method == "iir":
        coefficients = (5.2918, 0.3694, 0.06635, 0.006367, 2.33e-4, 2.97e-6)
    elif method == "adj":
        coefficients = (5.37, 0.4025, 0.0847, 0.007182, 2.39e-4, 2.87e-6)
    else:
        raise ValueError(f"Unknown cloud-phase method: {method!r}")

    c0, c1, c2, c3, c4, c5 = coefficients
    polynomial = (
        c0
        + c1 * temperature
        + c2 * temperature**2
        + c3 * temperature**3
        + c4 * temperature**4
        + c5 * temperature**5
    )
    with np.errstate(over="ignore"):
        fraction_liquid = 1.0 / (1.0 + np.exp(-polynomial))
    fraction_ice = 1.0 - fraction_liquid
    return fraction_ice, fraction_liquid


__all__ = ["CloudPhaseMethod", "cloud_phase_fractions"]

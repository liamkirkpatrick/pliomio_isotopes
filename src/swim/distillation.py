"""Stepwise Rayleigh distillation for the legacy SWIM pathway."""

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from swim.fractionation import EffectiveFractionationFactors
from swim.isotopes import (
    R17O_VSMOW,
    R18O_VSMOW,
    RD_VSMOW,
    linear_deuterium_excess,
    logarithmic_deuterium_excess,
    oxygen_17_excess,
    ratio_to_delta,
)


@dataclass(frozen=True)
class DistillationResult:
    """Vapor and precipitation isotope trajectories."""

    deuterium_vapor_ratio: NDArray[np.float64]
    deuterium_precipitation_ratio: NDArray[np.float64]
    oxygen_18_vapor_ratio: NDArray[np.float64]
    oxygen_18_precipitation_ratio: NDArray[np.float64]
    oxygen_17_vapor_ratio: NDArray[np.float64]
    oxygen_17_precipitation_ratio: NDArray[np.float64]
    delta_d_precipitation: NDArray[np.float64]
    delta_18o_precipitation: NDArray[np.float64]
    delta_17o_precipitation: NDArray[np.float64]
    delta_d_log: NDArray[np.float64]
    delta_18o_log: NDArray[np.float64]
    delta_17o_log: NDArray[np.float64]
    deuterium_excess: NDArray[np.float64]
    logarithmic_deuterium_excess: NDArray[np.float64]
    oxygen_17_excess_per_meg: NDArray[np.float64]


def _as_trajectory(value: ArrayLike, size: int, name: str) -> NDArray[np.float64]:
    array = np.asarray(value, dtype=np.float64)
    if array.shape != (size,):
        message = f"{name} must be a one-dimensional trajectory of length {size}"
        raise ValueError(message)
    return array


def _integrate_isotope_ratio(
    fraction_vapor_remaining: NDArray[np.float64],
    effective_fractionation: NDArray[np.float64],
    initial_vapor_ratio: float,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    vapor = np.full(fraction_vapor_remaining.shape, np.nan, dtype=np.float64)
    precipitation = np.full(fraction_vapor_remaining.shape, np.nan, dtype=np.float64)
    vapor[0] = initial_vapor_ratio
    precipitation[0] = vapor[0] * effective_fractionation[0]

    for index in range(1, fraction_vapor_remaining.size):
        change_log_ratio = (effective_fractionation[index] - 1.0) * (
            np.log(fraction_vapor_remaining[index])
            - np.log(fraction_vapor_remaining[index - 1])
        )
        vapor[index] = np.exp(np.log(vapor[index - 1]) + change_log_ratio)
        if change_log_ratio != 0.0:
            precipitation[index] = vapor[index] * effective_fractionation[index]

    return vapor, precipitation


def rayleigh_distillation(
    fraction_vapor_remaining: ArrayLike,
    effective_fractionation: EffectiveFractionationFactors,
    *,
    initial_delta_d_permil: float,
    initial_delta_18o_permil: float,
    initial_oxygen_17_excess_log: float,
) -> DistillationResult:
    """Run the active MATLAB stepwise Rayleigh update.

    ``initial_oxygen_17_excess_log`` is dimensionless, matching the value
    returned by ``evaporation_2021.m``.  Precipitation 17O excess in the result
    is reported in per meg, matching ``distillation_2020.m``.
    """
    fraction = np.asarray(fraction_vapor_remaining, dtype=np.float64)
    if fraction.ndim != 1 or fraction.size == 0:
        raise ValueError("fraction_vapor_remaining must be a nonempty 1-D array")
    size = fraction.size
    alpha_d = _as_trajectory(effective_fractionation.deuterium, size, "deuterium")
    alpha_18o = _as_trajectory(
        effective_fractionation.oxygen_18, size, "oxygen_18"
    )
    alpha_17o = _as_trajectory(
        effective_fractionation.oxygen_17, size, "oxygen_17"
    )

    initial_d_ratio = (1.0 + initial_delta_d_permil / 1000.0) * RD_VSMOW
    initial_18o_ratio = (1.0 + initial_delta_18o_permil / 1000.0) * R18O_VSMOW
    initial_delta_17o = 1000.0 * (
        np.exp(
        initial_oxygen_17_excess_log
            + 0.528 * np.log(1.0 + initial_delta_18o_permil / 1000.0)
        )
        - 1.0
    )
    initial_17o_ratio = (1.0 + initial_delta_17o / 1000.0) * R17O_VSMOW

    d_vapor, d_precipitation = _integrate_isotope_ratio(
        fraction, alpha_d, initial_d_ratio
    )
    oxygen_18_vapor, oxygen_18_precipitation = _integrate_isotope_ratio(
        fraction, alpha_18o, initial_18o_ratio
    )
    oxygen_17_vapor, oxygen_17_precipitation = _integrate_isotope_ratio(
        fraction, alpha_17o, initial_17o_ratio
    )

    delta_d = np.asarray(ratio_to_delta(d_precipitation, RD_VSMOW))
    delta_18o = np.asarray(ratio_to_delta(oxygen_18_precipitation, R18O_VSMOW))
    delta_17o = np.asarray(ratio_to_delta(oxygen_17_precipitation, R17O_VSMOW))
    delta_d_log = 1000.0 * np.log(d_precipitation / RD_VSMOW)
    delta_18o_log = 1000.0 * np.log(oxygen_18_precipitation / R18O_VSMOW)
    delta_17o_log = 1000.0 * np.log(oxygen_17_precipitation / R17O_VSMOW)

    return DistillationResult(
        deuterium_vapor_ratio=d_vapor,
        deuterium_precipitation_ratio=d_precipitation,
        oxygen_18_vapor_ratio=oxygen_18_vapor,
        oxygen_18_precipitation_ratio=oxygen_18_precipitation,
        oxygen_17_vapor_ratio=oxygen_17_vapor,
        oxygen_17_precipitation_ratio=oxygen_17_precipitation,
        delta_d_precipitation=delta_d,
        delta_18o_precipitation=delta_18o,
        delta_17o_precipitation=delta_17o,
        delta_d_log=delta_d_log,
        delta_18o_log=delta_18o_log,
        delta_17o_log=delta_17o_log,
        deuterium_excess=np.asarray(linear_deuterium_excess(delta_d, delta_18o)),
        logarithmic_deuterium_excess=np.asarray(
            logarithmic_deuterium_excess(delta_d, delta_18o)
        ),
        oxygen_17_excess_per_meg=np.asarray(
            oxygen_17_excess(delta_17o, delta_18o)
        ),
    )


__all__ = ["DistillationResult", "rayleigh_distillation"]

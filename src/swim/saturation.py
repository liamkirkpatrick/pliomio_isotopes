"""Saturation vapor-pressure formulas used by the legacy SWIM model."""

from typing import TypeAlias

import numpy as np
from numpy.typing import ArrayLike, NDArray

FloatResult: TypeAlias = np.float64 | NDArray[np.float64]

CELSIUS_TO_KELVIN = 273.15
WATER_VAPOR_TO_DRY_AIR_GAS_CONSTANT_RATIO = 0.622


def saturated_mixing_ratio(
    saturation_vapor_pressure_kpa: ArrayLike, pressure_kpa: ArrayLike
) -> FloatResult:
    """Return saturated water-vapor mixing ratio (kg/kg dry air).

    Both pressure inputs must use the same units; SWIM conventionally uses
    kPa.  The expression preserves paper Eq. (A2) and the legacy value
    ``R_d / R_wv = 0.622``.
    """
    vapor_pressure = np.asarray(saturation_vapor_pressure_kpa, dtype=np.float64)
    pressure = np.asarray(pressure_kpa, dtype=np.float64)
    return WATER_VAPOR_TO_DRY_AIR_GAS_CONSTANT_RATIO * vapor_pressure / (
        pressure - vapor_pressure
    )


def saturation_vapor_pressure_liquid(
    temperature_c: ArrayLike,
) -> FloatResult:
    """Return saturation vapor pressure over liquid water, in kPa.

    ``temperature_c`` is in degrees Celsius.  The expression is the Murphy and
    Koop parameterization copied from ``pseudo_adiabat_function.m``.
    """
    temperature_k = np.asarray(temperature_c, dtype=np.float64) + CELSIUS_TO_KELVIN
    log_temperature = np.log(temperature_k)
    return 1.0e-3 * np.exp(
        54.842763
        - 6763.22 / temperature_k
        - 4.21 * log_temperature
        + 0.000367 * temperature_k
        + np.tanh(0.0415 * (temperature_k - 218.8))
        * (
            53.878
            - 1331.22 / temperature_k
            - 9.44523 * log_temperature
            + 0.014025 * temperature_k
        )
    )


def saturation_vapor_pressure_ice(
    temperature_c: ArrayLike,
) -> FloatResult:
    """Return saturation vapor pressure over ice, in kPa.

    ``temperature_c`` is in degrees Celsius.  The expression is the Murphy and
    Koop parameterization copied from ``pseudo_adiabat_function.m``.
    """
    temperature_k = np.asarray(temperature_c, dtype=np.float64) + CELSIUS_TO_KELVIN
    return 1.0e-3 * np.exp(
        9.550426
        - 5723.265 / temperature_k
        + 3.53068 * np.log(temperature_k)
        - 0.00728332 * temperature_k
    )


__all__ = [
    "CELSIUS_TO_KELVIN",
    "WATER_VAPOR_TO_DRY_AIR_GAS_CONSTANT_RATIO",
    "saturated_mixing_ratio",
    "saturation_vapor_pressure_ice",
    "saturation_vapor_pressure_liquid",
]

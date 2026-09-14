"""Pseudo-adiabatic thermodynamics for the legacy SWIM pathway."""

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from swim.saturation import (
    WATER_VAPOR_TO_DRY_AIR_GAS_CONSTANT_RATIO,
    saturated_mixing_ratio,
    saturation_vapor_pressure_ice,
    saturation_vapor_pressure_liquid,
)

LATENT_HEAT_VAPORIZATION_J_KG = 2.501e6
DRY_AIR_GAS_CONSTANT_J_KG_K = 287.053
DRY_AIR_HEAT_CAPACITY_J_KG_K = 1004.0
DRY_AIR_TO_WATER_VAPOR_HEAT_CAPACITY_RATIO = 0.5427
WATER_MOLES_PER_KILOGRAM = 55.5084350618
CELSIUS_TO_KELVIN = 273.15


@dataclass(frozen=True)
class PseudoAdiabatResult:
    """Outputs from the legacy pseudo-adiabatic Euler integration."""

    fraction_vapor_remaining: NDArray[np.float64]
    pressure_kpa: NDArray[np.float64]
    saturation_vapor_pressure_kpa: NDArray[np.float64]
    saturated_mixing_ratio: NDArray[np.float64]
    pressure_gradient_kpa_per_k: NDArray[np.float64]
    saturation_over_ice: NDArray[np.float64]
    saturation_over_liquid: NDArray[np.float64]


def _validated_trajectory_inputs(
    temperature_c: ArrayLike,
    fraction_ice: ArrayLike,
    fraction_liquid: ArrayLike,
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    temperature = np.asarray(temperature_c, dtype=np.float64)
    ice = np.asarray(fraction_ice, dtype=np.float64)
    liquid = np.asarray(fraction_liquid, dtype=np.float64)
    if temperature.ndim != 1:
        raise ValueError("temperature_c must be one-dimensional")
    if temperature.size < 2:
        raise ValueError("temperature_c must contain at least two points")
    if ice.shape != temperature.shape or liquid.shape != temperature.shape:
        raise ValueError("temperature and cloud-phase arrays must have equal shapes")
    return temperature, ice, liquid


def _effective_thermodynamic_properties(
    temperature_k: NDArray[np.float64],
    fraction_ice: NDArray[np.float64],
    fraction_liquid: NDArray[np.float64],
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    latent_heat_ice = (
        46782.5
        + 35.8925 * temperature_k
        - 0.07414 * temperature_k**2
        + 541.5 * np.exp(-(temperature_k / 123.75) ** 2)
    ) * WATER_MOLES_PER_KILOGRAM
    heat_capacity_ice = (
        -2.0572
        + 0.14644 * temperature_k
        + 0.06163 * temperature_k * np.exp(-(temperature_k / 125.1) ** 2)
    )
    heat_capacity_vapor = (
        DRY_AIR_HEAT_CAPACITY_J_KG_K
        / DRY_AIR_TO_WATER_VAPOR_HEAT_CAPACITY_RATIO
    )
    effective_latent_heat = (
        fraction_liquid * LATENT_HEAT_VAPORIZATION_J_KG
        + fraction_ice * latent_heat_ice
    )
    effective_heat_capacity_ratio = DRY_AIR_HEAT_CAPACITY_J_KG_K / (
        fraction_liquid * heat_capacity_vapor
        + fraction_ice * heat_capacity_ice
    )
    return effective_latent_heat, effective_heat_capacity_ratio


def _pressure_gradient(
    pressure_kpa: np.float64,
    temperature_k: np.float64,
    mixing_ratio: np.float64,
    effective_latent_heat: np.float64,
    effective_heat_capacity_ratio: np.float64,
) -> np.float64:
    pressure_coefficient = (
        1.0 + mixing_ratio / WATER_VAPOR_TO_DRY_AIR_GAS_CONSTANT_RATIO
    ) / (1.0 + mixing_ratio / effective_heat_capacity_ratio)
    return (pressure_kpa / pressure_coefficient) * (
        DRY_AIR_HEAT_CAPACITY_J_KG_K
        + effective_latent_heat**2
        * mixing_ratio
        * WATER_VAPOR_TO_DRY_AIR_GAS_CONSTANT_RATIO
        * pressure_coefficient
        / (DRY_AIR_GAS_CONSTANT_J_KG_K * temperature_k**2)
    ) / (
        DRY_AIR_GAS_CONSTANT_J_KG_K * temperature_k
        + effective_latent_heat * mixing_ratio
    )


def mixed_phase_supersaturation(
    temperature_c: ArrayLike,
    initial_pressure_kpa: float,
    fraction_ice: ArrayLike,
    fraction_liquid: ArrayLike,
) -> NDArray[np.float64]:
    """Return mixed-phase saturation relative to ice along an Euler path.

    This is a direct port of ``mixed_phased_supersaturation.m``.  Its pressure
    integration is intentionally separate from :func:`pseudo_adiabat` because
    the MATLAB baseline performs both integrations sequentially.
    """
    temperature, ice, liquid = _validated_trajectory_inputs(
        temperature_c, fraction_ice, fraction_liquid
    )
    temperature_k = temperature + CELSIUS_TO_KELVIN
    effective_latent_heat, effective_heat_capacity_ratio = (
        _effective_thermodynamic_properties(temperature_k, ice, liquid)
    )
    liquid_vapor_pressure = np.asarray(
        saturation_vapor_pressure_liquid(temperature), dtype=np.float64
    )
    ice_vapor_pressure = np.asarray(
        saturation_vapor_pressure_ice(temperature), dtype=np.float64
    )
    mixed_vapor_pressure = (
        liquid * liquid_vapor_pressure + ice * ice_vapor_pressure
    )

    pressure = np.full(temperature.shape, np.nan, dtype=np.float64)
    mixing_ratio = np.full(temperature.shape, np.nan, dtype=np.float64)
    ice_mixing_ratio = np.full(temperature.shape, np.nan, dtype=np.float64)
    pressure[0] = initial_pressure_kpa

    for index in range(temperature.size - 1):
        mixing_ratio[index] = saturated_mixing_ratio(
            mixed_vapor_pressure[index], pressure[index]
        )
        ice_mixing_ratio[index] = saturated_mixing_ratio(
            ice_vapor_pressure[index], pressure[index]
        )
        gradient = _pressure_gradient(
            pressure[index],
            temperature_k[index],
            mixing_ratio[index],
            effective_latent_heat[index],
            effective_heat_capacity_ratio[index],
        )
        pressure[index + 1] = pressure[index] + gradient * (
            temperature_k[index + 1] - temperature_k[index]
        )

    mixing_ratio[-1] = saturated_mixing_ratio(
        mixed_vapor_pressure[-1], pressure[-1]
    )
    ice_mixing_ratio[-1] = saturated_mixing_ratio(
        ice_vapor_pressure[-1], pressure[-1]
    )
    return mixing_ratio / ice_mixing_ratio


def pseudo_adiabat(
    temperature_c: ArrayLike,
    initial_pressure_kpa: float,
    fraction_ice: ArrayLike,
    fraction_liquid: ArrayLike,
    mixed_phase_saturation: ArrayLike,
    supersaturation_a: float,
    supersaturation_b: float,
    supersaturation_c: float,
) -> PseudoAdiabatResult:
    """Integrate the legacy mixed-phase pseudo-adiabatic trajectory."""
    temperature, ice, liquid = _validated_trajectory_inputs(
        temperature_c, fraction_ice, fraction_liquid
    )
    mixed_saturation = np.asarray(mixed_phase_saturation, dtype=np.float64)
    if mixed_saturation.shape != temperature.shape:
        raise ValueError("mixed_phase_saturation must match temperature shape")

    temperature_k = temperature + CELSIUS_TO_KELVIN
    effective_latent_heat, effective_heat_capacity_ratio = (
        _effective_thermodynamic_properties(temperature_k, ice, liquid)
    )
    liquid_vapor_pressure = np.asarray(
        saturation_vapor_pressure_liquid(temperature), dtype=np.float64
    )
    ice_vapor_pressure = np.asarray(
        saturation_vapor_pressure_ice(temperature), dtype=np.float64
    )
    mixed_vapor_pressure = (
        liquid * liquid_vapor_pressure + ice * ice_vapor_pressure
    )
    prescribed_saturation = np.maximum(
        supersaturation_a
        - supersaturation_b * temperature
        - supersaturation_c * temperature**2,
        1.0,
    )

    pressure = np.full(temperature.shape, np.nan, dtype=np.float64)
    mixing_ratio = np.full(temperature.shape, np.nan, dtype=np.float64)
    liquid_mixing_ratio = np.full(temperature.shape, np.nan, dtype=np.float64)
    ice_mixing_ratio = np.full(temperature.shape, np.nan, dtype=np.float64)
    pressure_gradient = np.full(temperature.shape, np.nan, dtype=np.float64)
    pressure[0] = initial_pressure_kpa

    for index in range(temperature.size - 1):
        liquid_mixing_ratio[index] = saturated_mixing_ratio(
            liquid_vapor_pressure[index], pressure[index]
        )
        ice_mixing_ratio[index] = saturated_mixing_ratio(
            ice_vapor_pressure[index], pressure[index]
        )
        mixed_mixing_ratio = (
            liquid[index] * liquid_mixing_ratio[index]
            + ice[index] * ice_mixing_ratio[index]
        )
        target_saturation = (
            mixed_saturation[index]
            if temperature[index] > 0.0
            else prescribed_saturation[index]
        )
        factor = target_saturation / (
            mixed_mixing_ratio / ice_mixing_ratio[index]
        )
        mixing_ratio[index] = factor * mixed_mixing_ratio
        pressure_gradient[index] = _pressure_gradient(
            pressure[index],
            temperature_k[index],
            mixing_ratio[index],
            effective_latent_heat[index],
            effective_heat_capacity_ratio[index],
        )
        pressure[index + 1] = pressure[index] + pressure_gradient[index] * (
            temperature_k[index + 1] - temperature_k[index]
        )

    liquid_mixing_ratio[-1] = saturated_mixing_ratio(
        liquid_vapor_pressure[-1], pressure[-1]
    )
    ice_mixing_ratio[-1] = saturated_mixing_ratio(
        ice_vapor_pressure[-1], pressure[-1]
    )
    mixed_mixing_ratio = (
        liquid[-1] * liquid_mixing_ratio[-1] + ice[-1] * ice_mixing_ratio[-1]
    )
    # MATLAB tests T(i) here after the loop, where i is the penultimate index.
    # Preserve that endpoint quirk until a deliberate post-parity decision.
    target_saturation = (
        mixed_saturation[-1]
        if temperature[-2] > 0.0
        else prescribed_saturation[-1]
    )
    factor = target_saturation / (mixed_mixing_ratio / ice_mixing_ratio[-1])
    mixing_ratio[-1] = factor * mixed_mixing_ratio
    pressure_gradient[-1] = _pressure_gradient(
        pressure[-1],
        temperature_k[-1],
        mixing_ratio[-1],
        effective_latent_heat[-1],
        effective_heat_capacity_ratio[-1],
    )

    return PseudoAdiabatResult(
        fraction_vapor_remaining=mixing_ratio / mixing_ratio[0],
        pressure_kpa=pressure,
        saturation_vapor_pressure_kpa=mixed_vapor_pressure,
        saturated_mixing_ratio=mixing_ratio,
        pressure_gradient_kpa_per_k=pressure_gradient,
        saturation_over_ice=mixing_ratio / ice_mixing_ratio,
        saturation_over_liquid=mixing_ratio / liquid_mixing_ratio,
    )


__all__ = [
    "PseudoAdiabatResult",
    "mixed_phase_supersaturation",
    "pseudo_adiabat",
]

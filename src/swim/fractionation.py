"""Equilibrium isotope fractionation factors used by legacy SWIM."""

from dataclasses import dataclass
from typing import TypeAlias

import numpy as np
from numpy.typing import ArrayLike, NDArray

FloatResult: TypeAlias = np.float64 | NDArray[np.float64]

CELSIUS_TO_KELVIN = 273.15
O17_LIQUID_EXPONENT = 0.529
O17_ICE_EXPONENT_EVAPORATION = 0.529
O17_ICE_EXPONENT_DISTILLATION = 0.531


@dataclass(frozen=True)
class EquilibriumFractionationFactors:
    """Heavy/light equilibrium fractionation factors for liquid and ice."""

    deuterium_liquid: FloatResult
    deuterium_ice: FloatResult
    oxygen_18_liquid: FloatResult
    oxygen_18_ice: FloatResult
    oxygen_17_liquid: FloatResult
    oxygen_17_ice: FloatResult


@dataclass(frozen=True)
class TransportDiffusivityRatios:
    """Heavy/light molecular diffusivity ratios from Hellmann and Harvey."""

    hdo_over_h2o: FloatResult
    h217o_over_h216o: FloatResult
    h218o_over_h216o: FloatResult


@dataclass(frozen=True)
class KineticCondensationFactors:
    """Ice-vapor kinetic fractionation factors for the three isotopes."""

    deuterium_ice: FloatResult
    oxygen_18_ice: FloatResult
    oxygen_17_ice: FloatResult


@dataclass(frozen=True)
class EffectiveFractionationFactors:
    """Cloud-phase-weighted total fractionation factors."""

    deuterium: FloatResult
    oxygen_18: FloatResult
    oxygen_17: FloatResult


def equilibrium_fractionation_factors(
    temperature_c: ArrayLike, *, oxygen_17_ice_exponent: float
) -> EquilibriumFractionationFactors:
    """Calculate the active equilibrium fractionation factors.

    Temperatures are in degrees Celsius.  The 17O ice exponent is deliberately
    required because the legacy evaporation path uses 0.529 while the active
    distillation path uses 0.531.  Requiring it prevents that unresolved
    scientific discrepancy from being hidden behind a default.
    """
    temperature_k = np.asarray(temperature_c, dtype=np.float64) + CELSIUS_TO_KELVIN

    deuterium_liquid = np.exp(
        (
            52.612
            - 76.248e3 / temperature_k
            + 24.844e6 / temperature_k**2
        )
        / 1000.0
    )
    deuterium_ice = np.exp(-0.0559 + 13525.0 / temperature_k**2)
    oxygen_18_liquid = np.exp(
        (-2.0667 - 0.4156e3 / temperature_k + 1.137e6 / temperature_k**2)
        / 1000.0
    )
    oxygen_18_ice = np.exp(
        (-28.224 + 11.839e3 / temperature_k) / 1000.0
    )
    oxygen_17_liquid = oxygen_18_liquid**O17_LIQUID_EXPONENT
    oxygen_17_ice = oxygen_18_ice**oxygen_17_ice_exponent

    return EquilibriumFractionationFactors(
        deuterium_liquid=deuterium_liquid,
        deuterium_ice=deuterium_ice,
        oxygen_18_liquid=oxygen_18_liquid,
        oxygen_18_ice=oxygen_18_ice,
        oxygen_17_liquid=oxygen_17_liquid,
        oxygen_17_ice=oxygen_17_ice,
    )


def transport_diffusivity_ratios(
    temperature_c: ArrayLike,
) -> TransportDiffusivityRatios:
    """Return active heavy/light diffusivity ratios for water in air.

    These are the temperature-dependent Hellmann--Harvey ratios assigned to
    ``D_r_HDO``, ``D_r_17``, and ``D_r_18`` in ``distillation_2020.m``.
    """
    temperature_k = np.asarray(temperature_c, dtype=np.float64) + CELSIUS_TO_KELVIN
    scaled_temperature = temperature_k / 100.0
    hdo_over_h2o = (
        0.98258
        - 0.02546 / scaled_temperature
        + 0.02421 / scaled_temperature ** (5.0 / 2.0)
    )
    h217o_over_h216o = (
        0.98284
        + 0.003517 / scaled_temperature ** (1.0 / 2.0)
        - 0.001996 / scaled_temperature ** (5.0 / 2.0)
    )
    h218o_over_h216o = (
        0.96671
        + 0.007406 / scaled_temperature ** (1.0 / 2.0)
        - 0.004861 / scaled_temperature**3
    )
    return TransportDiffusivityRatios(
        hdo_over_h2o=hdo_over_h2o,
        h217o_over_h216o=h217o_over_h216o,
        h218o_over_h216o=h218o_over_h216o,
    )


def kinetic_condensation_factors(
    saturation_over_ice: ArrayLike,
    equilibrium: EquilibriumFractionationFactors,
    diffusivity: TransportDiffusivityRatios,
) -> KineticCondensationFactors:
    """Return ice-vapor kinetic factors using paper Eq. (A12).

    The liquid kinetic factors in the selected legacy pathway are exactly one,
    so only the ice factors are represented here.
    """
    saturation = np.asarray(saturation_over_ice, dtype=np.float64)

    def ice_factor(equilibrium_factor: FloatResult, heavy_over_light: FloatResult):
        equilibrium_array = np.asarray(equilibrium_factor, dtype=np.float64)
        heavy_over_light_array = np.asarray(heavy_over_light, dtype=np.float64)
        light_over_heavy = 1.0 / heavy_over_light_array
        return saturation / (
            equilibrium_array * light_over_heavy * (saturation - 1.0) + 1.0
        )

    return KineticCondensationFactors(
        deuterium_ice=ice_factor(
            equilibrium.deuterium_ice, diffusivity.hdo_over_h2o
        ),
        oxygen_18_ice=ice_factor(
            equilibrium.oxygen_18_ice, diffusivity.h218o_over_h216o
        ),
        oxygen_17_ice=ice_factor(
            equilibrium.oxygen_17_ice, diffusivity.h217o_over_h216o
        ),
    )


def mixed_phase_effective_fractionation(
    equilibrium: EquilibriumFractionationFactors,
    kinetic: KineticCondensationFactors,
    fraction_ice: ArrayLike,
    fraction_liquid: ArrayLike,
) -> EffectiveFractionationFactors:
    """Combine ice and liquid total fractionation using paper Eq. (A13)."""
    ice = np.asarray(fraction_ice, dtype=np.float64)
    liquid = np.asarray(fraction_liquid, dtype=np.float64)

    return EffectiveFractionationFactors(
        deuterium=(
            np.asarray(equilibrium.deuterium_ice) * kinetic.deuterium_ice * ice
            + np.asarray(equilibrium.deuterium_liquid) * liquid
        ),
        oxygen_18=(
            np.asarray(equilibrium.oxygen_18_ice) * kinetic.oxygen_18_ice * ice
            + np.asarray(equilibrium.oxygen_18_liquid) * liquid
        ),
        oxygen_17=(
            np.asarray(equilibrium.oxygen_17_ice) * kinetic.oxygen_17_ice * ice
            + np.asarray(equilibrium.oxygen_17_liquid) * liquid
        ),
    )


__all__ = [
    "O17_ICE_EXPONENT_DISTILLATION",
    "O17_ICE_EXPONENT_EVAPORATION",
    "O17_LIQUID_EXPONENT",
    "EffectiveFractionationFactors",
    "EquilibriumFractionationFactors",
    "KineticCondensationFactors",
    "TransportDiffusivityRatios",
    "equilibrium_fractionation_factors",
    "kinetic_condensation_factors",
    "mixed_phase_effective_fractionation",
    "transport_diffusivity_ratios",
]

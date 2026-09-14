from pathlib import Path

import numpy as np
import pytest

from swim.cloud_phase import cloud_phase_fractions
from swim.fractionation import (
    O17_ICE_EXPONENT_DISTILLATION,
    equilibrium_fractionation_factors,
    kinetic_condensation_factors,
    mixed_phase_effective_fractionation,
    transport_diffusivity_ratios,
)

TRAJECTORY_PATH = (
    Path(__file__).parents[1]
    / "fixtures"
    / "matlab"
    / "port_baseline_v1"
    / "allan_hills"
    / "trajectory_Tsource_10_Tcond_m30.csv"
)


@pytest.mark.parity
def test_equilibrium_factors_match_matlab_where_kinetic_factors_are_one() -> None:
    trajectory = np.genfromtxt(TRAJECTORY_PATH, delimiter=",", names=True)
    # MATLAB initializes all kinetic factors to one at the first point. At
    # 0°C its prescribed saturation is one, which also makes them exactly one.
    indices = np.array([0, np.flatnonzero(trajectory["temperature_C"] == 0.0)[0]])
    temperature_c = trajectory["temperature_C"][indices]
    fraction_ice, fraction_liquid = cloud_phase_fractions(
        temperature_c, method="adj"
    )
    factors = equilibrium_fractionation_factors(
        temperature_c,
        oxygen_17_ice_exponent=O17_ICE_EXPONENT_DISTILLATION,
    )

    effective_deuterium = (
        fraction_ice * factors.deuterium_ice
        + fraction_liquid * factors.deuterium_liquid
    )
    effective_oxygen_18 = (
        fraction_ice * factors.oxygen_18_ice
        + fraction_liquid * factors.oxygen_18_liquid
    )
    effective_oxygen_17 = (
        fraction_ice * factors.oxygen_17_ice
        + fraction_liquid * factors.oxygen_17_liquid
    )

    np.testing.assert_allclose(
        effective_deuterium,
        trajectory["alpha_D_effective"][indices],
        rtol=1.0e-13,
        atol=1.0e-14,
    )
    np.testing.assert_allclose(
        effective_oxygen_18,
        trajectory["alpha_18O_effective"][indices],
        rtol=1.0e-13,
        atol=1.0e-14,
    )
    np.testing.assert_allclose(
        effective_oxygen_17,
        trajectory["alpha_17O_effective"][indices],
        rtol=1.0e-13,
        atol=1.0e-14,
    )


@pytest.mark.parity
def test_effective_fractionation_matches_full_matlab_trajectory() -> None:
    trajectory = np.genfromtxt(TRAJECTORY_PATH, delimiter=",", names=True)
    temperature_c = trajectory["temperature_C"]
    fraction_ice, fraction_liquid = cloud_phase_fractions(
        temperature_c, method="adj"
    )
    equilibrium = equilibrium_fractionation_factors(
        temperature_c,
        oxygen_17_ice_exponent=O17_ICE_EXPONENT_DISTILLATION,
    )
    diffusivity = transport_diffusivity_ratios(temperature_c)
    kinetic = kinetic_condensation_factors(
        trajectory["saturation_used"], equilibrium, diffusivity
    )
    effective = mixed_phase_effective_fractionation(
        equilibrium, kinetic, fraction_ice, fraction_liquid
    )

    # The first precipitation point is initialized separately with kinetic
    # factors of one in MATLAB. The vector formulas govern all later points.
    comparison = slice(1, None)
    np.testing.assert_allclose(
        np.asarray(effective.deuterium)[comparison],
        trajectory["alpha_D_effective"][comparison],
        rtol=1.0e-12,
        atol=1.0e-14,
    )
    np.testing.assert_allclose(
        np.asarray(effective.oxygen_18)[comparison],
        trajectory["alpha_18O_effective"][comparison],
        rtol=1.0e-12,
        atol=1.0e-14,
    )
    np.testing.assert_allclose(
        np.asarray(effective.oxygen_17)[comparison],
        trajectory["alpha_17O_effective"][comparison],
        rtol=1.0e-12,
        atol=1.0e-14,
    )

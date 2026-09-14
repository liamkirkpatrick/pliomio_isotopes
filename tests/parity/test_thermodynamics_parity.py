from pathlib import Path

import numpy as np
import pytest

from swim.cloud_phase import cloud_phase_fractions
from swim.thermodynamics import mixed_phase_supersaturation, pseudo_adiabat

TRAJECTORY_PATH = (
    Path(__file__).parents[1]
    / "fixtures"
    / "matlab"
    / "port_baseline_v1"
    / "allan_hills"
    / "trajectory_Tsource_10_Tcond_m30.csv"
)


@pytest.mark.parity
def test_pseudo_adiabat_matches_full_matlab_trajectory() -> None:
    trajectory = np.genfromtxt(TRAJECTORY_PATH, delimiter=",", names=True)
    temperature_c = trajectory["temperature_C"]
    fraction_ice, fraction_liquid = cloud_phase_fractions(
        temperature_c, method="adj"
    )

    mixed_saturation = mixed_phase_supersaturation(
        temperature_c,
        initial_pressure_kpa=101.325,
        fraction_ice=fraction_ice,
        fraction_liquid=fraction_liquid,
    )
    result = pseudo_adiabat(
        temperature_c,
        initial_pressure_kpa=101.325,
        fraction_ice=fraction_ice,
        fraction_liquid=fraction_liquid,
        mixed_phase_saturation=mixed_saturation,
        supersaturation_a=1.0,
        supersaturation_b=0.00525,
        supersaturation_c=0.0,
    )

    np.testing.assert_allclose(
        mixed_saturation,
        trajectory["mixed_phase_supersaturation"],
        rtol=1.0e-12,
        atol=1.0e-14,
    )
    np.testing.assert_allclose(
        result.pressure_kpa, trajectory["pressure_kPa"], rtol=1.0e-12, atol=1e-11
    )
    np.testing.assert_allclose(
        result.saturation_vapor_pressure_kpa,
        trajectory["saturation_vapor_pressure_kPa"],
        rtol=1.0e-13,
        atol=1.0e-14,
    )
    np.testing.assert_allclose(
        result.saturated_mixing_ratio,
        trajectory["saturated_mixing_ratio"],
        rtol=1.0e-12,
        atol=1.0e-14,
    )
    np.testing.assert_allclose(
        result.fraction_vapor_remaining,
        trajectory["fraction_vapor_remaining"],
        rtol=1.0e-12,
        atol=1.0e-14,
    )
    np.testing.assert_allclose(
        result.saturation_over_ice,
        trajectory["saturation_used"],
        rtol=1.0e-12,
        atol=1.0e-14,
    )

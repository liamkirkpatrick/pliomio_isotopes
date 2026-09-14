from pathlib import Path

import numpy as np
import pytest

from swim.cloud_phase import cloud_phase_fractions
from swim.isotopes import (
    R17O_VSMOW,
    R18O_VSMOW,
    RD_VSMOW,
    linear_deuterium_excess,
    log_delta,
    logarithmic_deuterium_excess,
    oxygen_17_excess,
    ratio_to_delta,
)
from swim.saturation import (
    saturated_mixing_ratio,
    saturation_vapor_pressure_ice,
    saturation_vapor_pressure_liquid,
)

TRAJECTORY_PATH = (
    Path(__file__).parents[1]
    / "fixtures"
    / "matlab"
    / "port_baseline_v1"
    / "allan_hills"
    / "trajectory_Tsource_10_Tcond_m30.csv"
)


def load_trajectory() -> np.ndarray:
    return np.genfromtxt(TRAJECTORY_PATH, delimiter=",", names=True)


@pytest.mark.parity
def test_isotope_helpers_match_matlab_trajectory() -> None:
    trajectory = load_trajectory()

    np.testing.assert_allclose(
        ratio_to_delta(trajectory["R18O_precip"], R18O_VSMOW),
        trajectory["d18O_precip"],
        rtol=1.0e-13,
        atol=4.0e-12,
        equal_nan=True,
    )
    np.testing.assert_allclose(
        ratio_to_delta(trajectory["RD_precip"], RD_VSMOW),
        trajectory["dD_precip"],
        rtol=1.0e-13,
        atol=4.0e-12,
        equal_nan=True,
    )
    np.testing.assert_allclose(
        ratio_to_delta(trajectory["R17O_precip"], R17O_VSMOW),
        trajectory["d17O_precip"],
        rtol=1.0e-13,
        atol=4.0e-12,
        equal_nan=True,
    )
    np.testing.assert_allclose(
        log_delta(trajectory["d18O_precip"]),
        trajectory["d18O_ln"],
        rtol=1.0e-13,
        atol=1.0e-12,
        equal_nan=True,
    )
    np.testing.assert_allclose(
        log_delta(trajectory["dD_precip"]),
        trajectory["dD_ln"],
        rtol=1.0e-13,
        atol=1.0e-12,
        equal_nan=True,
    )
    np.testing.assert_allclose(
        linear_deuterium_excess(
            trajectory["dD_precip"], trajectory["d18O_precip"]
        ),
        trajectory["dxs"],
        rtol=1.0e-13,
        atol=1.0e-12,
        equal_nan=True,
    )
    np.testing.assert_allclose(
        logarithmic_deuterium_excess(
            trajectory["dD_precip"], trajectory["d18O_precip"]
        ),
        trajectory["dln"],
        rtol=1.0e-13,
        atol=1.0e-12,
        equal_nan=True,
    )
    np.testing.assert_allclose(
        oxygen_17_excess(
            trajectory["d17O_precip"], trajectory["d18O_precip"]
        ),
        trajectory["d17O_excess"],
        rtol=1.0e-12,
        atol=1.0e-8,
        equal_nan=True,
    )


@pytest.mark.parity
def test_saturation_pressure_matches_matlab_trajectory() -> None:
    trajectory = load_trajectory()
    liquid_kpa = saturation_vapor_pressure_liquid(trajectory["temperature_C"])
    ice_kpa = saturation_vapor_pressure_ice(trajectory["temperature_C"])
    mixed_kpa = (
        trajectory["fraction_liquid"] * liquid_kpa
        + trajectory["fraction_ice"] * ice_kpa
    )

    np.testing.assert_allclose(
        mixed_kpa,
        trajectory["saturation_vapor_pressure_kPa"],
        rtol=1.0e-13,
        atol=1.0e-14,
    )


@pytest.mark.parity
def test_active_cloud_phase_matches_matlab_trajectory() -> None:
    trajectory = load_trajectory()

    fraction_ice, fraction_liquid = cloud_phase_fractions(
        trajectory["temperature_C"], method="adj"
    )

    np.testing.assert_allclose(
        fraction_ice, trajectory["fraction_ice"], rtol=1.0e-12, atol=1.0e-15
    )
    np.testing.assert_allclose(
        fraction_liquid,
        trajectory["fraction_liquid"],
        rtol=2.0e-14,
        atol=1.0e-14,
    )


@pytest.mark.parity
def test_saturated_mixing_ratio_matches_matlab_trajectory() -> None:
    trajectory = load_trajectory()
    ice_vapor_pressure_kpa = saturation_vapor_pressure_ice(
        trajectory["temperature_C"]
    )
    ice_mixing_ratio = saturated_mixing_ratio(
        ice_vapor_pressure_kpa, trajectory["pressure_kPa"]
    )

    # pseudo_adiabat_function.m defines saturation_used as r_s / r_s_i.
    expected = trajectory["saturation_used"] * ice_mixing_ratio
    np.testing.assert_allclose(
        expected,
        trajectory["saturated_mixing_ratio"],
        rtol=1.0e-13,
        atol=1.0e-15,
    )

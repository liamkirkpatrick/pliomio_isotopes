from pathlib import Path

import numpy as np
import pytest

from swim.model import forward_trajectory

TRAJECTORY_PATH = (
    Path(__file__).parents[1]
    / "fixtures"
    / "matlab"
    / "port_baseline_v1"
    / "allan_hills"
    / "trajectory_Tsource_10_Tcond_m30.csv"
)


@pytest.mark.parity
def test_composed_forward_trajectory_matches_matlab() -> None:
    expected = np.genfromtxt(TRAJECTORY_PATH, delimiter=",", names=True)
    actual = forward_trajectory(10.0, -30.0)

    comparisons = {
        "temperature_C": actual.temperature_c,
        "pressure_kPa": actual.thermodynamics.pressure_kpa,
        "saturation_vapor_pressure_kPa": (
            actual.thermodynamics.saturation_vapor_pressure_kpa
        ),
        "saturated_mixing_ratio": actual.thermodynamics.saturated_mixing_ratio,
        "fraction_vapor_remaining": (
            actual.thermodynamics.fraction_vapor_remaining
        ),
        "saturation_used": actual.thermodynamics.saturation_over_ice,
        "mixed_phase_supersaturation": actual.mixed_phase_supersaturation,
        "fraction_ice": actual.fraction_ice,
        "fraction_liquid": actual.fraction_liquid,
        "d18O_precip": actual.distillation.delta_18o_precipitation,
        "dD_precip": actual.distillation.delta_d_precipitation,
        "d17O_precip": actual.distillation.delta_17o_precipitation,
        "d18O_ln": actual.distillation.delta_18o_log,
        "dD_ln": actual.distillation.delta_d_log,
        "d17O_ln": actual.distillation.delta_17o_log,
        "dxs": actual.distillation.deuterium_excess,
        "dln": actual.distillation.logarithmic_deuterium_excess,
        "d17O_excess": actual.distillation.oxygen_17_excess_per_meg,
        "RD_vapor": actual.distillation.deuterium_vapor_ratio,
        "RD_precip": actual.distillation.deuterium_precipitation_ratio,
        "R18O_vapor": actual.distillation.oxygen_18_vapor_ratio,
        "R18O_precip": actual.distillation.oxygen_18_precipitation_ratio,
        "R17O_vapor": actual.distillation.oxygen_17_vapor_ratio,
        "R17O_precip": actual.distillation.oxygen_17_precipitation_ratio,
        "alpha_D_effective": actual.effective_fractionation.deuterium,
        "alpha_18O_effective": actual.effective_fractionation.oxygen_18,
        "alpha_17O_effective": actual.effective_fractionation.oxygen_17,
    }
    for matlab_column, python_values in comparisons.items():
        absolute_tolerance = 5.0e-9 if matlab_column == "d17O_excess" else 1.0e-10
        np.testing.assert_allclose(
            python_values,
            expected[matlab_column],
            rtol=1.0e-10,
            atol=absolute_tolerance,
            equal_nan=True,
            err_msg=matlab_column,
        )

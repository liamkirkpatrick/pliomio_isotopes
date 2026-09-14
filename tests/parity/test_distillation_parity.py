import json
from pathlib import Path

import numpy as np
import pytest

from swim.cloud_phase import cloud_phase_fractions
from swim.distillation import rayleigh_distillation
from swim.fractionation import (
    O17_ICE_EXPONENT_DISTILLATION,
    EffectiveFractionationFactors,
    equilibrium_fractionation_factors,
    kinetic_condensation_factors,
    mixed_phase_effective_fractionation,
    transport_diffusivity_ratios,
)

FIXTURE_DIR = (
    Path(__file__).parents[1]
    / "fixtures"
    / "matlab"
    / "port_baseline_v1"
    / "allan_hills"
)


@pytest.mark.parity
def test_rayleigh_distillation_matches_full_matlab_trajectory() -> None:
    trajectory = np.genfromtxt(
        FIXTURE_DIR / "trajectory_Tsource_10_Tcond_m30.csv",
        delimiter=",",
        names=True,
    )
    metadata = json.loads((FIXTURE_DIR / "metadata.json").read_text())
    initial = metadata["result"]["trajectory_initial_conditions"]
    temperature_c = trajectory["temperature_C"]
    fraction_ice, fraction_liquid = cloud_phase_fractions(
        temperature_c, method="adj"
    )
    equilibrium = equilibrium_fractionation_factors(
        temperature_c,
        oxygen_17_ice_exponent=O17_ICE_EXPONENT_DISTILLATION,
    )
    kinetic = kinetic_condensation_factors(
        trajectory["saturation_used"],
        equilibrium,
        transport_diffusivity_ratios(temperature_c),
    )
    calculated = mixed_phase_effective_fractionation(
        equilibrium, kinetic, fraction_ice, fraction_liquid
    )

    # MATLAB calculates the first precipitation point before replacing its
    # scalar kinetic initialization with the vector kinetic formulas.
    alpha_d = np.asarray(calculated.deuterium).copy()
    alpha_18o = np.asarray(calculated.oxygen_18).copy()
    alpha_17o = np.asarray(calculated.oxygen_17).copy()
    alpha_d[0] = (
        fraction_ice[0] * np.asarray(equilibrium.deuterium_ice)[0]
        + fraction_liquid[0] * np.asarray(equilibrium.deuterium_liquid)[0]
    )
    alpha_18o[0] = (
        fraction_ice[0] * np.asarray(equilibrium.oxygen_18_ice)[0]
        + fraction_liquid[0] * np.asarray(equilibrium.oxygen_18_liquid)[0]
    )
    alpha_17o[0] = (
        fraction_ice[0] * np.asarray(equilibrium.oxygen_17_ice)[0]
        + fraction_liquid[0] * np.asarray(equilibrium.oxygen_17_liquid)[0]
    )
    effective = EffectiveFractionationFactors(alpha_d, alpha_18o, alpha_17o)

    result = rayleigh_distillation(
        trajectory["fraction_vapor_remaining"],
        effective,
        initial_delta_d_permil=initial["dD_v0"],
        initial_delta_18o_permil=initial["d18O_v0"],
        initial_oxygen_17_excess_log=initial["d17Oxs_v0"],
    )

    comparisons = {
        "RD_vapor": result.deuterium_vapor_ratio,
        "RD_precip": result.deuterium_precipitation_ratio,
        "R18O_vapor": result.oxygen_18_vapor_ratio,
        "R18O_precip": result.oxygen_18_precipitation_ratio,
        "R17O_vapor": result.oxygen_17_vapor_ratio,
        "R17O_precip": result.oxygen_17_precipitation_ratio,
        "dD_precip": result.delta_d_precipitation,
        "d18O_precip": result.delta_18o_precipitation,
        "d17O_precip": result.delta_17o_precipitation,
        "dD_ln": result.delta_d_log,
        "d18O_ln": result.delta_18o_log,
        "d17O_ln": result.delta_17o_log,
        "dxs": result.deuterium_excess,
        "dln": result.logarithmic_deuterium_excess,
        "d17O_excess": result.oxygen_17_excess_per_meg,
    }
    for matlab_column, python_values in comparisons.items():
        absolute_tolerance = 4.0e-9 if matlab_column == "d17O_excess" else 1.0e-10
        np.testing.assert_allclose(
            python_values,
            trajectory[matlab_column],
            rtol=1.0e-11,
            atol=absolute_tolerance,
            equal_nan=True,
            err_msg=matlab_column,
        )

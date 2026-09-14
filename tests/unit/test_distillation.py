import numpy as np

from swim.distillation import rayleigh_distillation
from swim.fractionation import EffectiveFractionationFactors


def test_zero_vapor_loss_produces_no_new_precipitation() -> None:
    fraction = np.array([1.0, 1.0, 0.9])
    effective = EffectiveFractionationFactors(
        deuterium=np.array([1.1, 1.1, 1.1]),
        oxygen_18=np.array([1.01, 1.01, 1.01]),
        oxygen_17=np.array([1.005, 1.005, 1.005]),
    )

    result = rayleigh_distillation(
        fraction,
        effective,
        initial_delta_d_permil=-100.0,
        initial_delta_18o_permil=-15.0,
        initial_oxygen_17_excess_log=5.0e-6,
    )

    assert np.isnan(result.deuterium_precipitation_ratio[1])
    assert np.isnan(result.oxygen_18_precipitation_ratio[1])
    assert np.isnan(result.oxygen_17_precipitation_ratio[1])
    assert np.isfinite(result.deuterium_precipitation_ratio[2])


def test_vapor_ratios_decrease_during_rayleigh_distillation() -> None:
    fraction = np.array([1.0, 0.9, 0.8])
    effective = EffectiveFractionationFactors(
        deuterium=np.full(3, 1.1),
        oxygen_18=np.full(3, 1.01),
        oxygen_17=np.full(3, 1.005),
    )
    result = rayleigh_distillation(
        fraction,
        effective,
        initial_delta_d_permil=-100.0,
        initial_delta_18o_permil=-15.0,
        initial_oxygen_17_excess_log=5.0e-6,
    )

    assert np.all(np.diff(result.deuterium_vapor_ratio) < 0.0)
    assert np.all(np.diff(result.oxygen_18_vapor_ratio) < 0.0)
    assert np.all(np.diff(result.oxygen_17_vapor_ratio) < 0.0)

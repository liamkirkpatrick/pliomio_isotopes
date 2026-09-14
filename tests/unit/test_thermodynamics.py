import numpy as np
import pytest

from swim.cloud_phase import cloud_phase_fractions
from swim.thermodynamics import (
    mixed_phase_supersaturation,
    prescribed_supersaturation,
    pseudo_adiabat,
)


def test_prescribed_supersaturation_clips_values_below_one() -> None:
    temperature_c = np.array([10.0, 0.0, -10.0, -30.0])

    result = prescribed_supersaturation(temperature_c, 1.0, 0.00525, 0.0)

    np.testing.assert_allclose(result, np.array([1.0, 1.0, 1.0525, 1.1575]))


def test_short_pseudo_adiabat_has_expected_physical_direction() -> None:
    temperature_c = np.array([10.0, 5.0, 0.0, -5.0])
    fraction_ice, fraction_liquid = cloud_phase_fractions(
        temperature_c, method="adj"
    )
    mixed_saturation = mixed_phase_supersaturation(
        temperature_c, 101.325, fraction_ice, fraction_liquid
    )

    result = pseudo_adiabat(
        temperature_c,
        101.325,
        fraction_ice,
        fraction_liquid,
        mixed_saturation,
        supersaturation_a=1.0,
        supersaturation_b=0.00525,
        supersaturation_c=0.0,
    )

    assert result.pressure_kpa[0] == 101.325
    assert result.fraction_vapor_remaining[0] == 1.0
    assert np.all(np.diff(result.pressure_kpa) < 0.0)
    assert np.all(np.diff(result.fraction_vapor_remaining) < 0.0)


def test_trajectory_inputs_must_have_matching_shapes() -> None:
    with pytest.raises(ValueError, match="equal shapes"):
        mixed_phase_supersaturation([10.0, 9.0], 101.325, [0.0], [1.0])

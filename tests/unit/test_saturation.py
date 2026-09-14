import numpy as np

from swim.saturation import (
    WATER_VAPOR_TO_DRY_AIR_GAS_CONSTANT_RATIO,
    saturated_mixing_ratio,
    saturation_vapor_pressure_ice,
    saturation_vapor_pressure_liquid,
)


def test_saturation_vapor_pressure_reference_points() -> None:
    temperature_c = np.array([0.0, -20.0, -40.0])

    liquid_kpa = saturation_vapor_pressure_liquid(temperature_c)
    ice_kpa = saturation_vapor_pressure_ice(temperature_c)

    np.testing.assert_allclose(
        liquid_kpa,
        np.array(
            [0.6112126978267947, 0.12550416935494008, 0.018912149430063536]
        ),
        rtol=1.0e-14,
    )
    np.testing.assert_allclose(
        ice_kpa,
        np.array(
            [0.6111535914292243, 0.1032524632801715, 0.012844281376147109]
        ),
        rtol=1.0e-14,
    )


def test_saturation_vapor_pressure_decreases_with_temperature() -> None:
    temperature_c = np.array([-40.0, -20.0, 0.0, 20.0])

    assert np.all(np.diff(saturation_vapor_pressure_liquid(temperature_c)) > 0.0)
    assert np.all(np.diff(saturation_vapor_pressure_ice(temperature_c)) > 0.0)


def test_saturated_mixing_ratio_matches_legacy_equation() -> None:
    vapor_pressure_kpa = np.array([1.0, 0.5])
    pressure_kpa = np.array([100.0, 80.0])

    expected = (
        WATER_VAPOR_TO_DRY_AIR_GAS_CONSTANT_RATIO
        * vapor_pressure_kpa
        / (pressure_kpa - vapor_pressure_kpa)
    )

    np.testing.assert_allclose(
        saturated_mixing_ratio(vapor_pressure_kpa, pressure_kpa), expected
    )

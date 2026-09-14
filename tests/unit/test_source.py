import numpy as np

from swim.source import (
    climatological_source_conditions,
    initial_vapor_from_climatology,
    seawater_delta_d_from_delta_18o,
)


def test_climatological_source_conditions_are_vectorized() -> None:
    source = climatological_source_conditions(np.array([0.0, 10.0, 20.0]))

    assert np.asarray(source.relative_humidity).shape == (3,)
    assert np.asarray(source.sea_surface_temperature_c).shape == (3,)
    assert np.all(np.isfinite(source.normalized_relative_humidity))


def test_corrected_source_conditions_are_the_default() -> None:
    corrected = climatological_source_conditions(10.0)
    legacy = climatological_source_conditions(10.0, version="2020")

    np.testing.assert_allclose(
        corrected.normalized_relative_humidity, 0.8025892975824845, atol=1e-14
    )
    np.testing.assert_allclose(
        legacy.normalized_relative_humidity, 0.8185693991506341, atol=1e-14
    )


def test_seawater_deuterium_fit_is_linear() -> None:
    delta_18o = np.array([-0.3, 0.0, 0.3])
    delta_d = seawater_delta_d_from_delta_18o(delta_18o)

    np.testing.assert_allclose(np.diff(delta_d), np.diff(delta_d)[0])


def test_global_and_local_evaporation_closures_are_available() -> None:
    local = initial_vapor_from_climatology(10.0, closure="local")
    global_result = initial_vapor_from_climatology(10.0, closure="global")

    assert np.isfinite(global_result.delta_d_permil)
    assert np.isfinite(global_result.delta_18o_permil)
    assert global_result.delta_d_permil != local.delta_d_permil


def test_corrected_2022_evaporation_is_the_default() -> None:
    corrected = initial_vapor_from_climatology(10.0)
    legacy = initial_vapor_from_climatology(10.0, evaporation_version="2021")

    np.testing.assert_allclose(corrected.delta_d_permil, -92.63024337997217)
    np.testing.assert_allclose(corrected.delta_18o_permil, -12.448385650174366)
    assert corrected.delta_d_permil != legacy.delta_d_permil
    assert corrected.delta_18o_permil != legacy.delta_18o_permil

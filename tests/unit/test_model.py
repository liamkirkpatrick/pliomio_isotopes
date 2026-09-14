import numpy as np

from swim.model import forward_trajectory, matlab_temperature_grid


def test_matlab_temperature_grid_is_descending_and_inclusive() -> None:
    result = matlab_temperature_grid(1.0, 0.0, 0.1)

    assert result.shape == (11,)
    assert result[0] == 1.0
    assert result[-1] == 0.0
    np.testing.assert_allclose(np.diff(result), -0.1)


def test_single_point_trajectory_matches_legacy_endpoint_branch() -> None:
    result = forward_trajectory(0.0, 0.0)

    assert result.temperature_c.shape == (1,)
    assert result.thermodynamics.pressure_kpa[0] == 101.325
    assert result.thermodynamics.fraction_vapor_remaining[0] == 1.0


def test_forward_model_defaults_to_corrected_2022_evaporation() -> None:
    corrected = forward_trajectory(10.0, 9.9)
    legacy = forward_trajectory(10.0, 9.9, evaporation_version="2021")

    np.testing.assert_allclose(
        corrected.source_vapor.delta_18o_permil, -12.448385650174366
    )
    assert (
        corrected.distillation.delta_18o_precipitation[-1]
        != legacy.distillation.delta_18o_precipitation[-1]
    )

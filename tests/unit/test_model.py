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

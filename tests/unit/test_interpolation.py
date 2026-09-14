import numpy as np

from swim.interpolation import natural_neighbor_interpolate


def test_natural_neighbor_reproduces_a_linear_plane() -> None:
    points = np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]])
    values = 2.0 * points[:, 0] - 3.0 * points[:, 1] + 4.0
    query = np.array([[0.25, 0.5], [0.75, 0.25]])

    actual = natural_neighbor_interpolate(points, values, query)

    np.testing.assert_allclose(actual, 2.0 * query[:, 0] - 3.0 * query[:, 1] + 4.0)


def test_natural_neighbor_returns_nan_outside_convex_hull() -> None:
    points = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
    result = natural_neighbor_interpolate(points, np.arange(3.0), [[2.0, 2.0]])

    assert np.isnan(result[0])

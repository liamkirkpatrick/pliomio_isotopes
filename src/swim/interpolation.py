"""Two-dimensional Sibson natural-neighbor interpolation."""

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.spatial import ConvexHull, Delaunay, QhullError


def _circumcenter(
    first: NDArray[np.float64],
    second: NDArray[np.float64],
    third: NDArray[np.float64],
) -> NDArray[np.float64]:
    ax, ay = first
    bx, by = second
    cx, cy = third
    denominator = 2.0 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
    return np.array(
        [
            (
                (ax**2 + ay**2) * (by - cy)
                + (bx**2 + by**2) * (cy - ay)
                + (cx**2 + cy**2) * (ay - by)
            )
            / denominator,
            (
                (ax**2 + ay**2) * (cx - bx)
                + (bx**2 + by**2) * (ax - cx)
                + (cx**2 + cy**2) * (bx - ax)
            )
            / denominator,
        ]
    )


def _polygon_area(vertices: NDArray[np.float64]) -> float:
    return float(
        0.5
        * abs(
            np.dot(vertices[:, 0], np.roll(vertices[:, 1], -1))
            - np.dot(vertices[:, 1], np.roll(vertices[:, 0], -1))
        )
    )


def _ordered_boundary_vertices(
    triangulation: Delaunay, neighboring_simplices: NDArray[np.int64]
) -> list[int]:
    edge_counts: dict[tuple[int, int], int] = {}
    for simplex_index in neighboring_simplices:
        p0, p1, p2 = triangulation.simplices[simplex_index]
        for start, end in ((p0, p1), (p1, p2), (p2, p0)):
            start_index = int(start)
            end_index = int(end)
            edge = (min(start_index, end_index), max(start_index, end_index))
            edge_counts[edge] = edge_counts.get(edge, 0) + 1
    boundary_edges = [edge for edge, count in edge_counts.items() if count == 1]
    if not boundary_edges:
        return []

    adjacency: dict[int, list[int]] = {}
    for start, end in boundary_edges:
        adjacency.setdefault(start, []).append(end)
        adjacency.setdefault(end, []).append(start)
    if any(len(neighbors) != 2 for neighbors in adjacency.values()):
        return []

    start = boundary_edges[0][0]
    ordered = [start]
    previous: int | None = None
    current = start
    while True:
        candidates = adjacency[current]
        following = candidates[0] if candidates[0] != previous else candidates[1]
        if following == start:
            break
        ordered.append(following)
        previous, current = current, following
        if len(ordered) > len(adjacency):
            return []
    return ordered


def natural_neighbor_interpolate(
    points: ArrayLike, values: ArrayLike, query_points: ArrayLike
) -> NDArray[np.float64]:
    """Interpolate values using Sibson natural-neighbor area weights.

    The construction follows the circumcircle/overlap-area approach of Liang
    and Hale (2010). Queries outside the convex hull return ``NaN``, matching
    MATLAB ``griddata(..., 'natural')``.
    """
    sample_points = np.asarray(points, dtype=np.float64)
    sample_values = np.asarray(values, dtype=np.float64)
    queries = np.asarray(query_points, dtype=np.float64)
    if sample_points.ndim != 2 or sample_points.shape[1] != 2:
        raise ValueError("points must have shape (n, 2)")
    if queries.ndim != 2 or queries.shape[1] != 2:
        raise ValueError("query_points must have shape (m, 2)")
    if sample_values.shape[0] != sample_points.shape[0]:
        raise ValueError("values and points must have the same first dimension")

    scalar_values = sample_values.ndim == 1
    values_2d = sample_values[:, None] if scalar_values else sample_values
    triangulation = Delaunay(sample_points)
    circumcenters = np.array(
        [_circumcenter(*sample_points[simplex]) for simplex in triangulation.simplices]
    )
    radius_squared = np.sum(
        (circumcenters - sample_points[triangulation.simplices[:, 0]]) ** 2,
        axis=1,
    )
    result = np.full((queries.shape[0], values_2d.shape[1]), np.nan)

    for query_index, query in enumerate(queries):
        if triangulation.find_simplex(query) < 0:
            continue
        exact = np.flatnonzero(np.all(sample_points == query, axis=1))
        if exact.size:
            result[query_index] = values_2d[exact[0]]
            continue

        distance_squared = np.sum((circumcenters - query) ** 2, axis=1)
        tolerance = 1.0e-12 * np.maximum(radius_squared, 1.0)
        neighboring_simplices = np.flatnonzero(
            distance_squared <= radius_squared + tolerance
        )
        boundary = _ordered_boundary_vertices(triangulation, neighboring_simplices)
        if len(boundary) < 3:
            continue

        neighbor_indices: list[int] = []
        areas: list[float] = []
        for boundary_index, middle_index in enumerate(boundary):
            previous_index = boundary[boundary_index - 1]
            following_index = boundary[(boundary_index + 1) % len(boundary)]
            try:
                polygon = [
                    _circumcenter(
                        query,
                        sample_points[previous_index],
                        sample_points[middle_index],
                    ),
                    _circumcenter(
                        query,
                        sample_points[following_index],
                        sample_points[middle_index],
                    ),
                ]
                polygon.extend(
                    circumcenters[simplex_index]
                    for simplex_index in neighboring_simplices
                    if middle_index in triangulation.simplices[simplex_index]
                )
                polygon_array = np.asarray(polygon)
                hull = ConvexHull(polygon_array)
                area = _polygon_area(polygon_array[hull.vertices])
            except (ZeroDivisionError, QhullError, FloatingPointError):
                area = 0.0
            if area > 0.0:
                neighbor_indices.append(middle_index)
                areas.append(area)

        total_area = sum(areas)
        if total_area > 0.0:
            weights = np.asarray(areas) / total_area
            result[query_index] = weights @ values_2d[neighbor_indices]

    return result[:, 0] if scalar_values else result


__all__ = ["natural_neighbor_interpolate"]

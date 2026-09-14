import json
from pathlib import Path

import numpy as np
import pytest

from swim.source import initial_vapor_from_climatology

FIXTURE_DIR = (
    Path(__file__).parents[1]
    / "fixtures"
    / "matlab"
    / "port_baseline_v1"
    / "allan_hills"
)


@pytest.mark.parity
def test_initial_vapor_matches_frozen_matlab_trajectory() -> None:
    metadata = json.loads((FIXTURE_DIR / "metadata.json").read_text())
    expected = metadata["result"]["trajectory_initial_conditions"]

    actual = initial_vapor_from_climatology(10.0)

    np.testing.assert_allclose(actual.delta_d_permil, expected["dD_v0"], atol=1e-11)
    np.testing.assert_allclose(
        actual.delta_18o_permil, expected["d18O_v0"], atol=1e-11
    )
    np.testing.assert_allclose(
        actual.oxygen_17_excess_log, expected["d17Oxs_v0"], atol=1e-14
    )
    np.testing.assert_allclose(
        actual.normalized_relative_humidity, expected["RHn0"], atol=1e-14
    )
    np.testing.assert_allclose(actual.relative_humidity, expected["RH0"], atol=1e-14)
    np.testing.assert_allclose(
        actual.sea_surface_temperature_c, expected["SST0_C"], atol=1e-13
    )

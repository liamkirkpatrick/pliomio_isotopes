from pathlib import Path

import numpy as np
import pytest

from swim.model import StateSpace, load_matlab_state_space
from swim.reconstruction import reconstruct_temperatures

FIXTURE_DIR = (
    Path(__file__).parents[1]
    / "fixtures"
    / "matlab"
    / "port_baseline_v1"
    / "allan_hills"
)


def _fixture_state_space() -> StateSpace:
    return load_matlab_state_space(FIXTURE_DIR / "state_space.mat")


@pytest.mark.parity
def test_allan_hills_reconstruction_matches_matlab_natural_neighbors() -> None:
    observations = np.genfromtxt(
        FIXTURE_DIR / "allan_hills_reconstruction.csv",
        delimiter=",",
        names=True,
        encoding=None,
    )
    result = reconstruct_temperatures(
        observations["dD_smow"], observations["d18O_smow"], _fixture_state_space()
    )

    np.testing.assert_allclose(
        result.condensation_temperature_c,
        observations["Tcond_reconstructed_C"],
        atol=2.0e-11,
        equal_nan=True,
    )
    np.testing.assert_allclose(
        result.source_temperature_c,
        observations["Tsource_reconstructed_C"],
        atol=2.0e-11,
        equal_nan=True,
    )
    np.testing.assert_allclose(
        result.surface_temperature_c,
        observations["Tsurface_reconstructed_C"],
        atol=3.0e-11,
        equal_nan=True,
    )
    np.testing.assert_allclose(
        result.saturated_mixing_ratio,
        observations["r_s_reconstructed"],
        atol=1.0e-15,
        equal_nan=True,
    )


@pytest.mark.parametrize("method", [2, 3])
def test_alternate_legacy_reconstruction_methods_execute(method: int) -> None:
    observations = np.genfromtxt(
        FIXTURE_DIR / "allan_hills_reconstruction.csv",
        delimiter=",",
        names=True,
        encoding=None,
    )

    result = reconstruct_temperatures(
        observations["dD_smow"][:20],
        observations["d18O_smow"][:20],
        _fixture_state_space(),
        method=method,  # type: ignore[arg-type]
    )

    assert result.condensation_temperature_c.shape == (20,)

from pathlib import Path

import pytest
from scipy.io import savemat

from swim.model import load_matlab_state_space

FIXTURE_PATH = (
    Path(__file__).parents[1]
    / "fixtures"
    / "matlab"
    / "port_baseline_v1"
    / "allan_hills"
    / "state_space.mat"
)


def test_load_matlab_state_space() -> None:
    state_space = load_matlab_state_space(FIXTURE_PATH)

    assert state_space.source_temperature_c.shape == (29,)
    assert state_space.condensation_temperature_c.shape == (71,)
    assert state_space.delta_18o.shape == (29, 71)


def test_load_matlab_state_space_accepts_string_path(tmp_path: Path) -> None:
    path = tmp_path / "not-a-state-space.mat"
    path.write_bytes(FIXTURE_PATH.read_bytes())

    assert load_matlab_state_space(str(path)).delta_d.shape == (29, 71)


def test_load_matlab_state_space_rejects_missing_arrays(tmp_path: Path) -> None:
    path = tmp_path / "incomplete.mat"
    savemat(path, {"T_source": [0.0, 1.0]})

    with pytest.raises(ValueError, match="missing"):
        load_matlab_state_space(path)


def test_load_matlab_state_space_reports_missing_file(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        load_matlab_state_space(tmp_path / "missing.mat")

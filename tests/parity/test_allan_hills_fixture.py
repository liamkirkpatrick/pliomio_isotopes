import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.io import loadmat

FIXTURE_DIR = (
    Path(__file__).parents[1]
    / "fixtures"
    / "matlab"
    / "port_baseline_v1"
    / "allan_hills"
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_fixture_metadata_and_checksums() -> None:
    metadata = json.loads((FIXTURE_DIR / "metadata.json").read_text())

    assert metadata["baseline_id"] == "matlab-port-baseline-v1"
    assert metadata["baseline_status"] == "frozen_for_python_parity"
    assert metadata["project"]["git_dirty"] is False

    for output in metadata["outputs"].values():
        output_path = FIXTURE_DIR / output["file"]
        assert output_path.is_file()
        assert _sha256(output_path) == output["sha256"]


def test_state_space_fixture_is_python_readable() -> None:
    state = loadmat(FIXTURE_DIR / "state_space.mat")

    assert state["T_source"].shape == (1, 29)
    assert state["T_site"].shape == (1, 71)
    assert state["d18O_site"].shape == (29, 71)
    assert state["dlnU_site"].shape == (29, 71)


def test_trajectory_fixture_is_python_readable() -> None:
    trajectory = loadmat(FIXTURE_DIR / "trajectory_Tsource_10_Tcond_m30.mat")

    np.testing.assert_allclose(trajectory["trajectory_T"][0, [0, -1]], [10, -30])
    assert trajectory["trajectory_T"].shape == (1, 401)
    assert trajectory["trajectory_P"].shape == (1, 401)
    assert trajectory["trajectory_f"].shape == (1, 401)


def test_allan_hills_reconstruction_coverage() -> None:
    with (FIXTURE_DIR / "allan_hills_reconstruction.csv").open(newline="") as file:
        rows = list(csv.DictReader(file))

    reconstructed = [
        row
        for row in rows
        if np.isfinite(float(row["Tcond_reconstructed_C"]))
        and np.isfinite(float(row["Tsource_reconstructed_C"]))
    ]

    assert len(rows) == 2326
    assert len(reconstructed) == 2248

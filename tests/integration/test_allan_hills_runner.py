import importlib.util
from pathlib import Path

import numpy as np

SCRIPT_PATH = Path(__file__).parents[2] / "scripts" / "python" / "run_allan_hills.py"


def _load_runner():  # type: ignore[no-untyped-def]
    spec = importlib.util.spec_from_file_location("run_allan_hills", SCRIPT_PATH)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_allan_hills_runner_writes_reconstruction(tmp_path: Path) -> None:
    runner = _load_runner()
    fixture_dir = (
        Path(__file__).parents[1]
        / "fixtures"
        / "matlab"
        / "port_baseline_v1"
        / "allan_hills"
    )
    output = tmp_path / "reconstruction.csv"

    count = runner.run(
        Path(__file__).parents[1] / "test_data.csv",
        fixture_dir / "state_space.mat",
        output,
    )
    values = np.genfromtxt(output, delimiter=",", names=True, encoding=None)

    assert count == 2248
    assert values.shape == (2326,)
    assert "Tsurface_reconstructed_C" in values.dtype.names

#!/usr/bin/env python3
"""Reconstruct Allan Hills temperatures with corrected SWIM evaporation."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from swim.model import forward_state_space, load_matlab_state_space
from swim.reconstruction import reconstruct_temperatures

PROJECT_ROOT = Path(__file__).parents[2]
DEFAULT_INPUT = PROJECT_ROOT / "tests" / "test_data.csv"
DEFAULT_OUTPUT = PROJECT_ROOT / "tests" / "generated" / "python_allan_hills.csv"


def _arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument(
        "--state-space",
        type=Path,
        help=(
            "Load an existing MATLAB-format state space. If omitted, generate "
            "a new state space with the corrected 2022 evaporation default."
        ),
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def run(
    input_path: Path, state_space_path: Path | None, output_path: Path
) -> int:
    """Run the reconstruction and return the number of finite results written."""
    with input_path.open(newline="", encoding="utf-8-sig") as input_file:
        reader = csv.DictReader(input_file)
        rows = list(reader)
        fieldnames = reader.fieldnames
    if fieldnames is None:
        raise ValueError(f"Input CSV has no header: {input_path}")
    required = {"d18O_smow", "dD_smow"}
    missing = sorted(required.difference(fieldnames))
    if missing:
        raise ValueError(f"Input CSV is missing: {', '.join(missing)}")

    def numeric(column: str) -> np.ndarray:
        return np.asarray(
            [float(row[column]) if row[column].strip() else np.nan for row in rows]
        )

    if state_space_path is None:
        state_space = forward_state_space(
            np.arange(0.0, 29.0), np.arange(-70.0, 1.0)
        )
    else:
        state_space = load_matlab_state_space(state_space_path)
    result = reconstruct_temperatures(
        numeric("dD_smow"), numeric("d18O_smow"), state_space
    )
    output_columns = {
        "Tcond_reconstructed_C": result.condensation_temperature_c,
        "Tsource_reconstructed_C": result.source_temperature_c,
        "Tsurface_reconstructed_C": result.surface_temperature_c,
        "r_s_reconstructed": result.saturated_mixing_ratio,
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as output_file:
        writer = csv.DictWriter(
            output_file, fieldnames=list(fieldnames) + list(output_columns)
        )
        writer.writeheader()
        for index, row in enumerate(rows):
            writer.writerow(
                {
                    **row,
                    **{
                        name: float(values[index])
                        for name, values in output_columns.items()
                    },
                }
            )
    return int(
        np.count_nonzero(
            np.isfinite(result.condensation_temperature_c)
            & np.isfinite(result.source_temperature_c)
        )
    )


def main() -> None:
    args = _arguments()
    reconstructed = run(args.input, args.state_space, args.output)
    print(f"Wrote reconstruction to {args.output}")
    print(f"Finite temperature reconstructions: {reconstructed}")


if __name__ == "__main__":
    main()

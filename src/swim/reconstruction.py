"""Nonlinear inversion of the SWIM isotope state space."""

from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.interpolate import griddata
from scipy.io import loadmat

from swim.interpolation import natural_neighbor_interpolate
from swim.isotopes import log_delta, logarithmic_deuterium_excess
from swim.model import StateSpace
from swim.source import DEFAULT_LEGACY_DATA_DIR


@dataclass(frozen=True)
class ReconstructionResult:
    """Reconstructed source, condensation, and surface temperatures."""

    condensation_temperature_c: NDArray[np.float64]
    source_temperature_c: NDArray[np.float64]
    surface_temperature_c: NDArray[np.float64]
    saturated_mixing_ratio: NDArray[np.float64]


@dataclass(frozen=True)
class SeawaterCorrectionResult:
    """Log-space seawater-corrected isotope observations."""

    delta_18o_log: NDArray[np.float64]
    delta_d_log: NDArray[np.float64]
    logarithmic_deuterium_excess: NDArray[np.float64]
    delta_18o_permil: NDArray[np.float64]
    delta_d_permil: NDArray[np.float64]


def surface_temperature_from_condensation(
    condensation_temperature_c: ArrayLike,
) -> NDArray[np.float64]:
    """Invert the active legacy relation ``Tc = 0.69*Ts - 8.2``."""
    condensation = np.asarray(condensation_temperature_c, dtype=np.float64)
    return (condensation + 8.2) / 0.69


def reconstruct_temperatures(
    delta_d_permil: ArrayLike,
    delta_18o_permil: ArrayLike,
    state_space: StateSpace,
    *,
    method: Literal[1, 2, 3] = 1,
) -> ReconstructionResult:
    """Invert a legacy isotope coordinate representation.

    Method 1 is the preferred natural-neighbor log-18O/d-ln reconstruction.
    Methods 2 and 3 preserve the active linear interpolation statements and
    their inconsistent query coordinates from the legacy helper.
    """
    delta_d = np.asarray(delta_d_permil, dtype=np.float64)
    delta_18o = np.asarray(delta_18o_permil, dtype=np.float64)
    if delta_d.shape != delta_18o.shape:
        raise ValueError("delta-D and delta-18O observations must have equal shapes")
    original_shape = delta_d.shape
    delta_d_flat = delta_d.ravel()
    delta_18o_flat = delta_18o.ravel()
    valid_observations = np.isfinite(delta_d_flat) & np.isfinite(delta_18o_flat)

    state_valid = np.isfinite(state_space.delta_18o_log)
    source_grid = np.repeat(
        state_space.source_temperature_c[:, None],
        state_space.condensation_temperature_c.size,
        axis=1,
    )
    condensation_grid = np.repeat(
        state_space.condensation_temperature_c[None, :],
        state_space.source_temperature_c.size,
        axis=0,
    )
    observed_18o_log = np.asarray(log_delta(delta_18o_flat[valid_observations]))
    observed_d_log = np.asarray(log_delta(delta_d_flat[valid_observations]))
    observed_dln = np.asarray(
        logarithmic_deuterium_excess(
            delta_d_flat[valid_observations], delta_18o_flat[valid_observations]
        )
    )
    observed_dxs = (
        delta_d_flat[valid_observations]
        - 8.0 * delta_18o_flat[valid_observations]
    )

    if method == 1:
        points = np.column_stack(
            (
                state_space.delta_18o_log[state_valid],
                state_space.logarithmic_deuterium_excess[state_valid],
            )
        )
        query = np.column_stack((observed_18o_log, observed_dln))
    elif method == 2:
        points = np.column_stack(
            (
                state_space.delta_18o_log[state_valid],
                state_space.delta_d_log[state_valid],
            )
        )
        query = np.column_stack((observed_18o_log, observed_d_log))
    elif method == 3:
        points = np.column_stack(
            (
                state_space.delta_18o_log[state_valid],
                state_space.deuterium_excess[state_valid],
            )
        )
        query = np.column_stack(
            (delta_18o_flat[valid_observations], observed_dxs)
        )
    else:
        raise ValueError(f"Unknown reconstruction method: {method}")

    natural_values = np.column_stack(
        (condensation_grid[state_valid], source_grid[state_valid])
    )
    if method == 1:
        temperature_result = natural_neighbor_interpolate(points, natural_values, query)
    else:
        temperature_result = griddata(points, natural_values, query, method="linear")

    # These inconsistent r_s query coordinates are active in methods 2 and 3.
    mixing_query = query if method == 1 else np.column_stack(
        (observed_18o_log, observed_dln)
    )
    mixing_ratio = griddata(
        points,
        state_space.saturated_mixing_ratio[state_valid],
        mixing_query,
        method="linear",
    )

    condensation = np.full(delta_d_flat.shape, np.nan)
    source = np.full(delta_d_flat.shape, np.nan)
    reconstructed_mixing_ratio = np.full(delta_d_flat.shape, np.nan)
    condensation[valid_observations] = temperature_result[:, 0]
    source[valid_observations] = temperature_result[:, 1]
    reconstructed_mixing_ratio[valid_observations] = mixing_ratio
    surface = surface_temperature_from_condensation(condensation)
    return ReconstructionResult(
        condensation_temperature_c=condensation.reshape(original_shape),
        source_temperature_c=source.reshape(original_shape),
        surface_temperature_c=surface.reshape(original_shape),
        saturated_mixing_ratio=reconstructed_mixing_ratio.reshape(original_shape),
    )


def seawater_correct_isotopes(
    delta_18o_permil: ArrayLike,
    delta_d_permil: ArrayLike,
    age_years: ArrayLike,
    initial_seawater_delta_18o_permil: float,
    *,
    data_dir: Path = DEFAULT_LEGACY_DATA_DIR,
) -> SeawaterCorrectionResult:
    """Port the log-space ice-volume correction in ``seawater_cor_ln.m``."""
    delta_18o = np.asarray(delta_18o_permil, dtype=np.float64)
    delta_d = np.asarray(delta_d_permil, dtype=np.float64)
    age = np.asarray(age_years, dtype=np.float64)
    if delta_18o.shape != delta_d.shape or delta_18o.shape != age.shape:
        raise ValueError("isotopes and ages must have equal shapes")

    data = loadmat(data_dir / "Bintaja.mat", simplify_cells=True)["Bintaja"]
    model_age = np.asarray(data["time"], dtype=np.float64) * 1000.0
    model_ice = np.asarray(data["iso_ice"], dtype=np.float64)
    seawater_delta_18o = np.interp(age, model_age, model_ice)
    outside = (age < model_age[0]) | (age > model_age[-1])
    seawater_delta_18o = np.where(outside, np.nan, seawater_delta_18o)
    offset = model_ice[0] - initial_seawater_delta_18o_permil
    seawater_delta_18o = seawater_delta_18o - offset

    seawater_18o_log = np.asarray(log_delta(seawater_delta_18o))
    initial_seawater_18o_log = np.asarray(
        log_delta(initial_seawater_delta_18o_permil)
    )
    change_18o_log = seawater_18o_log - initial_seawater_18o_log
    seawater_d_log = -0.0285 * seawater_18o_log**2 + 8.47 * seawater_18o_log
    initial_seawater_d_log = (
        -0.0285 * initial_seawater_18o_log**2
        + 8.47 * initial_seawater_18o_log
    )
    change_d_log = seawater_d_log - initial_seawater_d_log

    corrected_18o_log = np.asarray(log_delta(delta_18o)) - change_18o_log
    corrected_d_log = np.asarray(log_delta(delta_d)) - change_d_log
    corrected_dln = corrected_d_log - (
        -0.0285 * corrected_18o_log**2 + 8.47 * corrected_18o_log
    )
    corrected_18o = (np.exp(corrected_18o_log / 1000.0) - 1.0) * 1000.0
    corrected_d = (np.exp(corrected_d_log / 1000.0) - 1.0) * 1000.0
    return SeawaterCorrectionResult(
        delta_18o_log=corrected_18o_log,
        delta_d_log=corrected_d_log,
        logarithmic_deuterium_excess=corrected_dln,
        delta_18o_permil=corrected_18o,
        delta_d_permil=corrected_d,
    )


__all__ = [
    "ReconstructionResult",
    "SeawaterCorrectionResult",
    "reconstruct_temperatures",
    "seawater_correct_isotopes",
    "surface_temperature_from_condensation",
]

"""Composed forward-model workflows for the frozen legacy SWIM baseline."""

from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import numpy as np
from numpy.typing import NDArray
from scipy.io import loadmat

from swim.cloud_phase import cloud_phase_fractions
from swim.distillation import DistillationResult, rayleigh_distillation
from swim.fractionation import (
    O17_ICE_EXPONENT_DISTILLATION,
    EffectiveFractionationFactors,
    EquilibriumFractionationFactors,
    equilibrium_fractionation_factors,
    kinetic_condensation_factors,
    mixed_phase_effective_fractionation,
    transport_diffusivity_ratios,
)
from swim.source import (
    DEFAULT_LEGACY_DATA_DIR,
    Hemisphere,
    InitialVapor,
    initial_vapor_from_climatology,
)
from swim.thermodynamics import (
    PseudoAdiabatResult,
    mixed_phase_supersaturation,
    pseudo_adiabat,
)


@dataclass(frozen=True)
class ForwardTrajectory:
    """Complete intermediate state for one source-to-condensation pathway."""

    temperature_c: NDArray[np.float64]
    fraction_ice: NDArray[np.float64]
    fraction_liquid: NDArray[np.float64]
    mixed_phase_supersaturation: NDArray[np.float64]
    source_vapor: InitialVapor
    thermodynamics: PseudoAdiabatResult
    equilibrium_fractionation: EquilibriumFractionationFactors
    effective_fractionation: EffectiveFractionationFactors
    distillation: DistillationResult


@dataclass(frozen=True)
class StateSpace:
    """Endpoint products over source and condensation temperature grids."""

    source_temperature_c: NDArray[np.float64]
    condensation_temperature_c: NDArray[np.float64]
    delta_18o: NDArray[np.float64]
    delta_d: NDArray[np.float64]
    delta_18o_log: NDArray[np.float64]
    delta_d_log: NDArray[np.float64]
    deuterium_excess: NDArray[np.float64]
    oxygen_17_excess_per_meg: NDArray[np.float64]
    logarithmic_deuterium_excess: NDArray[np.float64]
    saturated_mixing_ratio: NDArray[np.float64]
    pressure_kpa: NDArray[np.float64]


def matlab_temperature_grid(
    start_temperature_c: float,
    end_temperature_c: float,
    step_c: float,
) -> NDArray[np.float64]:
    """Construct the descending, inclusive grid used by MATLAB ``start:-dT:end``."""
    if step_c <= 0.0:
        raise ValueError("step_c must be positive")
    if end_temperature_c > start_temperature_c:
        raise ValueError("end temperature cannot exceed start temperature")
    intervals = int(np.floor((start_temperature_c - end_temperature_c) / step_c))
    temperature = start_temperature_c - np.arange(intervals + 1) * step_c
    tolerance = np.finfo(np.float64).eps * max(abs(end_temperature_c), 1.0) * 8.0
    if abs(temperature[-1] - end_temperature_c) <= tolerance:
        temperature[-1] = end_temperature_c
    return temperature


def _preserve_legacy_initial_fractionation(
    calculated: EffectiveFractionationFactors,
    equilibrium: EquilibriumFractionationFactors,
    fraction_ice: NDArray[np.float64],
    fraction_liquid: NDArray[np.float64],
) -> EffectiveFractionationFactors:
    alpha_d = np.asarray(calculated.deuterium).copy()
    alpha_18o = np.asarray(calculated.oxygen_18).copy()
    alpha_17o = np.asarray(calculated.oxygen_17).copy()
    alpha_d[0] = (
        fraction_ice[0] * np.asarray(equilibrium.deuterium_ice)[0]
        + fraction_liquid[0] * np.asarray(equilibrium.deuterium_liquid)[0]
    )
    alpha_18o[0] = (
        fraction_ice[0] * np.asarray(equilibrium.oxygen_18_ice)[0]
        + fraction_liquid[0] * np.asarray(equilibrium.oxygen_18_liquid)[0]
    )
    alpha_17o[0] = (
        fraction_ice[0] * np.asarray(equilibrium.oxygen_17_ice)[0]
        + fraction_liquid[0] * np.asarray(equilibrium.oxygen_17_liquid)[0]
    )
    return EffectiveFractionationFactors(alpha_d, alpha_18o, alpha_17o)


def forward_trajectory(
    source_temperature_c: float,
    condensation_temperature_c: float,
    *,
    step_c: float = 0.1,
    initial_pressure_kpa: float = 101.325,
    supersaturation_a: float = 1.0,
    supersaturation_b: float = 0.00525,
    supersaturation_c: float = 0.0,
    closure: Literal["local", "global"] = "local",
    hemisphere: Hemisphere = "south",
    reanalysis: Literal["ncep", "era"] = "ncep",
    data_dir: Path = DEFAULT_LEGACY_DATA_DIR,
) -> ForwardTrajectory:
    """Run one complete frozen-baseline SWIM forward trajectory."""
    temperature = matlab_temperature_grid(
        source_temperature_c, condensation_temperature_c, step_c
    )
    source_vapor = initial_vapor_from_climatology(
        source_temperature_c,
        closure=closure,
        hemisphere=hemisphere,
        reanalysis=reanalysis,
        data_dir=data_dir,
    )
    fraction_ice_raw, fraction_liquid_raw = cloud_phase_fractions(
        temperature, method="adj"
    )
    fraction_ice = np.asarray(fraction_ice_raw)
    fraction_liquid = np.asarray(fraction_liquid_raw)
    mixed_saturation = mixed_phase_supersaturation(
        temperature,
        initial_pressure_kpa,
        fraction_ice,
        fraction_liquid,
    )
    thermodynamics = pseudo_adiabat(
        temperature,
        initial_pressure_kpa,
        fraction_ice,
        fraction_liquid,
        mixed_saturation,
        supersaturation_a,
        supersaturation_b,
        supersaturation_c,
    )
    equilibrium = equilibrium_fractionation_factors(
        temperature,
        oxygen_17_ice_exponent=O17_ICE_EXPONENT_DISTILLATION,
    )
    kinetic = kinetic_condensation_factors(
        thermodynamics.saturation_over_ice,
        equilibrium,
        transport_diffusivity_ratios(temperature),
    )
    calculated_effective = mixed_phase_effective_fractionation(
        equilibrium, kinetic, fraction_ice, fraction_liquid
    )
    effective = _preserve_legacy_initial_fractionation(
        calculated_effective, equilibrium, fraction_ice, fraction_liquid
    )
    distillation = rayleigh_distillation(
        thermodynamics.fraction_vapor_remaining,
        effective,
        initial_delta_d_permil=float(source_vapor.delta_d_permil),
        initial_delta_18o_permil=float(source_vapor.delta_18o_permil),
        initial_oxygen_17_excess_log=float(source_vapor.oxygen_17_excess_log),
    )
    return ForwardTrajectory(
        temperature_c=temperature,
        fraction_ice=fraction_ice,
        fraction_liquid=fraction_liquid,
        mixed_phase_supersaturation=mixed_saturation,
        source_vapor=source_vapor,
        thermodynamics=thermodynamics,
        equilibrium_fractionation=equilibrium,
        effective_fractionation=effective,
        distillation=distillation,
    )


def forward_state_space(
    source_temperature_c: NDArray[np.float64],
    condensation_temperature_c: NDArray[np.float64],
    *,
    trajectory_step_c: float = 0.1,
    supersaturation_a: float = 1.0,
    supersaturation_b: float = 0.00525,
    supersaturation_c: float = 0.0,
    initial_pressure_kpa: float = 101.325,
    closure: Literal["local", "global"] = "local",
    hemisphere: Hemisphere = "south",
    reanalysis: Literal["ncep", "era"] = "ncep",
    data_dir: Path = DEFAULT_LEGACY_DATA_DIR,
) -> StateSpace:
    """Build the endpoint state space from ``simple_water_isotope_model_2020``."""
    source_grid = np.asarray(source_temperature_c, dtype=np.float64)
    condensation_grid = np.asarray(condensation_temperature_c, dtype=np.float64)
    if source_grid.ndim != 1 or condensation_grid.ndim != 1:
        raise ValueError("state-space temperature grids must be one-dimensional")
    shape = (source_grid.size, condensation_grid.size)
    outputs = [np.full(shape, np.nan, dtype=np.float64) for _ in range(9)]
    (
        delta_18o,
        delta_d,
        delta_18o_log,
        delta_d_log,
        deuterium_excess,
        oxygen_17_excess,
        logarithmic_deuterium_excess,
        mixing_ratio,
        pressure,
    ) = outputs

    for source_index, source_temperature in enumerate(source_grid):
        for site_index, condensation_temperature in enumerate(condensation_grid):
            if source_temperature < condensation_temperature:
                continue
            trajectory = forward_trajectory(
                float(source_temperature),
                float(condensation_temperature),
                step_c=trajectory_step_c,
                initial_pressure_kpa=initial_pressure_kpa,
                supersaturation_a=supersaturation_a,
                supersaturation_b=supersaturation_b,
                supersaturation_c=supersaturation_c,
                closure=closure,
                hemisphere=hemisphere,
                reanalysis=reanalysis,
                data_dir=data_dir,
            )
            isotope = trajectory.distillation
            thermodynamics = trajectory.thermodynamics
            delta_18o[source_index, site_index] = isotope.delta_18o_precipitation[-1]
            delta_d[source_index, site_index] = isotope.delta_d_precipitation[-1]
            delta_18o_log[source_index, site_index] = isotope.delta_18o_log[-1]
            delta_d_log[source_index, site_index] = isotope.delta_d_log[-1]
            deuterium_excess[source_index, site_index] = isotope.deuterium_excess[-1]
            oxygen_17_excess[source_index, site_index] = (
                isotope.oxygen_17_excess_per_meg[-1]
            )
            logarithmic_deuterium_excess[source_index, site_index] = (
                isotope.logarithmic_deuterium_excess[-1]
            )
            mixing_ratio[source_index, site_index] = (
                thermodynamics.saturated_mixing_ratio[-1]
            )
            pressure[source_index, site_index] = thermodynamics.pressure_kpa[-1]

    return StateSpace(
        source_temperature_c=source_grid,
        condensation_temperature_c=condensation_grid,
        delta_18o=delta_18o,
        delta_d=delta_d,
        delta_18o_log=delta_18o_log,
        delta_d_log=delta_d_log,
        deuterium_excess=deuterium_excess,
        oxygen_17_excess_per_meg=oxygen_17_excess,
        logarithmic_deuterium_excess=logarithmic_deuterium_excess,
        saturated_mixing_ratio=mixing_ratio,
        pressure_kpa=pressure,
    )


def load_matlab_state_space(path: Path | str) -> StateSpace:
    """Load a state space exported by the MATLAB validation runner."""
    values = loadmat(str(path), simplify_cells=True)
    names = {
        "source_temperature_c": "T_source",
        "condensation_temperature_c": "T_site",
        "delta_18o": "d18O_site",
        "delta_d": "dD_site",
        "delta_18o_log": "d18Oln_site",
        "delta_d_log": "dDln_site",
        "deuterium_excess": "dxs_site",
        "oxygen_17_excess_per_meg": "d17O_xs_site",
        "logarithmic_deuterium_excess": "dlnU_site",
        "saturated_mixing_ratio": "r_s_site",
        "pressure_kpa": "P_site",
    }
    missing = sorted(
        matlab_name for matlab_name in names.values() if matlab_name not in values
    )
    if missing:
        raise ValueError(f"MATLAB state-space file is missing: {', '.join(missing)}")
    arrays = {
        python_name: np.asarray(values[matlab_name], dtype=np.float64)
        for python_name, matlab_name in names.items()
    }
    if arrays["source_temperature_c"].ndim != 1:
        arrays["source_temperature_c"] = arrays["source_temperature_c"].ravel()
    if arrays["condensation_temperature_c"].ndim != 1:
        arrays["condensation_temperature_c"] = arrays[
            "condensation_temperature_c"
        ].ravel()
    expected_shape = (
        arrays["source_temperature_c"].size,
        arrays["condensation_temperature_c"].size,
    )
    for name, array in arrays.items():
        if name.endswith("temperature_c"):
            continue
        if array.shape != expected_shape:
            raise ValueError(
                f"MATLAB state-space array {names[name]!r} has shape {array.shape}; "
                f"expected {expected_shape}"
            )
    return StateSpace(**arrays)


__all__ = [
    "ForwardTrajectory",
    "StateSpace",
    "forward_state_space",
    "forward_trajectory",
    "load_matlab_state_space",
    "matlab_temperature_grid",
]

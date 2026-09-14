"""Source-region climatology and ocean evaporation for legacy SWIM."""

from dataclasses import dataclass
from functools import cache
from pathlib import Path
from typing import Literal, TypeAlias

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.interpolate import PPoly
from scipy.io import loadmat

from swim.fractionation import (
    O17_ICE_EXPONENT_EVAPORATION,
    equilibrium_fractionation_factors,
    transport_diffusivity_ratios,
)
from swim.isotopes import R18O_VSMOW, RD_VSMOW

FloatResult: TypeAlias = np.float64 | NDArray[np.float64]
Hemisphere: TypeAlias = Literal["south", "north", "all"]

DEFAULT_LEGACY_DATA_DIR = Path(__file__).parents[2] / "legacy_matlab" / "data"


@dataclass(frozen=True)
class SourceConditions:
    """Climatological source environment for a surface-air temperature."""

    relative_humidity: FloatResult
    relative_humidity_uncertainty_percent: FloatResult
    sea_surface_temperature_c: FloatResult
    sea_surface_temperature_uncertainty_c: FloatResult
    normalized_relative_humidity: FloatResult
    normalized_relative_humidity_uncertainty: FloatResult


@dataclass(frozen=True)
class InitialVapor:
    """Initial vapor isotope composition and its source environment."""

    delta_d_permil: FloatResult
    delta_18o_permil: FloatResult
    oxygen_17_excess_log: FloatResult
    normalized_relative_humidity: FloatResult
    relative_humidity: FloatResult
    sea_surface_temperature_c: FloatResult


def _spline_filename(reanalysis: str, hemisphere: Hemisphere) -> str:
    suffix = {"south": "SH", "north": "NH", "all": "all"}[hemisphere]
    if reanalysis == "ncep":
        return f"ncep_data/ncep_spline_model_{suffix}.mat"
    if reanalysis == "era":
        return f"era_data/era_spline_model_{suffix}.mat"
    raise ValueError(f"Unknown reanalysis: {reanalysis!r}")


@cache
def _load_spline_models(path: str) -> tuple[PPoly, PPoly, NDArray, NDArray, NDArray]:
    data = loadmat(path, simplify_cells=True)

    def convert(name: str) -> PPoly:
        spline = data[name]
        return PPoly(
            np.asarray(spline["coefs"], dtype=np.float64).T,
            np.asarray(spline["breaks"], dtype=np.float64),
            extrapolate=True,
        )

    return (
        convert("sp_sst"),
        convert("sp_rh"),
        np.asarray(data["T_model_spline"], dtype=np.float64),
        np.asarray(data["delta_sst_smooth"], dtype=np.float64),
        np.asarray(data["delta_rh_smooth"], dtype=np.float64),
    )


def _matlab_interp1_no_extrapolation(
    x: NDArray[np.float64], y: NDArray[np.float64], values: NDArray[np.float64]
) -> NDArray[np.float64]:
    result = np.interp(values, x, y)
    return np.where((values < x[0]) | (values > x[-1]), np.nan, result)


def _legacy_liquid_vapor_pressure_from_bad_kelvin(
    temperature: NDArray[np.float64],
) -> NDArray[np.complex128]:
    """Evaluate the 2020 routine's liquid formula after its C-to-K sign bug."""
    temperature_complex = temperature.astype(np.complex128)
    log_temperature = np.log(temperature_complex)
    return 1.0e-3 * np.exp(
        54.842763
        - 6763.22 / temperature_complex
        - 4.21 * log_temperature
        + 0.000367 * temperature_complex
        + np.tanh(0.0415 * (temperature_complex - 218.8))
        * (
            53.878
            - 1331.22 / temperature_complex
            - 9.44523 * log_temperature
            + 0.014025 * temperature_complex
        )
    )


def climatological_source_conditions(
    source_air_temperature_c: ArrayLike,
    *,
    hemisphere: Hemisphere = "south",
    reanalysis: Literal["ncep", "era"] = "ncep",
    data_dir: Path = DEFAULT_LEGACY_DATA_DIR,
) -> SourceConditions:
    """Evaluate the selected legacy smoothing-spline climatology.

    This deliberately preserves the erroneous Celsius-to-Kelvin subtraction
    in ``T_RH_RHn_2020.m`` because the frozen MATLAB baseline used it.
    Relative humidity is returned as a fraction, not MATLAB's intermediate
    percent value.
    """
    temperature = np.asarray(source_air_temperature_c, dtype=np.float64)
    spline_path = data_dir / _spline_filename(reanalysis, hemisphere)
    sst_spline, rh_spline, model_temperature, delta_sst, delta_rh = (
        _load_spline_models(str(spline_path))
    )
    sea_surface_temperature = np.asarray(sst_spline(temperature))
    relative_humidity_percent = np.asarray(rh_spline(temperature))
    sea_surface_uncertainty = _matlab_interp1_no_extrapolation(
        model_temperature, delta_sst, temperature
    )
    relative_humidity_uncertainty = _matlab_interp1_no_extrapolation(
        model_temperature, delta_rh, temperature
    )

    skin_vapor_pressure = _legacy_liquid_vapor_pressure_from_bad_kelvin(
        temperature - 273.15
    )
    sst_vapor_pressure = _legacy_liquid_vapor_pressure_from_bad_kelvin(
        sea_surface_temperature - 273.15
    )
    normalized_percent = np.real(
        relative_humidity_percent * skin_vapor_pressure / sst_vapor_pressure
    )
    normalized_uncertainty = (
        np.real(
            (relative_humidity_percent + relative_humidity_uncertainty)
            * skin_vapor_pressure
            / sst_vapor_pressure
        )
        - normalized_percent
    ) / 100.0

    return SourceConditions(
        relative_humidity=relative_humidity_percent / 100.0,
        relative_humidity_uncertainty_percent=relative_humidity_uncertainty,
        sea_surface_temperature_c=sea_surface_temperature,
        sea_surface_temperature_uncertainty_c=sea_surface_uncertainty,
        normalized_relative_humidity=normalized_percent / 100.0,
        normalized_relative_humidity_uncertainty=normalized_uncertainty,
    )


@cache
def _seawater_deuterium_fit(data_path: str) -> NDArray[np.float64]:
    observations = np.loadtxt(data_path, comments="%")
    delta_18o = observations[:, 4]
    delta_d = observations[:, 5]
    valid = (delta_d != -99.9) & (delta_18o != -99.9)
    return np.polyfit(delta_18o[valid], delta_d[valid], 1)


def seawater_delta_d_from_delta_18o(
    delta_18o_permil: ArrayLike,
    *,
    data_dir: Path = DEFAULT_LEGACY_DATA_DIR,
) -> FloatResult:
    """Evaluate the legacy linear seawater D--18O observational fit."""
    coefficients = _seawater_deuterium_fit(
        str(data_dir / "NASA_GISS_sw_d18O_dD_short.txt")
    )
    return np.polyval(coefficients, np.asarray(delta_18o_permil, dtype=np.float64))


def initial_vapor_from_climatology(
    source_air_temperature_c: ArrayLike,
    *,
    closure: Literal["local", "global"] = "local",
    hemisphere: Hemisphere = "south",
    reanalysis: Literal["ncep", "era"] = "ncep",
    seawater_delta_18o_permil: float = -0.3,
    oxygen_18_diffusive_fractionation: float = 1.009,
    data_dir: Path = DEFAULT_LEGACY_DATA_DIR,
) -> InitialVapor:
    """Calculate initial vapor using the active ``evaporation_2021.m`` path."""
    temperature = np.asarray(source_air_temperature_c, dtype=np.float64)
    source = climatological_source_conditions(
        temperature,
        hemisphere=hemisphere,
        reanalysis=reanalysis,
        data_dir=data_dir,
    )
    equilibrium = equilibrium_fractionation_factors(
        temperature,
        oxygen_17_ice_exponent=O17_ICE_EXPONENT_EVAPORATION,
    )
    diffusivity = transport_diffusivity_ratios(temperature)
    hdo_ratio = np.asarray(diffusivity.hdo_over_h2o)
    oxygen_18_ratio = np.asarray(diffusivity.h218o_over_h216o)
    phi_diffusivity = (1.0 - hdo_ratio) / (1.0 - oxygen_18_ratio)

    # The argument remains for traceability to the legacy signature, but the
    # selected active branch overwrites its initial 1.009 value using N=0.302.
    _ = oxygen_18_diffusive_fractionation
    diffusivity_exponent = 0.302
    alpha_18_diffusive = (1.0 / oxygen_18_ratio) ** diffusivity_exponent
    alpha_d_diffusive = phi_diffusivity * (alpha_18_diffusive - 1.0) + 1.0
    alpha_18_equilibrium = np.asarray(equilibrium.oxygen_18_liquid)
    alpha_d_equilibrium = np.asarray(equilibrium.deuterium_liquid)
    normalized_humidity = np.asarray(source.normalized_relative_humidity)

    seawater_delta_d = seawater_delta_d_from_delta_18o(
        seawater_delta_18o_permil, data_dir=data_dir
    )
    seawater_18o_ratio = (1.0 + seawater_delta_18o_permil / 1000.0) * R18O_VSMOW
    seawater_d_ratio = (1.0 + seawater_delta_d / 1000.0) * RD_VSMOW

    if closure == "local":
        vapor_18o_ratio = seawater_18o_ratio / (
            alpha_18_equilibrium
            * (
                alpha_18_diffusive
                + normalized_humidity * (1.0 - alpha_18_diffusive)
            )
        )
        vapor_d_ratio = seawater_d_ratio / (
            alpha_d_equilibrium
            * (
                alpha_d_diffusive
                + normalized_humidity * (1.0 - alpha_d_diffusive)
            )
        )
    elif closure == "global":
        alpha_18_evaporation = 1.0045
        alpha_d_evaporation = 1.0267
        vapor_18o_ratio = seawater_18o_ratio * (
            1.0
            - alpha_18_equilibrium
            * alpha_18_diffusive
            * (1.0 - normalized_humidity)
            / alpha_18_evaporation
        ) / (alpha_18_equilibrium * normalized_humidity)
        vapor_d_ratio = seawater_d_ratio * (
            1.0
            - alpha_d_equilibrium
            * alpha_d_diffusive
            * (1.0 - normalized_humidity)
            / alpha_d_evaporation
        ) / (alpha_d_equilibrium * normalized_humidity)
    else:
        raise ValueError(f"Unknown closure: {closure!r}")

    delta_18o_vapor = (vapor_18o_ratio / R18O_VSMOW - 1.0) * 1000.0
    delta_d_vapor = (vapor_d_ratio / RD_VSMOW - 1.0) * 1000.0
    oxygen_17_excess_log = -np.log(
        alpha_18_equilibrium**0.529
        * (
            alpha_18_diffusive**0.518 * (1.0 - normalized_humidity)
            + normalized_humidity
        )
    ) + 0.528 * np.log(
        alpha_18_equilibrium
        * (
            alpha_18_diffusive * (1.0 - normalized_humidity)
            + normalized_humidity
        )
    )

    return InitialVapor(
        delta_d_permil=delta_d_vapor,
        delta_18o_permil=delta_18o_vapor,
        oxygen_17_excess_log=oxygen_17_excess_log,
        normalized_relative_humidity=source.normalized_relative_humidity,
        relative_humidity=source.relative_humidity,
        sea_surface_temperature_c=source.sea_surface_temperature_c,
    )


__all__ = [
    "InitialVapor",
    "SourceConditions",
    "climatological_source_conditions",
    "initial_vapor_from_climatology",
    "seawater_delta_d_from_delta_18o",
]

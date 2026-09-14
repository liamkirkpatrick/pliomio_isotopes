# SWIM scientific and implementation traceability

## Purpose

This is a living traceability matrix connecting:

```text
paper / scientific concept
        ↓
MATLAB implementation
        ↓
Python implementation
        ↓
parity / validation test
```

Paper:

> Markle, B. R., & Steig, E. J. (2022). Improving temperature reconstructions from ice-core water-isotope records. *Climate of the Past*, 18, 1321–1368. https://doi.org/10.5194/cp-18-1321-2022

MATLAB repository snapshot inspected:

- `bradley-markle/simple_water_isotope_model`
- commit `4db23b79aab8f2154254d9111f77e70a5ce8427d`

The frozen core baseline is now implemented. Remaining open items concern
publication provenance, optional historical branches, or later application and
uncertainty workflows rather than the validated core forward/inverse path.

---

## Status legend

- **VERIFIED-PAPER** — directly stated in the publication.
- **VERIFIED-CODE** — directly observed in active checked-in MATLAB statements.
- **LIKELY-ACTIVE** — inferred from the checked-in example/call chain, but publication provenance still needs confirmation.
- **OPEN** — discrepancy or provenance question requiring explicit resolution.
- **TBD** — Python implementation/test does not yet exist.

---

## 1. Core forward-model traceability

| Scientific concept | Paper source | MATLAB implementation | Code status | Python target | Parity / validation target | Notes |
|---|---|---|---|---|---|---|
| delta notation | Eq. (1) | `evaporation_2021.m`; `distillation_2020.m` | VERIFIED-CODE | `swim.isotopes.delta_to_ratio`; `ratio_to_delta` | unit tests + frozen trajectory ratios | MATLAB uses per-mil values and VSMOW ratios |
| VSMOW isotope ratios | background / standard | `evaporation_2021.m`; `distillation_2020.m` | VERIFIED-CODE | `swim.isotopes.R18O_VSMOW`; `RD_VSMOW`; `R17O_VSMOW` | exact-constant test | constants preserved exactly during parity |
| linear deuterium excess \(d_{xs}\) | Sect. 1.1 | `distillation_2020.m`; wrapper | VERIFIED-CODE | `swim.isotopes.linear_deuterium_excess` | unit test + frozen trajectory | \(dD-8d18O\) |
| logarithmic isotope transform | Sect. 1.2, Eq. (4) | `distillation_2020.m`; `Tsite_Tsource_reconstruction_quick.m` | VERIFIED-CODE | `swim.isotopes.log_delta` | unit test + frozen trajectory | MATLAB uses \(1000\ln(1+\delta/1000)\) |
| logarithmic deuterium excess \(d_{ln}\) | Eq. (4) | `simple_water_isotope_model_2020.m`; `Tsite_Tsource_reconstruction_quick.m` | VERIFIED-CODE | `swim.isotopes.logarithmic_deuterium_excess` | unit test + frozen trajectory | paper and code use different numerical scaling conventions but equivalent polynomial coefficients after scaling |
| \(^{17}O\) excess diagnostic | Sect. 1.2 | `distillation_2020.m`; wrapper | VERIFIED-CODE | `swim.isotopes.oxygen_17_excess` | unit test + frozen trajectory | returns per meg and preserves the legacy 0.528 exponent |
| source \(T_0\) → SST/RH climatology | Appendix A1.1 | `T_RH_RHn_2020.m`, `_2022.m`; data spline `.mat` files | VERIFIED-CODE | `swim.source.climatological_source_conditions` | frozen and corrected T0=10 cases + vector unit test | evaluates checked-in MATLAB pp-form splines; 2022 is default |
| normalized RH | Eq. (A1) | `T_RH_RHn_2020.m`, `_2022.m` | VERIFIED-CODE | `swim.source.climatological_source_conditions` | direct MATLAB 2022 value + frozen 2020 value | corrected Kelvin conversion is default; `version="2020"` preserves the bug |
| saturation mixing ratio | Eq. (A2) | `mixed_phased_supersaturation.m`; `pseudo_adiabat_function.m` | VERIFIED-CODE | `swim.saturation.saturated_mixing_ratio` | unit test + frozen trajectory | \(\eta=0.622\) preserved exactly |
| saturation vapor pressure, liquid | Appendix A1.2 / Murphy & Koop | `pseudo_adiabat_function.m`; `evaporation_2021.m` | VERIFIED-CODE | `swim.saturation.saturation_vapor_pressure_liquid` | fixed-T unit values + frozen mixed-phase trajectory | exact formula preserved; Celsius input and kPa output |
| saturation vapor pressure, ice | Appendix A1.2 / Murphy & Koop | same | VERIFIED-CODE | `swim.saturation.saturation_vapor_pressure_ice` | fixed-T unit values + frozen mixed-phase trajectory | exact formula preserved; Celsius input and kPa output |
| mixed cloud liquid/ice fraction | Appendix A1.2; Fig. A7 | `fraction_il_brm_H10.m` | VERIFIED-CODE | `swim.cloud_phase.cloud_phase_fractions` | unit tests + frozen full fraction curve | all three legacy fits retained; active method in `distillation_2020.m` is `adj` |
| mixed-phase supersaturation | Appendix A1.2 | `mixed_phased_supersaturation.m` | VERIFIED-CODE | `swim.thermodynamics.mixed_phase_supersaturation` | frozen full trajectory | separate first Euler integration preserved |
| pseudo-adiabatic pressure trajectory | Appendix A1.2 | `mixed_phased_supersaturation.m`; `pseudo_adiabat_function.m` | VERIFIED-CODE | `swim.thermodynamics.pseudo_adiabat` | frozen full P(T) trajectory | second explicit Euler integration preserved |
| fraction of vapor remaining \(f\) | Eq. (A10) | `pseudo_adiabat_function.m` | VERIFIED-CODE | `swim.thermodynamics.PseudoAdiabatResult` | frozen full f(T) trajectory | `f = r_s/r_s(1)` |
| fractionation factor definition | Eq. (A3) | embedded in evaporation/distillation formulas | VERIFIED-CODE | `swim.fractionation.EquilibriumFractionationFactors` | unit tests on α | factors are heavy-phase/light-vapor ratios |
| equilibrium D fractionation, liquid | Appendix A2 | `evaporation_2021.m`; `distillation_2020.m` | VERIFIED-CODE | `swim.fractionation.equilibrium_fractionation_factors` | fixed-T unit values + frozen effective α | Criss coefficients |
| equilibrium D fractionation, ice | Appendix A2; Lamb et al. update | same | VERIFIED-CODE | same | fixed-T unit values + frozen effective α | Lamb formula active |
| equilibrium \(^{18}O\), liquid | Appendix A2 | same | VERIFIED-CODE | same | fixed-T unit values + frozen effective α | Criss coefficients |
| equilibrium \(^{18}O\), ice | Appendix A2 | same | VERIFIED-CODE | same | fixed-T unit values + frozen effective α | Criss coefficients |
| \(^{17}O\) equilibrium scaling | Appendix A2.1 | `evaporation_2021.m`; `distillation_2020.m` | OPEN | explicit exponent argument to `equilibrium_fractionation_factors` | fixed-T unit values + frozen effective α | unresolved choice remains explicit: evaporation uses 0.529; distillation ice uses 0.531 |
| evaporation diffusive fractionation | Eq. (A4), Eq. (A5) | `evaporation_2021.m`, `_2022.m` | VERIFIED-CODE / OPEN-PROVENANCE | `swim.source.initial_vapor_from_climatology` | direct MATLAB 2022 initial vapor + frozen 2021 initial vapor | default 2022 uses N=0.27 and direct HDO ratio; legacy option preserves N=0.302 approximation |
| source seawater \(\delta^{18}O\) | Appendix A2.1 | `evaporation_2021.m` | VERIFIED-CODE | `swim.source.initial_vapor_from_climatology` | frozen initial vapor | default −0.3‰ preserved |
| source seawater \(\delta D\) relation | Appendix A2.1 | `d18Osw_to_dDsw.m` | VERIFIED-CODE | `swim.source.seawater_delta_d_from_delta_18o` | linear-fit unit test + frozen initial vapor | fit is recalculated from checked-in observations |
| local evaporation closure | Eq. (A7) | `evaporation_2021.m` | VERIFIED-CODE | `swim.source.initial_vapor_from_climatology` | frozen T0=10 initial vapor | top-level default |
| global evaporation closure | Eq. (A8) | `evaporation_2021.m` | VERIFIED-CODE | same, `closure="global"` | formula unit coverage pending | sensitivity/end-member configuration |
| supersaturation \(S_i\) | Eq. (A14); Appendix A4 | `pseudo_adiabat_function.m` | VERIFIED-CODE | `swim.thermodynamics.prescribed_supersaturation` | exact S(T) unit test + frozen trajectory | code supports \(a-bT-cT^2\) |
| tuned \(b=0.00525\) | Appendix A4 | `run_SWIM_example.m`; saved result filenames | VERIFIED-PAPER / CODE | default in `swim.model.forward_trajectory` | frozen trajectory and state space | local closure base |
| mixed-phase effective α | Eq. (A13) | `distillation_2020.m` | VERIFIED-CODE | `swim.fractionation.mixed_phase_effective_fractionation` | frozen full αeff(T) curves | weighted ice + liquid |
| kinetic condensation α | Eq. (A12) | `distillation_2020.m` | VERIFIED-CODE | `swim.fractionation.kinetic_condensation_factors` | unit tests + frozen effective α curves | selected liquid kinetic factors remain one |
| transport diffusivity ratios | Appendix A2.2 | `distillation_2020.m`; `distillation_2022.m` | OPEN | `swim.fractionation.transport_diffusivity_ratios` | unit tests + frozen effective α curves | active checked-in versions use HH temperature dependence; publication provenance remains open |
| Rayleigh differential equation | Eq. (A9), Eq. (A11) | `distillation_2020.m` | VERIFIED-CODE | `swim.distillation.rayleigh_distillation` | frozen stepwise ratio trajectories | direct ln(R) update |
| precipitation isotope ratio | Appendix A2.2 | `distillation_2020.m` | VERIFIED-CODE | `swim.distillation.DistillationResult` | frozen Rp trajectories | \(R_p=\alpha R_v\) |
| no precipitation if no loss of vapor | implementation behavior | `distillation_2020.m` | VERIFIED-CODE | `swim.distillation.rayleigh_distillation` | zero-Δf unit case | MATLAB `NaN` behavior preserved |
| complete forward trajectory | Sect. 3; Appendix A | active evaporation + distillation call chain | VERIFIED-CODE | `swim.model.forward_trajectory` | frozen 401-step trajectory with explicit 2021 option | corrected 2022 evaporation is the Python default |
| forward state-space grid | Sect. 3 | `simple_water_isotope_model_2020.m` | LIKELY-ACTIVE | `swim.model.forward_state_space` | frozen 29×71 state space with explicit 2021 option | newly generated Python state spaces default to corrected 2022 evaporation |
| full state-space products | Sect. 3 | `SWIM_results/*.mat` | VERIFIED-CODE/DATA | `swim.model.StateSpace` | frozen Allan Hills generation grid | broader publication-file provenance remains open |

---

## 2. Reconstruction traceability

| Scientific concept | Paper source | MATLAB implementation | Code status | Python target | Test target | Notes |
|---|---|---|---|---|---|---|
| nonlinear state-space inversion | Sect. 4.1 | `Tsite_Tsource_reconstruction_quick.m` | VERIFIED-CODE | `swim.reconstruction.reconstruct_temperatures` | all 2,326 Allan Hills rows | interpolates modeled state space |
| preferred \((\delta^{18}O,d_{ln})\) coordinates | Sect. 4.1; Appendix A6 | method 1 in `Tsite_Tsource_reconstruction_quick.m` | VERIFIED-PAPER/CODE | same | frozen Allan Hills reconstruction | MATLAB log-scaled \(\delta^{18}O\) preserved |
| alternate \((\delta^{18}O,\delta D)\) coordinates | Appendix A6 | method 2 | VERIFIED-CODE | `reconstruct_temperatures(..., method=2)` | smoke coverage | active linear interpolation and inconsistent r_s query preserved |
| alternate \((\delta^{18}O,d_{xs})\) coordinates | Appendix A6 | method 3 | VERIFIED-CODE | `reconstruct_temperatures(..., method=3)` | smoke coverage | active inconsistent raw/log x query and r_s query preserved |
| natural-neighbor interpolation behavior | implementation detail | MATLAB `griddata(...,'natural')` | VERIFIED-CODE | `swim.interpolation.natural_neighbor_interpolate` | frozen Allan Hills reconstruction | Sibson area weights reproduce MATLAB to ~1e-12 °C |
| condensation temperature meaning | Appendix A3.2 | state-space \(T_{\rm site}\) coordinate | VERIFIED-PAPER | `swim.model.StateSpace.condensation_temperature_c` | full state-space + reconstruction parity | weighted condensation temperature, not surface T |
| surface ↔ condensation relation | Appendix A3.2 | `Ts_to_Tc_2020.m` | VERIFIED-PAPER/CODE | `swim.reconstruction.surface_temperature_from_condensation` | frozen Allan Hills reconstruction | \(T_c=0.69T_s-8.2\) |
| moisture-source temperature meaning | Appendix A3.1 | state-space \(T_{\rm source}\) coordinate | VERIFIED-PAPER | `swim.model.StateSpace.source_temperature_c` | full state-space + reconstruction parity | moisture-weighted source temperature, not fixed geographic SST |
| seawater correction of ice-core records | Sect. 4.2; Appendix | `seawater_cor_ln.m`; `reconstruction_2020.m` | LIKELY-ACTIVE | `swim.reconstruction.seawater_correct_isotopes` | initial-age identity unit test | uses checked-in Bintanja ice-volume series |
| publication multi-core application | Sect. 4–5 | `reconstruction_2020.m` + external data compilation | PARTIAL / OPEN | application layer, not core | reproduce one core first | script depends on unavailable author-local paths/files |
| reconstruction uncertainty | Appendix A9 | `Tsite_Tsource_reconstruction_2020_comb_unc.m`, `_ensemble.m`; multiple SWIM result files | PARTIAL | TBD | WDC uncertainty benchmark | map only after base inversion parity |

---

## 3. Tuning and sensitivity traceability

| Scientific choice | Paper source | MATLAB implementation / artifact | Status | Python target | Validation |
|---|---|---|---|---|---|
| tune supersaturation against observed δD–δ18O relationship | Appendix A4 | `tuning_SWIM_2020.m` | LIKELY-ACTIVE | later analysis module | reproduce cost curve / selected b |
| tuned linear \(S_i=1-0.00525T\) | Appendix A4 | example + saved result names | VERIFIED | model config | example trajectory |
| nonlinear \(S_i=a-bT-cT^2\) sensitivity | Appendix A4 | `pseudo_adiabat_function.m`; saved state spaces | VERIFIED-CODE | config option | selected c cases |
| local vs global closure sensitivity | Appendix A2.1, A9 | evaporation routine + multiple SWIM result files | VERIFIED | config option | compare selected state-space points |
| NCEP vs ERA source climatology | Appendix A1.1, A9 | `T_RH_RHn_2020.m`; saved state spaces | VERIFIED | source-condition module | cross-dataset reference values |
| cloud precipitation scheme sensitivity | Appendix A1.2, A9 | alternate pseudo-adiabat / RH routines | PARTIAL | later | only after base parity |
| transport mixing sensitivity | Appendix A5, A9 | not yet fully mapped in this draft | OPEN | later | reproduce one published mixing test |
| evaporation mixing / closure | Appendix A2.1 | local/global and paper mixing experiments | PARTIAL | later | not core parity unless needed |

---

## 4. MATLAB version/provenance traceability

This section exists because the repository contains parallel dated files.

| Logical role | Example / wrapper-selected file | Other variants present | Current assessment |
|---|---|---|---|
| top-level state-space model | `simple_water_isotope_model_2020.m` | `simple_water_isotope_model_2019.m`, `_2019_v2.m` | example points to 2020 |
| evaporation, climatological RH | `evaporation_2021.m` | `evaporation.m`, `_2020`, `_2022`, `_sw` | selected by wrapper |
| evaporation, explicit RH | `evaporation_2020.m` | same | branch-dependent difference |
| distillation | `distillation_2020.m` | base, `_2018`, `_2019`, `_2022` | selected by wrapper |
| source climatology | `T_RH_RHn_2020.m` | `_2022.m` | 2020 is called by evaporation_2021 |
| cloud phase | `fraction_il_brm_H10.m`, method `adj` | `fraction_il_brm.m` | selected by distillation |
| pseudo-adiabat | `pseudo_adiabat_function.m` | `function2`, `alt`, RH variants | selected by distillation |
| inversion | `Tsite_Tsource_reconstruction_quick.m` in `reconstruction_2020.m` | older / uncertainty variants | stronger publication-analysis candidate |
| surface conversion | `Ts_to_Tc_2020.m` | no same-role modern variant seen | active formula matches paper |

---

## 5. Open discrepancy register

Use stable IDs so issues, commits, and tests can refer to these items.

| ID | Issue | Evidence | Risk | Resolution needed before |
|---|---|---|---|---|
| SWIM-D001 | publication baseline versions are not explicitly identified | wrapper uses mixed date suffixes; 2022 variants coexist | very high | exact publication reproduction |
| SWIM-D003 | explicit-RH branch calls `evaporation_2020`, climatology branch calls `_2021` | top-level wrapper | medium/high | RH-dimension support |
| SWIM-D004 | \(^{17}O\) ice equilibrium exponent 0.529 vs 0.531 | evaporation vs distillation | medium | triple-isotope parity |
| SWIM-D005 | transport diffusivity base configuration needs provenance | paper discusses fixed + HH; active code uses HH | medium/high | distillation parity |
| SWIM-D006 | `T_RH_RHn_2022.m` file declares function `T_RH_RHn_2020` | file header | medium | use of 2022 source-condition file |
| SWIM-D007 | older reconstruction helper references stale/missing function names | `Tsite_Tsource_reconstruction.m` | medium | choosing inverse baseline |
| SWIM-D008 | publication reconstruction script has external absolute paths and missing compilation scripts | `reconstruction_2020.m` | high for full paper reproduction, low for core model | full ice-core application |
| SWIM-D010 | MATLAB spline objects are loaded from intentionally untracked legacy data | `data/*_spline_model_*.mat` | medium | standalone package distribution |
| SWIM-D011 | pseudo-adiabat endpoint phase branch tests the penultimate temperature | `pseudo_adiabat_function.m` uses `T(i)` after its loop | low/medium for trajectories ending at 0°C | preserved in Python; reconsider only after parity phase |

Do not close a discrepancy merely because one interpretation is more scientifically plausible. Close it when the selected reference behavior is established and the intended future behavior is separately documented.

---

### Resolved during the baseline port

| ID | Resolution |
|---|---|
| SWIM-D002 | Both behaviors are explicit: corrected `T_RH_RHn_2022.m` conversion is the Python default, while `evaporation_version="2021"` retains the frozen subtraction behavior. |
| SWIM-D009 | `swim.interpolation.natural_neighbor_interpolate` implements Sibson area weights and matches all 2,326 frozen Allan Hills reconstruction rows to about \(10^{-12}\) °C. |

---

## 6. Implemented parity fixture coverage

The checked-in `matlab-port-baseline-v1` fixture covers:

- helper isotope, saturation, cloud-phase, thermodynamic, and fractionation
  calculations;
- source evaporation initial conditions;
- every step of the 10 °C to −30 °C diagnostic trajectory;
- the complete 29 × 71 forward state space; and
- all 2,326 Allan Hills inverse rows, including points between model-grid
  coordinates and observations outside the interpolation domain.

The fixture includes MATLAB v7 files for array-level comparisons, CSV files
for transparent inspection, and JSON generation provenance. Python coverage is
split among `tests/unit/`, `tests/integration/`, and `tests/parity/`.

---

## 8. Maintenance rule

Whenever an agent discovers an important scientific or implementation fact:

1. update the relevant row here;
2. add or update the discrepancy register if needed;
3. link the MATLAB file / paper equation;
4. link the Python implementation once it exists;
5. link a parity test or explain why one is not practical.

The purpose of this file is to prevent scientific understanding from living only in agent conversations.

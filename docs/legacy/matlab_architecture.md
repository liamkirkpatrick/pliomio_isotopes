# MATLAB repository architecture

## Purpose

This document maps the checked-in MATLAB implementation of the Simple Water Isotope Model (SWIM) so that the code can be ported to Python without repeatedly rediscovering the repository.

Paper:

> Markle, B. R., & Steig, E. J. (2022). Improving temperature reconstructions from ice-core water-isotope records. *Climate of the Past*, 18, 1321–1368. https://doi.org/10.5194/cp-18-1321-2022

Repository inspected:

- `bradley-markle/simple_water_isotope_model`
- commit: `4db23b79aab8f2154254d9111f77e70a5ce8427d`

This is a **research-code repository with multiple dated versions**, not a single clean application/package. The architecture below therefore distinguishes the observed example call chain from related historical or experimental files.

---

## 1. High-level repository character

The repository is primarily flat at the top level and mixes:

- forward-model code;
- alternate and historical model versions;
- reconstruction scripts;
- tuning / sensitivity scripts;
- observational and reanalysis data;
- precomputed model state spaces;
- final or near-final analysis outputs.

The root `README.md` contains only a title and the phrase "simple water isotope model"; it does not document execution, dependencies, or version provenance.

This means the first Python task should **not** be "translate every `.m` file." The correct first task is to pin a reproducible MATLAB baseline.

---

## 2. Candidate publication / example execution path

The checked-in example establishes the clearest runnable intent:

```text
run_SWIM_example.m
        |
        v
simple_water_isotope_model_2020.m
        |
        +-----------------------------+
        |                             |
        v                             v
evaporation_2021.m              distillation_2020.m
        |                             |
        v                             +--> fraction_il_brm_H10.m
T_RH_RHn_2020.m                      |
        |                             +--> mixed_phased_supersaturation.m
        v                             |
data/{ncep,era}_data/                 +--> pseudo_adiabat_function.m
                                      |
                                      +--> Rayleigh isotope loop
        \                             /
         \                           /
          +------ endpoint values --+
                    |
                    v
         SWIM isotope state spaces
```

### Important qualification

That graph is the active **climatological-RH branch** of `simple_water_isotope_model_2020.m`.

If the caller explicitly supplies a source-RH grid, the wrapper calls `evaporation_2020.m` instead of `evaporation_2021.m`.

The repository also contains 2022 versions of several lower-level functions. They are **not automatically used** by the checked-in example wrapper.

---

## 3. Top-level example: `run_SWIM_example.m`

This is the clearest user-facing example.

It defines:

- a range of condensation/site temperatures, `Tsite=[lo hi step]`;
- a range of source temperatures, `Tsource=[lo hi step]`;
- optional source RH (`[]` means infer from climatology);
- hemisphere;
- closure assumption;
- reanalysis dataset;
- season;
- supersaturation coefficients `a`, `b`, and `c`.

The example sets:

```matlab
RHsource = [];
SH = 1;
closure = 'local';
reanalysis = 'ncep';
season = 'annual';

a = 1;
b = 0.00525;
c = 0;
```

and calls `simple_water_isotope_model_2020(...)`.

The script also declares `pathway='adiabatic'`, but this variable is **not passed** to the top-level function. `simple_water_isotope_model_2020.m` independently hard-codes the adiabatic path.

---

## 4. Forward-model orchestrator: `simple_water_isotope_model_2020.m`

### Role

This function constructs a grid of source and condensation temperatures, executes one complete evaporation + distillation trajectory for every valid pair, and stores the endpoint precipitation isotope values as state-space arrays.

### Inputs

```text
Tsite       = [Tsite_lo, Tsite_hi, dTsite]
Tsource     = [Tsource_lo, Tsource_hi, dTsource]
RHsource    = [] or [RH_lo, RH_hi, dRH]
a, b, c     = supersaturation parameters
closure     = 'local' or 'global'
reanalysis  = 'ncep' or 'era'
SH          = hemisphere selector
season      = annual / seasonal option where supported
```

### Defaults / forced behavior

- closure: local
- reanalysis: NCEP
- hemisphere: Southern
- season: annual
- NCEP forces annual
- pathway: hard-coded `adiabatic`
- internal trajectory step: `dT = 0.1`

### Core loops

For climatological RH:

```text
for each Tsource (i)
    for each Tsite (j)
        if Tsource < Tsite:
            mark trajectory invalid
        else:
            evaporation
            distillation
            save trajectory endpoint
```

For explicit RH:

```text
for each RH (k)
    for each Tsource (i)
        for each Tsite (j)
            ...
```

### State-space orientation

Observed array convention:

```text
dimension 1: source temperature
dimension 2: condensation/site temperature
dimension 3: explicit source RH, if present
```

Example:

```matlab
d18O_site(i,j)
d18O_site(i,j,k)
```

The Python version should name or encapsulate these coordinates rather than rely only on positional array dimensions.

### Major outputs

- `T_site`
- `T_source`
- `RH_source`
- `d18O_site`
- `dD_site`
- `d18Oln_site`
- `dDln_site`
- `dxs_site`
- `d17O_xs_site`
- `dlnU_site`
- `r_s_site`
- `P_site`

Only the **final point** of each full distillation trajectory is saved into these state-space arrays.

---

## 5. Source-environment module: `T_RH_RHn_2020.m`

### Role

Given source air temperature \(T_0\), returns climatologically related:

- SST;
- RH;
- normalized RH;
- uncertainty/error estimates for the fitted relationships.

### Supported choices

- hemisphere: Southern, Northern, or combined
- reanalysis: NCEP or ERA
- fit method: spline, polynomial, or binned
- season where supported
- optionally regenerate fits rather than load precomputed ones

### Base path

`evaporation_2021.m` calls:

```matlab
T_RH_RHn_2020(T0, SH, reanalysis, 'spline', season, 0)
```

so the normal model loads precomputed spline fits rather than rebuilding them.

### Data dependencies

Examples include:

```text
data/ncep_data/ncep_spline_model_SH.mat
data/ncep_data/ncep_spline_model_NH.mat
data/ncep_data/ncep_spline_model_all.mat

data/era_data/era_spline_model_SH.mat
data/era_data/era_spline_model_NH.mat
data/era_data/era_spline_model_all.mat
```

### Important provenance issue

The checked-in `T_RH_RHn_2020.m` appears to calculate Kelvin temperatures using subtraction (`T0 - 273.15`) in its RHn calculation. `T_RH_RHn_2022.m` changes this to addition, consistent with its °C input documentation.

Do not decide which behavior is "correct for parity" until the exact MATLAB run used to create the selected reference state space is pinned.

---

## 6. Initial evaporation: `evaporation_2021.m`

### Role

Calculates the initial isotopic composition of vapor evaporated from the ocean.

### Inputs / environmental controls

- source air temperature \(T_0\)
- optional SST\(_0\)
- optional RH\(_0\)
- closure assumption
- reanalysis / hemisphere / season
- optional seawater isotope value
- optional diffusivity / uncertainty configuration

When SST and RH are omitted, it calls the climatological environment routine.

### Major scientific responsibilities

- convert \(T_0\) to Kelvin;
- obtain SST, RH, and RHn;
- compute equilibrium fractionation factors;
- establish seawater isotope ratios;
- estimate source seawater \(\delta D\) from \(\delta^{18}O\);
- compute diffusive fractionation;
- apply either local or global closure;
- output initial vapor \(\delta D\), \(\delta^{18}O\), and \(^{17}O\)-excess.

### Key helper

```text
d18Osw_to_dDsw.m
```

maps source seawater \(\delta^{18}O\) to \(\delta D\).

### Version note

The active climatological branch of the top-level wrapper uses `_2021`; the explicit-RH branch uses `_2020`.

---

## 7. Cloud-phase partition: `fraction_il_brm_H10.m`

### Role

Returns temperature-dependent liquid and ice fractions using polynomial/logistic approximations based on Hu et al. (2010).

Supported modes:

```text
'mid'
'iir'
'adj'
```

The active `distillation_2020.m` path uses:

```text
'adj'
```

This file is small and should be one of the earliest components ported and parity-tested.

---

## 8. Mixed-phase supersaturation: `mixed_phased_supersaturation.m`

### Role

Calculates the supersaturation implied by the mixed liquid/ice thermodynamic trajectory.

It:

1. calculates liquid and ice saturation vapor pressures;
2. forms mixed-phase thermodynamic properties;
3. integrates the pseudo-adiabatic pressure trajectory;
4. computes the mixed saturated vapor mixing ratio;
5. returns supersaturation relative to ice:

```matlab
ss = r_s ./ r_s_i;
```

The function shares much thermodynamic machinery with `pseudo_adiabat_function.m`.

### Porting concern

There is duplicated physics between this routine and the pseudo-adiabat routine. During parity, port the duplication faithfully. Consolidating the implementation should happen only after equivalence tests exist.

---

## 9. Pseudo-adiabatic thermodynamics: `pseudo_adiabat_function.m`

### Role

Computes the pressure and vapor-content trajectory for the prescribed temperature path.

### Inputs

- temperature vector
- initial pressure
- liquid/ice fractions
- mixed-phase supersaturation
- supersaturation parameters `a`, `b`, `c`

### Outputs

- `f`: fraction of initial vapor remaining
- `P`: pressure
- `e_s`: saturation vapor pressure
- `r_s`: vapor mixing ratio
- `dPdT`
- `ss`: supersaturation relative to ice
- in some versions, `ss_l`: saturation relative to liquid

### Core numerical behavior

- fixed temperature grid;
- explicit Euler update of pressure;
- temperature-dependent latent heat and heat capacities;
- mixed liquid/ice state;
- prescribed low-temperature supersaturation;
- `f = r_s / r_s(1)`.

### Endpoint branch quirk

After its Euler loop, `pseudo_adiabat_function.m` chooses the final-point
supersaturation branch using `T(i)`. At that point `i` is still the
penultimate MATLAB index, so a trajectory whose last step crosses 0 °C uses
the penultimate temperature for that decision. The Python parity port
preserves this behavior explicitly; it should not be corrected until after
legacy parity is established.

The precise ordering of these calculations may affect parity and should not be "cleaned up" prematurely.

---

## 10. Distillation: `distillation_2020.m`

### Role

Combines thermodynamics and isotope fractionation along one full \(T_0\rightarrow T_c\) pathway.

### Internal flow

```text
make T grid
    |
    v
cloud liquid/ice fractions
    |
    v
thermodynamic pathway and f(T)
    |
    v
equilibrium fractionation factors
    |
    v
kinetic fractionation factors
    |
    v
mixed-phase effective alpha(T)
    |
    v
Euler / stepwise Rayleigh integration in ln(R)
    |
    v
precipitation isotope ratios at each step
    |
    v
delta values, dxs, logarithmic excesses
```

### Key dependencies

Active adiabatic path:

```text
fraction_il_brm_H10.m
mixed_phased_supersaturation.m
pseudo_adiabat_function.m
```

### Rayleigh update

For each isotope at each step, the code calculates:

```text
dlnR = (alpha - 1) * (ln(f_i) - ln(f_i-1))
R_v,i = exp(ln(R_v,i-1) + dlnR)
R_p,i = alpha_i * R_v,i
```

It sets precipitation to `NaN` when \(f\) does not change.

### Important implementation detail

The file contains many commented experimental formulations. The Python port should follow **active statements**, not comments describing older candidate configurations.

---

## 11. Derived isotope coordinates in the wrapper

After each trajectory, `simple_water_isotope_model_2020.m` stores the final isotope values and calculates / stores:

- \(d_{xs}\)
- \(d_{ln}\) using the Uemura coefficients
- \(^{17}O\)-excess
- final `r_s`
- final pressure

The log-isotope arrays use the code's \(1000\ln(R/R_{\rm std})\) scaling.

---

## 12. Saved forward-model state spaces: `SWIM_results/`

The repository includes many large `.mat` files with filenames encoding assumptions, for example:

```text
SWIM_results_202006_ncep_annual_adiabat_local_initsat_
a1_b00525_c0_adiffevap_009_SOURCE_0_28_SITE_m70_0.mat
```

Filename components appear to record:

- date/version;
- reanalysis;
- averaging season;
- pathway;
- closure assumption;
- initial saturation / RH treatment;
- supersaturation coefficients;
- evaporation diffusivity choice;
- source temperature domain;
- site temperature domain.

These saved state spaces are extremely valuable as potential golden-reference data for the Python port.

### Critical next step

Before porting, select one exact central/reference `.mat` state space and identify the code configuration that generated it.

Do not assume the current active wrapper recreates every saved file exactly.

---

## 13. Nonlinear inversion: `Tsite_Tsource_reconstruction_quick.m`

### Role

Maps measured isotope pairs back to source and condensation temperature by interpolating a precomputed forward-model state space.

### Default/preferred coordinate pair

Method 1:

```text
x: transformed/log δ18O
y: dln
```

The function uses:

```matlab
griddata(..., 'natural')
```

to interpolate:

- \(T_c\)
- \(T_0\)
- optionally `r_s`

after removing `NaN` cells from the state space.

### Alternate methods

Method 2:

```text
log δ18O + log δD
```

Method 3:

```text
δ18O + dxs
```

The paper identifies the log-\(\delta^{18}O\) + \(d_{ln}\) representation as the preferred reconstruction coordinate system.

---

## 14. Broader reconstruction scripts

### `reconstruction_2020.m`

This is a publication-analysis style script, not a reusable library function.

It:

- changes into author-specific absolute filesystem paths;
- loads / compiles multiple deep ice-core records using external scripts not present in this repository;
- applies seawater corrections;
- chooses a precomputed SWIM result file;
- runs `Tsite_Tsource_reconstruction_quick.m` for many cores;
- creates figures and further analysis.

This file should **not** be ported wholesale into the core Python package.

Instead, separate:

1. reusable isotope correction logic;
2. reusable inversion logic;
3. project-specific dataset loading;
4. figure/reconstruction application scripts.

### Other reconstruction files

The repository also contains:

- `Tsite_Tsource_reconstruction.m`
- `Tsite_Tsource_reconstruction_2020_comb_unc.m`
- `Tsite_Tsource_reconstruction_2020_ensemble.m`

These represent older, uncertainty, or ensemble workflows. Their exact roles should be mapped when the reconstruction/uncertainty phase is ported.

---

## 15. Condensation-to-surface conversion: `Ts_to_Tc_2020.m`

Implements the paper's Antarctic base relationship:

```text
Tc = 0.69 * Ts - 8.2
```

and the inverse:

```text
Ts = (Tc + 8.2) / 0.69
```

The file contains many commented historical candidate fits. Only the active `a=0.69`, `b=8.2` branch should define parity behavior.

---

## 16. Tuning and uncertainty code

### `tuning_SWIM_2020.m`

Large analysis script associated with choosing / evaluating supersaturation parameterizations against observational isotope relationships.

This should be treated as a **calibration/reproduction workflow**, not a dependency of the core forward model unless inspection proves otherwise.

### Uncertainty/reconstruction files

Relevant files include:

```text
Tsite_Tsource_reconstruction_2020_comb_unc.m
Tsite_Tsource_reconstruction_2020_ensemble.m
```

and multiple saved SWIM state spaces under `SWIM_results/` representing perturbed assumptions.

A clean Python design should eventually separate:

```text
core forward model
state-space generation
inverse reconstruction
calibration/tuning
uncertainty ensemble
paper/application analysis
```

but only after MATLAB behavior has been pinned.

---

## 17. Data layout

### Reanalysis-derived environmental relationships

```text
data/ncep_data/
data/era_data/
```

include both raw / derived climate data and precomputed `.mat` fit objects.

### Isotope observational datasets

Examples include:

```text
data/GNIP/
data/Masson_Delmotte_surface_short_2.txt
data/Uemura_2008_vapor.txt
data/Uemura_2010_vapor.txt
data/vapor_data/
data/Ship_*.csv
data/d17Oexcess_comp.txt
```

### Saved model products

```text
SWIM_results/
```

contains large precomputed model state spaces.

### Other outputs

```text
All_Ice_Core_Reconstructions.xlsx
MSD_results_2019.mat
ice_core_Ts_Tc.mat
```

appear to be higher-level products rather than core model dependencies.

---

## 18. File classification

| Category | Files / examples | Port priority |
|---|---|---:|
| runnable example | `run_SWIM_example.m` | high |
| state-space orchestrator | `simple_water_isotope_model_2020.m` | high |
| source environment | `T_RH_RHn_2020.m` | high |
| evaporation physics | `evaporation_2021.m` | high |
| cloud phase | `fraction_il_brm_H10.m` | high |
| thermodynamic pathway | `mixed_phased_supersaturation.m`, `pseudo_adiabat_function.m` | high |
| isotope distillation | `distillation_2020.m` | high |
| inverse reconstruction | `Tsite_Tsource_reconstruction_quick.m` | high after forward parity |
| surface/condensation relation | `Ts_to_Tc_2020.m` | medium |
| seawater corrections | `d18Osw_to_dDsw.m`, `seawater_cor_ln.m` | medium |
| tuning | `tuning_SWIM_2020.m` | later |
| publication analysis | `reconstruction_2020.m` | later / application layer |
| dated alternatives | `_2018`, `_2019`, `_2022`, `*_alt*` | inspect selectively |
| raw/reference data | `data/` | preserve / document |
| precomputed state spaces | `SWIM_results/` | critical validation fixtures |

---

## 19. Known architectural hazards

### 19.1 Date suffix does not equal active implementation

The top-level `_2020` function calls an `_2021` evaporation routine and a `_2020` distillation routine. Newer `_2022` files exist but are not called by the example.

### 19.2 Branch-dependent version selection

Supplying explicit RH changes which evaporation function is called.

### 19.3 Large amount of dead/commented experimental code

Core files contain many candidate formulas and older parameterizations left in comments. A mechanical code translator may treat these as equally important context.

### 19.4 Duplicated physics

Saturation vapor pressure, fractionation constants, and pseudo-adiabatic calculations appear in multiple files.

Do not deduplicate during the initial parity port.

### 19.5 External path dependencies in application scripts

`reconstruction_2020.m` contains absolute paths under the original author's filesystem and calls data-compilation scripts not present in this repository.

### 19.6 MATLAB-specific interpolation

The inverse model uses `griddata(...,'natural')`. A Python replacement must be tested numerically against MATLAB; merely choosing a roughly similar SciPy interpolator does not guarantee identical behavior.

### 19.7 Saved `.mat` spline objects

The climatology code loads MATLAB spline objects and evaluates them with `fnval`. Decide whether to:

1. reproduce the spline from original data in Python; or
2. export evaluated reference curves / knots from MATLAB.

For parity, the second option may be simpler initially.

---

## 20. Recommended port order

The safest order is dependency-first, with a MATLAB reference test for each stage.

### Stage A — pin reference configuration

- identify one exact baseline state-space `.mat`;
- identify the generating functions / parameters;
- run a very small set of trajectories in MATLAB;
- save intermediate arrays, not only final isotope values.

### Stage B — pure helper physics

Port and test:

1. isotope conversions;
2. saturation vapor pressures;
3. liquid/ice fraction;
4. equilibrium fractionation formulas;
5. supersaturation formula.

### Stage C — environmental source conditions

Port:

```text
T_RH_RHn_*
```

Prefer exporting known climatological curves from MATLAB for initial parity rather than immediately reproducing the whole fit-generation pipeline.

### Stage D — pseudo-adiabat

Port:

```text
mixed_phased_supersaturation.m
pseudo_adiabat_function.m
```

Compare full arrays:

```text
T, P, e_s, r_s, f, ss
```

### Stage E — evaporation

Port the exact selected evaporation baseline and compare:

```text
SST0, RH0, RHn0
δ18O_v0, δD_v0, 17Oxs_v0
```

### Stage F — distillation

Compare full trajectory arrays, especially:

```text
f
alpha_eq
alpha_k
alpha_eff
vapor R values
precipitation R values
δ18O, δD, dln
```

### Stage G — state-space wrapper

Reproduce a small \(T_0\times T_c\) grid before attempting the full ~24k-trajectory state space.

### Stage H — inverse reconstruction

Only after forward parity, port `Tsite_Tsource_reconstruction_quick.m` and test known isotope pairs against MATLAB.

### Stage I — uncertainty and scientific application

Port uncertainty ensembles, ice-core preprocessing, and project-specific analyses after the core forward/inverse model is trustworthy.

---

## 21. Suggested eventual Python boundaries

This is a target architecture, not a requirement for the first parity commit.

```text
src/swim/
    isotopes.py
    saturation.py
    cloud_phase.py
    thermodynamics.py
    source_conditions.py
    evaporation.py
    distillation.py
    model.py
    inversion.py

tests/
    unit/
    parity/
    fixtures/

analysis/
    tuning/
    ice_core_reconstruction/

legacy_matlab/
    ...
```

Dependency direction should remain:

```text
core model  ---> reusable inversion
    ^
    |
analysis / applications
```

The core package should never import application-specific ice-core analysis code.

---

## 22. Immediate archaeology questions to resolve before implementation

1. Which saved `SWIM_results/*.mat` file should be the canonical central/base state space?
2. Which exact versions of `evaporation_*`, `distillation_*`, and `T_RH_RHn_*` generated it?
3. Was the publication run performed before or after the Kelvin conversion correction visible in `T_RH_RHn_2022.m`?
4. Which diffusivity formulation was active for the publication's central model?
5. Was the 0.531 \(^{17}O\)-ice exponent intentional for the publication baseline?
6. Which reconstruction helper/script produced the published ice-core reconstruction files?
7. Can MATLAB reference runs be generated on the user's machine from the checked-in repository without unavailable external files?

These questions should be answered by reference execution and provenance inspection, not by guessing from the newest filename.

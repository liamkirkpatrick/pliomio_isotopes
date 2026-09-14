# SWIM legacy reference baseline

## Purpose

This document defines the MATLAB implementation that will be treated as the **legacy numerical reference** during the Python port of the Simple Water Isotope Model (SWIM).

The purpose of the reference baseline is not to declare that every detail of the MATLAB implementation is scientifically correct. Its purpose is to define one reproducible implementation whose behavior the Python port can initially match.

During the parity phase:

> MATLAB reference behavior is preserved unless a deliberate deviation is documented separately.

Scientific improvements, bug fixes, refactoring, and updated parameterizations should be made only after the corresponding Python component has an established parity baseline.

---

## Scientific reference

Primary model description:

Markle, B. R., & Steig, E. J. (2022). Improving temperature reconstructions from ice-core water-isotope records. *Climate of the Past*, 18, 1321–1368. https://doi.org/10.5194/cp-18-1321-2022

The primary model description is Appendix A of the paper.

Relevant sections include:

- Appendix A1 — environmental trajectory
- Appendix A1.1 — source-region conditions
- Appendix A1.2 — transport
- Appendix A2 — isotope fractionation
- Appendix A2.1 — evaporation from the ocean
- Appendix A2.2 — distillation
- Appendix A2.3 — supersaturation
- Appendix A3 — application to Antarctic ice-core sites
- Appendix A4 — SWIM tuning
- Appendix A6 — reconstruction coordinates
- Appendix A9 — reconstruction uncertainty

The paper describes the intended scientific formulation. The MATLAB code is the numerical reference for initial parity. Where they differ, the difference must be documented rather than silently reconciled.

---

## MATLAB source repository

Repository:

```text
https://github.com/bradley-markle/simple_water_isotope_model
```

Repository snapshot used during the initial archaeology:

```text
commit: 4db23b79aab8f2154254d9111f77e70a5ce8427d
```

This commit should be treated as the initial source snapshot unless the reference configuration is deliberately repinned.

If the legacy MATLAB source is copied into this repository, record the exact source commit here and do not modify the copied reference implementation without documenting why.

---

## Baseline status

**Status: PROVISIONAL — publication-generating configuration not yet fully pinned**

**Initial Python porting baseline: FROZEN as `matlab-port-baseline-v1`**

**Python core port status: COMPLETE for this frozen baseline**

The frozen porting baseline is an executable numerical target for the first
Python implementation. It does not claim to identify the exact configuration
used for every published result. Its fixture is stored at:

```text
tests/fixtures/matlab/port_baseline_v1/allan_hills/
```

It was generated from a clean repository at commit
`5c0102b3b319f302c056e3dcb2fe92083158c2c7` with MATLAB R2026a Update 4 on
Apple silicon. The fixture metadata records the exact source-file hashes,
configuration, toolbox inventory, input checksum, and output checksums.

The Python implementation now reproduces the helper physics, complete
10 °C-to-−30 °C trajectory, complete 29 × 71 state space, and preferred Allan
Hills inversion stored in this fixture. This completion statement is limited
to the frozen baseline: the exact state-space provenance for every published
result and the unavailable publication application inputs remain unresolved.

The checked-in repository contains multiple dated versions of important model components. The numerically newest filename is not necessarily the active implementation.

The checked-in example currently implies the following forward-model call chain for the climatological-relative-humidity case:

```text
run_SWIM_example.m
        |
        v
simple_water_isotope_model_2020.m
        |
        +--> evaporation_2021.m
        |       |
        |       +--> T_RH_RHn_2020.m
        |
        +--> distillation_2020.m
                |
                +--> fraction_il_brm_H10.m
                +--> mixed_phased_supersaturation.m
                +--> pseudo_adiabat_function.m
```

This is currently the strongest candidate for the baseline forward-model execution path, but it must be confirmed against the saved state spaces and publication analysis scripts before being declared final.

---

## Candidate base configuration

The following configuration is supported by both the checked-in example and the base assumptions described in the paper.

```text
hemisphere       Southern Hemisphere
closure          local
reanalysis       NCEP/NCAR
season           annual
transport        pseudo-adiabatic
trajectory dT    0.1 deg C
initial pressure 101.325 kPa

supersaturation:
    a = 1
    b = 0.00525 deg C^-1
    c = 0
```

The checked-in example uses climatologically inferred relative humidity:

```matlab
RHsource = [];
```

and therefore follows the climatological source-condition branch.

This configuration should not be called the final baseline until the publication state-space provenance is verified.

---

## Candidate canonical saved state space

The repository contains precomputed model state spaces under:

```text
SWIM_results/
```

A prominent candidate base file referenced by the publication reconstruction scripts is:

```text
SWIM_results_202006_ncep_annual_adiabat_local_initsat_a1_b00525_c0_adiffevap_009_SOURCE_0_28_SITE_m70_0.mat
```

The exact relative path and file contents should be verified before this file is adopted as the canonical reference fixture.

Once selected, record here:

```text
Canonical state-space file:
    TBD

SHA-256:
    TBD

Generated by MATLAB version:
    TBD

Generated from Git commit:
    TBD

Generation script:
    TBD

Generation date:
    TBD
```

---

## Baseline forward-model API

The candidate top-level MATLAB function is:

```matlab
simple_water_isotope_model_2020(
    Tsite,
    Tsource,
    RHsource,
    a,
    b,
    c,
    closure,
    reanalysis,
    SH,
    season
)
```

The wrapper returns state-space arrays including:

```text
T_site
T_source
RH_source
d18O_site
dD_site
d18Oln_site
dDln_site
dxs_site
d17O_xs_site
dlnU_site
r_s_site
P_site
```

For climatological RH, the output arrays are indexed as:

```text
dimension 1 = source temperature
dimension 2 = condensation/site temperature
```

For explicitly swept RH:

```text
dimension 1 = source temperature
dimension 2 = condensation/site temperature
dimension 3 = source relative humidity
```

The Python implementation should make these dimensions explicit rather than relying only on positional array orientation.

---

## Candidate inverse-model baseline

The strongest candidate for the publication-era inverse path is:

```text
Tsite_Tsource_reconstruction_quick.m
```

as used by:

```text
reconstruction_2020.m
```

For the preferred reconstruction method, the function uses:

```text
log-transformed delta18O
+
logarithmic deuterium excess (dln)
```

and interpolates the precomputed SWIM state space using MATLAB:

```matlab
griddata(..., 'natural')
```

The Python equivalent of this interpolation must be parity-tested. A scientifically similar interpolation method is not automatically a numerical match.

---

## Known unresolved baseline questions

The following items must be resolved before declaring this file final.

### BASE-001 — exact publication state space

Which `SWIM_results/*.mat` file generated the paper's central nonlinear reconstructions?

### BASE-002 — source-condition routine

Did the publication baseline use:

```text
T_RH_RHn_2020.m
```

or a corrected/later source-condition implementation?

The checked-in 2020 function appears to contain a Celsius-to-Kelvin conversion using subtraction in one calculation, whereas a later file changes the conversion to addition.

This should not be silently corrected until reference behavior is established.

### BASE-003 — evaporation implementation

The top-level wrapper calls:

```text
evaporation_2021.m
```

for climatological RH, but:

```text
evaporation_2020.m
```

for explicitly specified RH.

Determine whether this difference was intentional and whether both branches need to be preserved.

### BASE-004 — distillation implementation

Confirm whether:

```text
distillation_2020.m
```

generated the paper's base state space, or whether a later configuration was used.

### BASE-005 — transport diffusivity parameterization

The paper discusses both fixed Jouzel-Merlivat transport diffusivities and temperature-dependent Hellmann-Harvey diffusivities.

Confirm which implementation was active for the canonical state space.

### BASE-006 — triple-oxygen fractionation

The checked-in evaporation and distillation routines do not use identical 17O equilibrium exponents for ice.

Confirm the publication configuration before changing either behavior.

### BASE-007 — MATLAB environment

Record:

```text
MATLAB release:
required toolboxes:
operating system:
```

if known.

This is especially important for interpolation and spline behavior.

---

## Criteria for freezing the baseline

The baseline may be changed from **PROVISIONAL** to **FROZEN** when all of the following are true:

1. The canonical MATLAB source commit is recorded.
2. The exact forward-model function versions are recorded.
3. The canonical parameter set is recorded.
4. The canonical saved state-space file is identified.
5. At least one MATLAB run reproduces selected values from that state space.
6. At least one complete trajectory has been exported with intermediate values.
7. The inverse interpolation method is identified.
8. Known paper/code discrepancies are recorded in `docs/legacy/known_quirks.md` or `docs/porting/traceability.md`.
9. Reference fixtures have checksums and provenance metadata.

---

## Baseline change policy

Once frozen, the baseline should not be changed casually.

If a different MATLAB configuration becomes the preferred reference:

1. preserve the old fixtures;
2. document why the baseline changed;
3. assign a new baseline identifier;
4. regenerate fixtures under a new directory;
5. do not overwrite historical reference outputs.

Suggested identifiers:

```text
matlab-baseline-v1
matlab-baseline-v2
```

Each parity test should make clear which baseline it targets.

---

## Current provisional baseline summary

Until the archaeology is complete, agents should assume:

```text
Repository commit:
    4db23b79aab8f2154254d9111f77e70a5ce8427d

Top-level wrapper:
    simple_water_isotope_model_2020.m

Climatological evaporation:
    evaporation_2021.m

Source climatology:
    T_RH_RHn_2020.m

Distillation:
    distillation_2020.m

Cloud phase:
    fraction_il_brm_H10.m
    method = 'adj'

Thermodynamic pathway:
    mixed_phased_supersaturation.m
    pseudo_adiabat_function.m

Default scientific configuration:
    local closure
    NCEP annual
    Southern Hemisphere
    pseudo-adiabatic
    a = 1
    b = 0.00525
    c = 0
    dT = 0.1 deg C
```

**Important:** this summary is a working hypothesis for archaeology and test generation, not permission to silently resolve the open discrepancies above.

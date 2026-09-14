# MATLAB reference fixtures

## Purpose

This directory contains fixed numerical outputs generated from the legacy MATLAB implementation of SWIM.

These fixtures are the primary numerical reference for MATLAB-to-Python parity tests.

They are **test data**, not analysis outputs and not manually curated model results.

---

## Core rule

> Reference fixtures must reflect the selected MATLAB baseline exactly.

Do not edit a fixture simply because the Python implementation disagrees with it.

A parity failure should trigger investigation of the Python implementation, the fixture provenance, or the baseline definition.

If a fixture is genuinely wrong or generated from the wrong MATLAB configuration, replace it only through the documented fixture-generation workflow and record the reason.

---

## Provenance requirement

Every fixture must have enough metadata to determine:

- MATLAB source repository
- MATLAB source commit
- baseline identifier
- MATLAB function versions
- model configuration
- model inputs
- MATLAB version, if known
- generation script
- generation date
- file checksum

The preferred way to store this is a small metadata file adjacent to the fixture.

Example:

```text
tests/fixtures/matlab/
    baseline_v1/
        helper_physics/
            reference.npz
            metadata.yml

        source_evaporation/
            T0_10_local/
                reference.npz
                metadata.yml

        trajectories/
            T0_10_Tc_m30/
                reference.npz
                metadata.yml

        state_spaces/
            tiny_3x3/
                reference.npz
                metadata.yml

        inversion/
            interior_points/
                reference.npz
                metadata.yml
```

---

## Recommended fixture hierarchy

```text
tests/fixtures/
    README.md

    matlab/
        baseline_v1/

            helper_physics/
                saturation/
                cloud_phase/
                fractionation/

            source_conditions/
                climatology/

            evaporation/
                local/
                global/

            trajectories/
                T0_10_Tc_m30/

            state_spaces/
                tiny_3x3/

            inversion/
                preferred_coordinates/

    synthetic/
```

`synthetic/` may contain hand-constructed inputs used for general unit tests.

`matlab/` should contain only outputs generated from MATLAB.

---

## Baseline identifier

All initial fixtures should use a single baseline identifier once the legacy reference configuration is frozen.

Suggested initial identifier:

```text
matlab-baseline-v1
```

The publication-generating baseline remains provisional. The separately named
`port_baseline_v1` directory is frozen as the initial executable target for
Python parity and must not be described as the confirmed publication baseline.

Before that point, temporary archaeology exports should live outside the committed canonical fixture directory or be explicitly labeled provisional.

The initial frozen porting fixture is:

```text
tests/fixtures/matlab/port_baseline_v1/allan_hills/
```

The post-parity corrected evaporation reference is:

```text
tests/fixtures/matlab/evaporation_2022_Tsource_10.json
```

It records a direct MATLAB `evaporation_2022.m` run and does not replace or
modify the frozen 2021-compatible fixture.

---

## Preferred storage formats

Where practical, prefer portable formats.

### `.npz`

Good default for numerical arrays:

```text
reference.npz
```

Advantages:

- preserves array shapes and numeric types;
- loads directly with NumPy;
- compact;
- easy to compare in pytest.

### `.csv`

Useful for small human-readable tables such as:

- saturation vapor pressure reference values;
- fractionation factors at selected temperatures;
- source-condition lookup tables.

### `.mat`

Keep the original MATLAB export where it adds value, especially during archaeology, but avoid making Python tests depend exclusively on proprietary MATLAB object serialization.

For a canonical case it may be useful to keep both:

```text
reference.mat
reference.npz
```

where the `.npz` file is generated directly from the `.mat` export by a documented conversion script.

---

## Floating-point precision

Do not round fixture arrays for convenience.

Store numerical arrays at the precision produced by MATLAB unless there is a specific reason not to.

Human-readable summaries may contain rounded values, but parity tests should use the full-precision fixture.

---

## Initial fixture set

The first fixture suite should be small and diagnostic.

### 1. Helper physics

Evaluate several temperatures spanning the model range, for example:

```text
20
10
0
-10
-30
-50 deg C
```

Export:

- saturation vapor pressure over liquid
- saturation vapor pressure over ice
- liquid fraction
- ice fraction
- equilibrium fractionation factors
- prescribed supersaturation

### 2. Source conditions and evaporation

Example source temperatures:

```text
0
5
10
20
28 deg C
```

Export:

- SST0
- RH0
- RHn0
- initial vapor delta18O
- initial vapor deltaD
- initial vapor 17O excess

Include both:

```text
local closure
global closure
```

where possible.

### 3. Complete trajectory

Recommended first full trajectory:

```text
T0 = 10 deg C
Tc = -30 deg C
dT = 0.1 deg C
```

Export every temperature step of all practical intermediate quantities, including:

```text
T
P
e_s
r_s
f
supersaturation
fraction_liquid
fraction_ice

equilibrium fractionation factors
kinetic fractionation factors
effective fractionation factors

vapor isotope ratios
precipitation isotope ratios

delta18O
deltaD
delta17O
log-transformed isotope values
dxs
dln
17O excess
```

The purpose of this fixture is to identify the **first point of divergence** between MATLAB and Python.

### 4. Tiny state space

Recommended grid:

```text
Tsource = [5, 10, 15] deg C
Tsite   = [-40, -30, -20] deg C
```

Use the same wrapper behavior and output orientation as the canonical MATLAB baseline.

Export all state-space arrays.

### 5. Inverse interpolation

Choose several isotope pairs that are inside the model state-space domain, including points that do **not** fall exactly on forward-model grid points.

Export:

- input delta18O
- input deltaD
- transformed delta18O
- dln
- reconstructed Tsource
- reconstructed Tcond
- reconstructed Tsurface if applicable

These cases should be used to validate the Python replacement for MATLAB `griddata(..., 'natural')`.

---

## Metadata example

Suggested `metadata.yml`:

```yaml
fixture_id: trajectory_T0_10_Tc_m30

baseline:
  id: matlab-baseline-v1
  repository: bradley-markle/simple_water_isotope_model
  commit: 4db23b79aab8f2154254d9111f77e70a5ce8427d

matlab:
  version: unknown
  generation_script: scripts/matlab/export_reference_case.m

configuration:
  closure: local
  reanalysis: ncep
  hemisphere: southern
  season: annual
  pathway: adiabatic
  supersaturation:
    a: 1.0
    b: 0.00525
    c: 0.0

inputs:
  T0_C: 10.0
  Tc_C: -30.0
  dT_C: 0.1

notes:
  - Generated before any changes to legacy MATLAB behavior.
```

---

## Checksums

Canonical fixture files should have recorded SHA-256 checksums.

A simple manifest may live at:

```text
tests/fixtures/matlab/baseline_v1/SHA256SUMS
```

This makes accidental modification obvious.

---

## Relationship to parity tests

Parity tests should be located under:

```text
tests/parity/
```

A parity test should:

1. load a committed fixture;
2. run the corresponding Python calculation;
3. compare meaningful intermediate or final arrays;
4. use a documented tolerance;
5. report the first or largest useful difference.

Numerical tolerances belong in:

```text
docs/validation/numerical_tolerances.md
```

Do not choose tolerances independently in individual tests unless necessary.

---

## Regenerating fixtures

Canonical fixtures should only be regenerated when:

- the baseline is intentionally changed;
- a fixture is proven not to represent the stated baseline;
- the fixture format changes without changing numerical content;
- additional diagnostic outputs are needed.

Fixture regeneration must not be used as a way to make Python tests pass.

When the baseline changes, preserve the previous version:

```text
baseline_v1/
baseline_v2/
```

rather than overwriting it.

---

## What does not belong here

Do not store the following as parity fixtures unless explicitly needed for a test:

- publication figures;
- exploratory notebook outputs;
- complete raw climate datasets;
- user-specific ice-core analysis outputs;
- large intermediate datasets that can be regenerated cheaply;
- arbitrary snapshots without provenance.

The fixture suite should be **small, diagnostic, immutable, and reproducible**.

# Simple Water Isotope Model (SWIM) — Python Reimplementation

## Overview

This repository is a Python reimplementation and future extension of the
**Simple Water Isotope Model (SWIM)** described by Markle and Steig (2022).

The project has two major phases:

1. **Reproduce the legacy MATLAB implementation faithfully in Python.**
2. **After numerical and scientific parity is established, improve the software
   design and extend the model for new scientific applications.**

During the porting phase, preserving the behavior of the selected MATLAB
reference implementation is more important than refactoring, optimization, or
scientific modification.

Primary scientific reference:

> Markle, B. R., & Steig, E. J. (2022). Improving temperature reconstructions
> from ice-core water-isotope records. *Climate of the Past*, 18, 1321–1368.  
> https://doi.org/10.5194/cp-18-1321-2022

Legacy MATLAB repository:

> https://github.com/bradley-markle/simple_water_isotope_model

---

## Project goals

The immediate goal is to create a transparent, testable Python implementation
of SWIM that can reproduce the legacy MATLAB model.

The longer-term goals are to:

- make the scientific calculations easier to understand and maintain;
- improve numerical testing and reproducibility;
- separate reusable model code from dataset-specific analysis;
- preserve clear traceability from the paper to MATLAB to Python;
- support future changes to model assumptions and parameterizations;
- apply the model to new water-isotope datasets.

Scientific changes should be clearly distinguished from the initial MATLAB
port.

---

## Repository structure

The planned repository structure is:

```text
project/
├── README.md
├── AGENTS.md
├── LICENSE
├── pyproject.toml
├── environment.yml              # optional if Conda is used
├── .gitignore
│
├── src/
│   └── swim/
│       ├── __init__.py
│       ├── isotopes.py
│       ├── saturation.py
│       ├── cloud_phase.py
│       ├── thermodynamics.py
│       ├── source_conditions.py
│       ├── evaporation.py
│       ├── distillation.py
│       ├── model.py
│       ├── inversion.py
│       └── temperature.py
│
├── tests/
│   ├── unit/
│   ├── integration/
│   ├── parity/
│   └── fixtures/
│       ├── README.md
│       ├── matlab/
│       │   └── baseline_v1/
│       └── synthetic/
│
├── docs/
│   ├── science/
│   │   ├── model_overview.md
│   │   ├── equations.md
│   │   ├── assumptions.md
│   │   └── terminology.md
│   │
│   ├── legacy/
│   │   ├── matlab_architecture.md
│   │   ├── execution_flow.md
│   │   ├── function_map.md
│   │   ├── known_quirks.md
│   │   └── reference_baseline.md
│   │
│   ├── porting/
│   │   ├── strategy.md
│   │   ├── traceability.md
│   │   └── decisions/
│   │       └── README.md
│   │
│   └── validation/
│       ├── reference_cases.md
│       └── numerical_tolerances.md
│
├── scripts/
│   ├── matlab/
│   │   ├── export_reference_case.m
│   │   └── README.md
│   └── python/
│
├── analysis/
│   ├── exploratory/
│   ├── notebooks/
│   └── figures/
│
├── applications/
│   ├── ice_core_reconstruction/
│   └── future_studies/
│
├── data/
│   └── README.md
│
└── legacy_matlab/
    └── ...
```

Not every file or module needs to exist immediately. The structure is intended
to provide clear boundaries as the project grows.

---

## Directory responsibilities

### `src/swim/`

Contains the reusable Python implementation of the scientific model.

The core model package should contain generally useful scientific
functionality, not project-specific plotting, notebooks, or individual
dataset workflows.

Expected eventual responsibilities include:

```text
isotopes.py           isotope-ratio conversions and excess parameters
saturation.py         saturation vapor pressure calculations
cloud_phase.py        liquid/ice cloud partitioning
thermodynamics.py     pseudo-adiabatic trajectory calculations
source_conditions.py  source SST/RH climatological relationships
evaporation.py        initial ocean evaporation fractionation
distillation.py       Rayleigh distillation and isotope fractionation
model.py              forward state-space generation
inversion.py          nonlinear isotope-to-temperature reconstruction
temperature.py        condensation/surface temperature relationships
```

These modules should be created as their functionality is ported rather than
as empty placeholders.

---

### `tests/`

Contains automated tests.

```text
tests/unit/
```

Tests small calculations and individual functions.

```text
tests/integration/
```

Tests interactions among multiple model components.

```text
tests/parity/
```

Compares the Python implementation directly with frozen outputs from the
legacy MATLAB model.

```text
tests/fixtures/
```

Contains immutable reference inputs and outputs used by tests.

MATLAB-generated fixtures should include provenance describing the exact
legacy commit, configuration, and script used to generate them.

---

### `docs/science/`

Documents the science independently of a particular programming language.

Key files:

```text
model_overview.md
```

High-level description of what SWIM represents and how the major physical
processes fit together.

```text
equations.md
```

Equations, units, parameterizations, and their relationship to the MATLAB
implementation.

```text
assumptions.md
```

Scientific simplifications and structural assumptions of the model.

```text
terminology.md
```

Definitions of model-specific notation, abbreviations, and commonly used
variables.

---

### `docs/legacy/`

Documents the existing MATLAB implementation.

```text
matlab_architecture.md
```

Major MATLAB components and their responsibilities.

```text
execution_flow.md
```

The sequence of calculations during one forward model run.

```text
function_map.md
```

Mapping of MATLAB functions to their role, dependencies, Python replacement,
and porting status.

```text
known_quirks.md
```

Known implementation oddities, possible bugs, historical behavior, and
paper/code discrepancies.

```text
reference_baseline.md
```

Defines the exact MATLAB configuration used as the numerical reference for
Python parity.

---

### `docs/porting/`

Documents decisions associated with the MATLAB-to-Python transition.

```text
strategy.md
```

Overall porting sequence and rules.

```text
traceability.md
```

Maintains mappings among:

```text
scientific concept
    ↓
paper equation or section
    ↓
MATLAB implementation
    ↓
Python implementation
    ↓
validation / parity test
```

```text
decisions/
```

Stores records of deliberate implementation or scientific decisions made
during the port.

---

### `docs/validation/`

Documents how numerical correctness is established.

```text
reference_cases.md
```

Defines the canonical MATLAB runs used for comparison.

```text
numerical_tolerances.md
```

Records justified numerical tolerances for parity and regression tests.

Tolerances should not be loosened solely to make failing tests pass.

---

### `legacy_matlab/`

Contains the frozen MATLAB reference implementation if the legacy code is
copied into this repository.

The legacy source should remain unchanged during normal Python development.

Any instrumentation required to export additional reference values should be
minimal and documented.

---

### `scripts/matlab/`

Contains MATLAB utilities used specifically for validation and archaeology.

For example:

```text
export_reference_case.m
```

should generate complete numerical reference cases, including intermediate
states needed to diagnose MATLAB/Python divergence.

These scripts are separate from the frozen legacy implementation.

---

### `analysis/`

Contains exploratory scientific work.

Examples include:

- notebooks;
- diagnostic plots;
- sensitivity exploration;
- ad hoc comparisons;
- figure development.

Code that becomes generally useful should eventually move into the appropriate
reusable package or application layer.

---

### `applications/`

Contains reproducible scientific uses of SWIM on real datasets.

For example:

```text
applications/ice_core_reconstruction/
```

may eventually contain the data-processing and reconstruction workflow used
to apply the Python model to ice-core isotope records.

Application-specific behavior should not be placed inside the core model
package unless it represents reusable model functionality.

---

### `data/`

Documents the expected location and provenance of scientific datasets.

Large raw or derived datasets generally should not be committed directly to
Git.

A `data/README.md` should explain:

- where datasets originate;
- how to obtain them;
- which data are raw versus derived;
- expected local directory structure;
- preprocessing steps;
- which small files are intentionally tracked.

Small synthetic datasets and numerical test fixtures belong under `tests/`
instead.

---

## Architectural boundary

The reusable scientific model lives under:

```text
src/swim/
```

Analysis and applications may depend on that package:

```text
analysis/ --------                   applications/ -----> src/swim/
                   /
tests/ ------------/
```

The reverse dependency should not occur.

In particular:

```text
src/swim/
```

should not import code from:

```text
analysis/
applications/
```

The model should not contain functions that exist only to produce a particular
paper figure or operate on one particular dataset.

---

## Porting philosophy

### Establish parity before improvement

The first Python implementation should reproduce the selected MATLAB baseline.

Do not silently:

- fix suspected scientific errors;
- replace numerical algorithms;
- change fractionation parameterizations;
- alter interpolation behavior;
- change array conventions;
- simplify apparently redundant calculations;
- refactor behavior before it is tested.

If MATLAB behavior appears wrong or inconsistent with the publication,
document the discrepancy first.

The appropriate behavior during the parity phase is determined by the frozen
reference baseline.

---

## MATLAB reference baseline

The repository contains several dated MATLAB implementations of major model
components.

Therefore filenames alone are not sufficient to determine which version
represents the publication model.

The current baseline investigation is documented in:

```text
docs/legacy/reference_baseline.md
```

That document should eventually record:

- legacy Git commit;
- exact MATLAB function versions;
- model configuration;
- canonical saved SWIM state-space file;
- MATLAB version and relevant toolboxes;
- generation scripts;
- known unresolved discrepancies.

Until that baseline is frozen, MATLAB version provenance should be treated as
an active archaeology task.

---

## Validation strategy

The Python implementation will be validated incrementally.

The preferred workflow is:

```text
1. Select and freeze the MATLAB reference configuration.
        ↓
2. Generate MATLAB reference fixtures.
        ↓
3. Port one small scientific component.
        ↓
4. Compare Python against MATLAB.
        ↓
5. Identify the first numerical divergence.
        ↓
6. Resolve the difference.
        ↓
7. Continue to the next component.
```

Reference fixtures should include intermediate model states, not only final
outputs.

A full trajectory fixture should ideally include quantities such as:

```text
temperature
pressure
saturation vapor pressure
water-vapor mixing ratio
fraction of vapor remaining
supersaturation
ice/liquid fractions
equilibrium fractionation factors
kinetic fractionation factors
effective fractionation factors
vapor isotope ratios
precipitation isotope ratios
delta values
isotope-excess parameters
```

This makes parity failures much easier to localize.

---

## Current project phase

The project is currently in the **model archaeology and validation setup**
phase.

Current priorities are:

1. document the scientific formulation;
2. map the MATLAB architecture;
3. determine the exact publication/reference implementation;
4. record known paper/code discrepancies;
5. generate MATLAB reference outputs;
6. only then begin the Python port.

The Python package should remain intentionally small until the reference
behavior has been established.

---

## Development setup

The package configuration is defined in:

```text
pyproject.toml
```

A typical editable development installation will eventually be:

```bash
python -m pip install -e ".[dev]"
```

Run tests with:

```bash
pytest
```

Run only MATLAB parity tests with:

```bash
pytest -m parity
```

If Conda is used, environment setup should be documented in
`environment.yml` and this section should be updated with the exact commands.

---

## Documentation expectations

Important scientific or implementation discoveries should be written back
into the repository.

Do not leave critical knowledge only in:

- agent conversations;
- issue comments;
- commit messages;
- notebooks;
- personal notes.

When new behavior is discovered, update the appropriate document under
`docs/`.

Unresolved interpretations should be clearly identified as unresolved rather
than documented as fact.

---

## Citation

If this project is used scientifically, cite the original SWIM publication:

> Markle, B. R., & Steig, E. J. (2022). Improving temperature reconstructions
> from ice-core water-isotope records. *Climate of the Past*, 18, 1321–1368.
> https://doi.org/10.5194/cp-18-1321-2022

Additional citation information for the Python implementation should be added
once the software is released.

---

## License

TBD.

The license for this Python repository should be chosen explicitly before
public release. The licensing status of copied or redistributed legacy MATLAB
code and bundled datasets should also be checked separately.

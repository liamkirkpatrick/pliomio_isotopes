# Simple Water Isotope Model (SWIM)

This repository contains a parity-first Python reimplementation of the Simple
Water Isotope Model described by Markle and Steig (2022), together with the
frozen MATLAB source used as its numerical reference.

The port currently reproduces the frozen `matlab-port-baseline-v1` workflow:

- source-region climatology and ocean evaporation;
- mixed-phase pseudo-adiabatic transport;
- equilibrium and kinetic isotope fractionation;
- Rayleigh distillation and isotope diagnostics;
- forward trajectory and state-space generation;
- natural-neighbor isotope-to-temperature reconstruction; and
- the surface/condensation-temperature and seawater corrections used by the
  legacy reconstruction path.

The implementation deliberately preserves known MATLAB quirks where they
affect parity. Scientific corrections and broader extensions should be made
after the relevant baseline behavior has been documented and tested.

The Python forward model now defaults to the corrected MATLAB
`evaporation_2022.m` behavior. Pass `evaporation_version="2021"` to reproduce
the frozen legacy state-space fixture exactly.

## Scope

“Port complete” in this repository means the frozen core forward and inverse
model can be run and matches its checked-in MATLAB fixture. It does not mean
that every historical MATLAB variant or every publication analysis script has
been translated. In particular, the publication-wide multi-core workflow and
uncertainty ensembles depend on author-local inputs that are not present here.
Those remain later application work, separate from the core model.

The baseline and unresolved provenance questions are recorded in
`docs/legacy/reference_baseline.md` and `docs/porting/traceability.md`.

## Setup

Python 3.11 or newer is required. In the existing Conda environment:

```bash
conda activate mioplio
python -m pip install -e ".[dev]"
```

NumPy and SciPy are the only runtime dependencies. Pytest, Ruff, and mypy are
included in the `dev` extra. The larger `environment.yml` also contains the
interactive analysis tools used in this workspace, but they are not required
by the model package.

## Quick start

Run one baseline trajectory:

```python
from swim.model import forward_trajectory

trajectory = forward_trajectory(10.0, -30.0)
print(trajectory.distillation.delta_18o_precipitation[-1])
print(trajectory.distillation.delta_d_precipitation[-1])
```

The call above uses corrected 2022 source humidity and evaporation. For an
exact frozen-baseline comparison:

```python
legacy_trajectory = forward_trajectory(
    10.0,
    -30.0,
    evaporation_version="2021",
)
```

Generate a state space:

```python
import numpy as np

from swim.model import forward_state_space

state_space = forward_state_space(
    np.arange(0.0, 29.0, 1.0),
    np.arange(-70.0, 1.0, 1.0),
)
```

The source-climatology calculation reads the original spline and observational
files under `legacy_matlab/data/`. Those files are intentionally not tracked by
this repository's current Git policy. A fresh checkout therefore needs a local
copy of the upstream legacy data to regenerate a state space. The checked-in
MATLAB state-space fixture is self-contained and can be used for reconstruction
without those local source files.

Run the complete Allan Hills reconstruction with corrected 2022 evaporation:

```bash
python scripts/python/run_allan_hills.py
```

By default, this generates a new state space with corrected 2022 evaporation.
To reconstruct with the frozen 2021-compatible MATLAB fixture instead:

```bash
python scripts/python/run_allan_hills.py \
  --state-space tests/fixtures/matlab/port_baseline_v1/allan_hills/state_space.mat
```

Output is written to `tests/generated/python_allan_hills.csv`; generated runs
are ignored by Git. Generating the corrected state space requires the local
legacy climatology data described above.

Use the model directly for inversion:

```python
from swim.model import load_matlab_state_space
from swim.reconstruction import reconstruct_temperatures

state_space = load_matlab_state_space(
    "tests/fixtures/matlab/port_baseline_v1/allan_hills/state_space.mat"
)
result = reconstruct_temperatures(delta_d, delta_18o, state_space)
```

## Validation

Run everything:

```bash
conda run -n mioplio python -m pytest -q
conda run -n mioplio ruff check src tests scripts
conda run -n mioplio mypy src scripts
```

Run only MATLAB parity tests:

```bash
conda run -n mioplio python -m pytest -q tests/parity
```

The full state-space test explicitly selects 2021 evaporation and regenerates
the 29 × 71 frozen baseline grid. It is the slowest test. A separate direct
MATLAB reference validates the corrected 2022 initial evaporation. Fixtures,
their provenance, and numerical tolerances are documented under
`tests/fixtures/` and `docs/porting/traceability.md`.

To regenerate the MATLAB reference locally:

```matlab
run('tests/run_swim_allan_hills_test.m')
```

That runner writes timestamped output beneath `tests/generated/`; reference
fixtures should never be overwritten merely because Python differs.

## Repository layout

```text
src/swim/                   Python model and reconstruction package
tests/unit/                 focused behavior tests
tests/integration/          end-to-end Python workflows
tests/parity/               comparisons with frozen MATLAB output
tests/fixtures/             immutable MATLAB reference outputs
legacy_matlab/              frozen legacy source and local data
docs/science/               scientific equations and assumptions
docs/legacy/                MATLAB architecture and baseline provenance
docs/porting/               implementation traceability and discrepancies
scripts/python/             runnable Python application scripts
```

The central Python modules are:

```text
isotopes.py        isotope notation and excess diagnostics
saturation.py      saturation vapor pressure and mixing ratio
cloud_phase.py     liquid/ice cloud partitioning
thermodynamics.py  supersaturation and pseudo-adiabatic trajectories
fractionation.py   equilibrium, kinetic, and effective fractionation
source.py          source climatology and evaporation
distillation.py    Rayleigh isotope distillation
model.py           complete trajectories and state spaces
interpolation.py   MATLAB-compatible Sibson natural neighbors
reconstruction.py  isotope inversion and temperature corrections
```

## Porting rules

Scientific behavior is resolved in this order: the paper-derived science
documentation, legacy MATLAB, recorded porting decisions, Python code, then
tests and fixtures. A discrepancy must be documented before scientific
behavior changes. During parity work, correctness and traceability take
priority over refactoring or optimization.

See `AGENTS.md` for the detailed development rules.

## Citation

Markle, B. R., & Steig, E. J. (2022). Improving temperature reconstructions
from ice-core water-isotope records. *Climate of the Past*, 18, 1321–1368.
https://doi.org/10.5194/cp-18-1321-2022

## License

No license has yet been selected for this Python repository. The redistribution
status of the legacy MATLAB source and data should be resolved before public
release.

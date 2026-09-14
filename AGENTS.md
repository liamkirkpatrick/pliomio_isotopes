# AGENTS.md

## Repository-specific details

- Python package: `swim` under `src/swim/`
- Full tests: `conda run -n mioplio python -m pytest -q`
- MATLAB parity tests: `conda run -n mioplio python -m pytest -q tests/parity`
- Lint: `conda run -n mioplio ruff check src tests scripts`
- Type check: `conda run -n mioplio mypy src scripts`
- MATLAB baseline runner: `run('tests/run_swim_allan_hills_test.m')`
- Frozen fixtures: `tests/fixtures/matlab/port_baseline_v1/allan_hills/`
- Fixture formats: numeric MATLAB v7 `.mat`, CSV, and JSON provenance
- Numerical tolerances: use the explicit tolerances in `tests/parity/`; do not
  change them without documenting numerical justification
- No subdirectories currently require a more specialized `AGENTS.md`

--- 

## Project purpose

This repository is a Python reimplementation and future extension of Simple Water Isotope Model (SWIM) an existing scientific modeling framework originally implemented in MATLAB and described in Markle, Bradley R., and Eric J. Steig. "Improving temperature reconstructions from ice-core water-isotope records." Climate of the Past 18.6 (2022): 1321-1368.

The project has two distinct goals:

1. **Reproduce the published / legacy MATLAB model faithfully in Python.**
2. **Only after parity is established, improve the software design and extend the science. Here we will attempt to apply SWIM to produce a first-of-it's-kind reconstrunction of Pliocene and Miocene water isotopes.**

During the porting phase, preserving scientific and numerical behavior is more important than making the implementation elegant. Note that the current script often contains multiple versions of the same function or dataset - I will want to preserve this rough form of version control (e.g. be able to select older versions), but the outdated versions should more explicitly be documented as outdated and the differences should be readily accessible.

---

## Sources of truth

When determining intended model behavior, consult sources in this order:

1. **Scientific documentation derived from the paper and supplement**
   - `docs/science/`
2. **Original MATLAB implementation**
   - `legacy_matlab/`
3. **Recorded project decisions and known discrepancies**
   - `docs/porting/`
4. **Current Python implementation**
   - `src/`
5. **Tests and MATLAB reference outputs**
   - `tests/`

The paper describes the intended science, while the MATLAB code represents the behavior of the legacy implementation. These may not always agree.

**Never silently resolve a discrepancy between the paper, supplement, MATLAB code, documentation, and Python implementation.** Identify it explicitly and document it before changing scientific behavior.

---

## Required background before changing model behavior

Before modifying scientific calculations, read the relevant project documentation.

At minimum, consult:

- `docs/science/model_overview.md`
- `docs/science/equations.md`
- `docs/legacy/matlab_architecture.md`
- `docs/legacy/execution_flow.md`
- `docs/porting/traceability.md`

Also inspect any topic-specific documentation relevant to the code being changed.

If one of these files does not yet exist, do not invent its contents. Use the paper, supplement, MATLAB implementation, tests, and existing documentation as appropriate, and create or update documentation when useful discoveries are made.

---

## Core porting principles

### 1. Establish parity before improving the model

During the MATLAB-to-Python port:

- Preserve observable MATLAB behavior unless a deliberate deviation has been documented.
- Do not fix suspected scientific bugs, numerical quirks, unusual conventions, or awkward algorithms merely because they appear incorrect or inelegant. Do flag these for review.
- Do not combine a behavioral port with a scientific change.
- Prefer a faithful, testable implementation first; refactoring can happen after parity is established.

If the MATLAB code appears to contain a bug or conflicts with the paper:

1. Identify the discrepancy.
2. Determine whether it affects model output.
3. Document it.
4. Preserve legacy behavior during the parity phase unless the project explicitly decides otherwise.

### 2. Preserve scientific meaning, not MATLAB syntax

The Python implementation does not need to mimic MATLAB line-for-line.

Prefer clear Python and NumPy constructs when they preserve the same behavior.

Be especially careful with:

- MATLAB 1-based indexing versus Python 0-based indexing
- inclusive MATLAB ranges
- column-major array assumptions
- implicit array expansion / broadcasting
- matrix multiplication versus elementwise multiplication
- transpose behavior
- shape changes caused by indexing
- MATLAB logical indexing
- default floating-point behavior
- `NaN` handling
- boundary indexing
- interpolation and extrapolation behavior
- solver tolerances and stopping criteria
- ordering of operations
- hidden state, globals, persistent variables, or script workspace variables

Do not simplify these behaviors until equivalence has been demonstrated.

### 3. Make units explicit

Scientific quantities must have clear units.

When adding or modifying model variables:

- Preserve the units used by the original model unless deliberately changing them.
- Document unit conversions near the relevant calculation.
- Avoid unexplained numerical conversion factors.
- Prefer descriptive variable names over ambiguous abbreviations in new Python code.
- Where practical, include units in docstrings or comments for public model variables and parameters.

If units are uncertain, investigate rather than infer silently.

### 4. Preserve numerical reproducibility

Scientific code should be reproducible.

When relevant:

- Use deterministic random seeds in tests.
- Record solver tolerances.
- Avoid unnecessary changes to operation ordering during the parity phase.
- Do not replace an algorithm with a mathematically equivalent one without checking whether numerical output changes materially.
- Compare intermediate states, not only final outputs, when debugging parity failures.

---

## MATLAB reference implementation

The original MATLAB implementation lives under:

`legacy_matlab/`

Treat this directory as a reference implementation.

Unless a task explicitly requires changing the legacy code:

- Do not refactor the MATLAB implementation.
- Do not rename legacy files.
- Do not alter legacy calculations to make Python tests pass.
- Do not overwrite MATLAB-generated reference outputs.

If MATLAB needs to be instrumented to generate additional reference data, keep such changes minimal and clearly documented.

---

## Python implementation

Python source code lives under:

`src/`

Follow the existing repository structure and conventions. Do not create new architectural layers or abstractions unless they solve a concrete problem.

General expectations:

- Prefer small, composable functions.
- Avoid global mutable state.
- Make dependencies explicit.
- Add type annotations to new public functions where practical.
- Use NumPy idiomatically, but not at the expense of numerical parity.
- Keep scientific calculations separate from plotting, file I/O, and presentation logic where practical.
- Do not optimize performance before correctness and parity are established.

When porting a MATLAB routine, preserve a clear correspondence between the legacy routine and its Python replacement until the port is mature.

---

## Testing and validation

Tests live under:

`tests/`

Expected categories may include:

- `tests/unit/` — local behavior of individual functions
- `tests/integration/` — behavior across model components
- `tests/parity/` — comparison with MATLAB reference behavior
- `tests/fixtures/` — fixed inputs and reference outputs

### MATLAB parity tests

For ported scientific calculations, prefer tests against MATLAB-generated reference values.

A useful parity test should record, when available:

- model input
- relevant parameters
- intermediate state
- expected output
- tolerance used for comparison
- provenance of the MATLAB reference result

Do not loosen tolerances simply to make a failing test pass.

If a tolerance must change, explain why the new tolerance is scientifically and numerically justified.

### Debugging parity failures

When Python and MATLAB disagree:

1. Reproduce the failure with the smallest practical case.
2. Compare shapes, indexing, units, and boundary values.
3. Compare intermediate calculations.
4. Identify the earliest point of divergence.
5. Determine whether the difference is caused by:
   - translation error
   - indexing
   - shape / broadcasting behavior
   - floating-point ordering
   - solver behavior
   - MATLAB-specific semantics
   - undocumented legacy behavior
   - an actual discrepancy in the original model
6. Fix the cause rather than compensating for the final output.

Never alter a reference fixture merely because the Python result differs.

---

## Scientific traceability

Maintain the mapping between the scientific description, legacy implementation, Python implementation, and validation.

The primary traceability document is:

`docs/porting/traceability.md`

When porting or substantially modifying scientific functionality, update the traceability documentation where appropriate.

Useful mappings include:

- scientific concept
- paper / supplement section or equation
- MATLAB file and function
- Python module and function
- parity or validation test
- known discrepancy or implementation note

Important scientific reasoning should not live only in chat history, issue comments, or commit messages.

---

## Documentation as persistent project knowledge

This repository should become progressively easier for future humans and agents to understand.

If investigation reveals important information that is not already documented, update the appropriate file under `docs/`.

Examples include:

- meaning of poorly named MATLAB variables
- units
- hidden assumptions
- execution order
- important global state
- boundary conditions
- solver behavior
- unexpected indexing conventions
- undocumented parameter defaults
- differences between the publication and MATLAB implementation
- behavior that initially appears to be a bug but is required for parity
- decisions made during the Python port

Do not leave important discoveries only in the current conversation.

Do not add speculative documentation as fact. Clearly label unresolved interpretations or open questions.

---

## Expected task workflow

For nontrivial porting or scientific tasks, follow this sequence unless the task explicitly requires something else.

### Before editing

1. Read the relevant documentation.
2. Inspect the MATLAB implementation and its dependencies.
3. Inspect the current Python implementation and tests.
4. Identify relevant scientific equations, variables, units, and assumptions.
5. Note any ambiguity or discrepancy that may affect the implementation.

### While editing

1. Make the smallest coherent change.
2. Preserve legacy behavior during the parity phase.
3. Add or update tests.
4. Avoid unrelated refactors.
5. Document important newly discovered behavior.

### After editing

1. Run the relevant tests.
2. Run parity comparisons where available.
3. Review numerical tolerances.
4. Update traceability or other project documentation if needed.
5. Summarize:
   - what changed
   - what was validated
   - remaining discrepancies or uncertainties

---

## Scope discipline

Keep changes narrowly aligned with the requested task.

Do not:

- refactor unrelated modules
- rename scientific concepts merely for stylistic preference
- change numerical algorithms without justification
- replace working scientific calculations with more sophisticated methods during the parity phase
- remove apparently unused MATLAB behavior until its role has been investigated
- introduce new dependencies unless they provide clear value

If a larger redesign would be beneficial, describe it separately rather than bundling it into a behavioral port.

---

## Scientific changes after parity

Once a component has demonstrated adequate MATLAB parity, deliberate scientific changes may be made.

Such changes should be clearly distinguished from the legacy port.

For a scientific modification:

1. State what behavior is changing.
2. Identify the scientific rationale.
3. Identify which equations, assumptions, or parameterizations are affected.
4. Add tests for the new intended behavior.
5. Preserve the ability to compare against the legacy implementation when useful.
6. Update the relevant scientific and traceability documentation.

Do not describe a scientifically changed implementation as a faithful port of the MATLAB model.

---

## Code review priorities

When reviewing changes, prioritize:

1. Scientific correctness
2. Agreement with intended legacy behavior
3. Numerical correctness and stability
4. Units and dimensional consistency
5. Test coverage and reproducibility
6. Traceability to the paper and MATLAB implementation
7. Clarity and maintainability
8. Performance

During the parity phase, elegance and optimization are secondary to correctness and traceability.

---

## Commands

Use the Python 3.11 Conda environment and repository configuration:

```bash
conda run -n mioplio python -m pytest -q
conda run -n mioplio python -m pytest -q tests/parity
conda run -n mioplio ruff check src tests scripts
conda run -n mioplio mypy src scripts
```

Regenerate a local MATLAB reference run with:

```matlab
run('tests/run_swim_allan_hills_test.m')
```

Do not assume an unavailable tool or dependency exists; inspect `pyproject.toml`, environment files, CI configuration, and repository documentation first.

---



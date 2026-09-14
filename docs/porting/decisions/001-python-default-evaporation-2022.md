# Python default evaporation version

## Decision

The Python forward model uses the corrected `evaporation_2022.m` behavior by
default. The frozen `evaporation_2021.m` behavior remains available through
`evaporation_version="2021"` for exact legacy comparisons.

## Scientific and numerical changes

Relative to the frozen 2021 path, the 2022 MATLAB routine:

1. calls the source-condition routine with the physically correct conversion
   `T_K = T_C + 273.15`;
2. changes the Hellmann-Harvey diffusivity exponent from 0.302 to 0.27; and
3. calculates HDO diffusive fractionation directly as
   `(1 / D_r_HDO) ** 0.27` rather than deriving it from oxygen-18.

The Python implementation applies these changes together. It does not describe
the 2022 default as matching the frozen 2021 state-space fixture.

## Validation

A direct MATLAB R2026a run of `evaporation_2022.m` at a 10 °C source
temperature, local closure, Southern Hemisphere annual NCEP climatology was
recorded in `tests/fixtures/matlab/evaporation_2022_Tsource_10.json`. The Python
default matches all recorded source outputs. Existing forward trajectory and
state-space parity tests explicitly select `evaporation_version="2021"` and
continue to validate the frozen legacy baseline.

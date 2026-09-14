# SWIM model overview

SWIM maps a moisture-source air temperature and a condensation temperature to
the isotope composition of precipitation. The frozen baseline uses annual
Southern Hemisphere NCEP source climatology, 2021 local evaporation closure, a
pseudo-adiabatic transport path, and the tuned linear ice-supersaturation
coefficient `b = 0.00525`. The Python scientific default deliberately advances
the source calculation to corrected `evaporation_2022.m` behavior; the frozen
2021 behavior remains selectable for parity.

One forward trajectory has four stages:

1. Source-air temperature is mapped to sea-surface temperature and relative
   humidity using the legacy climatology splines.
2. Ocean evaporation determines the initial vapor isotope ratios.
3. A mixed-phase pseudo-adiabat determines pressure, saturation, and remaining
   vapor along a descending temperature grid.
4. Equilibrium and kinetic fractionation drive stepwise Rayleigh distillation,
   producing precipitation isotope ratios and excess diagnostics.

The state-space model repeats this trajectory across source and condensation
temperature grids and retains each endpoint. The preferred inverse model uses
the modeled log-delta-18O and logarithmic deuterium-excess coordinates to
recover source and condensation temperatures with Sibson natural-neighbor
interpolation. Condensation temperature is converted to surface temperature by
the active legacy Antarctic relationship.

The port preserves baseline MATLAB behavior, including known numerical and
scientific quirks. These are catalogued in `docs/porting/traceability.md` and
must be considered before the model is modified for new scientific work.

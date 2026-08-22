# SWIM scientific assumptions

## Purpose

This document records the scientific assumptions underlying the Simple Water Isotope Model (SWIM), with emphasis on assumptions relevant to the MATLAB-to-Python port.

It is intentionally separate from:

```text
docs/science/equations.md
```

which describes the mathematical formulation, and:

```text
docs/legacy/known_quirks.md
```

which should describe implementation-specific oddities or likely bugs.

An assumption belongs here when it is part of the model's conceptual or scientific simplification, even if it is physically imperfect.

A behavior that exists only because of MATLAB implementation details belongs elsewhere.

Primary scientific source:

Markle, B. R., & Steig, E. J. (2022). Improving temperature reconstructions from ice-core water-isotope records. *Climate of the Past*, 18, 1321–1368. https://doi.org/10.5194/cp-18-1321-2022

Most assumptions below are described in Appendix A.

---

## 1. Temperature is the model's primary transport coordinate

SWIM models isotope evolution along an atmospheric **temperature pathway**, not an explicit geographic trajectory or time-resolved parcel path.

A parcel begins at a source-region air temperature \(T_0\) and is cooled to a final condensation temperature \(T_c\).

The model therefore assumes that the isotope-relevant evolution of the parcel can be represented primarily as a function of temperature.

Implications:

- geographic distance is not a state variable;
- elapsed transport time is not a state variable;
- latitude and longitude are not required for an individual trajectory;
- multiple real atmospheric histories may map to similar temperature trajectories.

This is a deliberate model abstraction, not an implementation limitation.

---

## 2. Moisture transport is represented by idealized parcel trajectories

The base model treats an evaporated air parcel as an idealized parcel that cools along a prescribed thermodynamic pathway.

The central configuration is pseudo-adiabatic.

The parcel is not modeled with a complete atmospheric circulation model.

Processes such as synoptic dynamics, explicit frontal evolution, and three-dimensional transport are not directly resolved.

---

## 3. Base transport is pseudo-adiabatic

The base SWIM trajectory assumes pseudo-adiabatic cooling.

As the parcel cools:

- water condenses;
- latent heat affects the thermodynamic path;
- condensed water is removed from the parcel.

The model does not retain condensate and allow full reversible equilibration with the vapor.

This assumption is fundamental to the modeled relationship between temperature, remaining vapor fraction, and isotope ratios.

---

## 4. Condensate is removed immediately

Once condensation occurs, the resulting condensate is treated as precipitation and is removed from the parcel.

This is the Rayleigh-style removal assumption.

Consequences:

- previously condensed precipitation does not remain available to re-equilibrate with the vapor;
- falling precipitation is not part of the modeled parcel state;
- the vapor reservoir becomes progressively smaller and isotopically depleted during transport.

---

## 5. Re-evaporation of falling precipitation is not represented in the core transport model

The paper explicitly notes that falling precipitation is not re-evaporated within the idealized SWIM transport calculation.

Therefore the forward isotope trajectory does not include isotope exchange or fractionation caused by precipitation re-evaporation below the condensation level.

This assumption should not be changed during the parity port.

---

## 6. Initial source conditions may be inferred from modern climatological relationships

In the base model, source-region:

```text
sea-surface temperature
relative humidity
normalized relative humidity
```

can be inferred from the initial source air temperature \(T_0\) using relationships derived from modern reanalysis.

The paper uses modern climatological correlations rather than requiring source latitude or explicit atmospheric trajectories.

This assumes that these relationships are useful approximations for the climate states being reconstructed.

The model tests sensitivity to this assumption, but it remains part of the base framework.

---

## 7. Modern climatological source relationships are applied to paleoclimate reconstruction

When SWIM is used for paleoclimate reconstruction, modern relationships between source air temperature, SST, and RH are used as a baseline parameterization.

This does not imply that these relationships are invariant in all climates.

Rather, their possible changes contribute to model uncertainty.

The distinction between:

```text
base model assumption
```

and:

```text
uncertainty / sensitivity experiment
```

should remain explicit in the Python implementation.

---

## 8. Source-region temperature represents a moisture-weighted source condition

The reconstructed \(T_0\) should not be interpreted as the temperature of one unique geographic source point.

Moisture reaching an Antarctic site originates from a broad source distribution.

The reconstructed source temperature is therefore interpreted as a moisture-weighted effective evaporation/source temperature.

This distinction is scientifically important and should be reflected in API and documentation naming where practical.

---

## 9. A single modeled trajectory represents one end-member transport history

The paper interprets real precipitation as arising from a distribution of source and condensation conditions.

A single trajectory is therefore not assumed to represent every water molecule arriving at a site.

The full SWIM state space represents a family of possible trajectories.

The reconstruction framework uses this family to infer effective source and condensation temperatures from observed isotope values.

---

## 10. Base source evaporation uses a closure assumption

The isotopic composition of boundary-layer vapor cannot be determined from local evaporation alone without an assumption about how the boundary layer is supplied and mixed.

SWIM therefore includes closure assumptions.

### Local closure

The base configuration assumes local closure, in which boundary-layer vapor is linked directly to local evaporate.

This is an idealized end-member.

### Global closure

A global closure option represents a strongly mixed atmospheric end-member.

The paper treats closure choice as a source of uncertainty.

Neither closure option should be interpreted as a complete physical model of marine boundary-layer mixing.

---

## 11. Ocean isotopic composition is prescribed rather than dynamically simulated

The model requires an initial seawater isotopic composition.

The checked-in MATLAB evaporation implementation has a default Southern Ocean \(\delta^{18}O\) value and derives corresponding \(\delta D\) using an observational relationship.

The ocean isotope state is not dynamically evolved by SWIM.

When applying the model to past climates, seawater isotope changes are handled through prescribed corrections or sensitivity experiments.

---

## 12. Evaporation kinetic fractionation is parameterized

Kinetic isotope effects during ocean evaporation are not calculated from a fully resolved turbulent boundary layer.

Instead, the model uses effective diffusive fractionation relationships.

These relationships summarize unresolved effects of:

- molecular diffusion;
- turbulence;
- relative humidity;
- isotopologue diffusivities.

The effective parameters are constrained using laboratory, theoretical, and observational results.

---

## 13. Relative humidity relevant to evaporation is normalized to SST

The evaporation scheme uses normalized relative humidity \(RH_n\), not simply the reported atmospheric RH.

\(RH_n\) accounts for the difference between saturation vapor pressure at the air temperature and at the sea-surface temperature.

This assumes that the air-sea humidity disequilibrium relevant to isotope fractionation can be represented by this normalized quantity.

---

## 14. Cloud condensate can contain both liquid water and ice

The model does not use a single abrupt liquid-to-ice transition at 0 °C.

Instead, the cloud is treated as mixed phase over a temperature range.

Liquid and ice fractions vary smoothly with temperature using an empirical parameterization based on satellite observations.

This affects both:

- the thermodynamic trajectory;
- the effective isotope fractionation factor.

---

## 15. Cloud liquid/ice fraction is determined primarily by temperature

The mixed-phase parameterization represents liquid/ice partitioning as a function of temperature.

Other controls on cloud phase, such as:

- aerosol concentration;
- vertical velocity;
- cloud age;
- local microphysics;

are not explicit state variables in the base SWIM formulation.

Their effects are absorbed into the empirical parameterization and model uncertainty.

---

## 16. Supersaturation over ice is parameterized

At cold temperatures, vapor may be supersaturated with respect to ice.

SWIM does not attempt to explicitly simulate all microphysical processes that create this supersaturation.

Instead, supersaturation is parameterized as a function of temperature.

The paper's preferred simple form is:

\[
S_i = a - bT
\]

with the base local-closure tuning approximately:

```text
a = 1
b = 0.00525 deg C^-1
```

The MATLAB framework also supports a quadratic term.

The supersaturation relationship is a calibrated model parameterization, not a fundamental physical law.

---

## 17. The same supersaturation framework should control moisture removal and kinetic fractionation

A central improvement discussed in the paper is the use of a consistent view of supersaturation across:

- thermodynamic moisture removal;
- isotope kinetic fractionation.

The model assumes that these should be physically linked.

Using inconsistent supersaturation relationships in different parts of the model can create unphysical isotope-temperature behavior.

This coupling should be treated as a core scientific design constraint.

---

## 18. Equilibrium isotope fractionation depends on temperature and phase

Equilibrium fractionation factors are prescribed from empirical or experimentally constrained relationships.

The model assumes those published relationships adequately describe equilibrium partitioning between:

- vapor and liquid;
- vapor and ice;

over the relevant temperature range.

These parameterizations may be updated scientifically later, but should be preserved during parity work.

---

## 19. Kinetic fractionation during cold condensation is represented through effective diffusivity ratios

Kinetic fractionation during ice formation is parameterized from:

- supersaturation;
- equilibrium fractionation;
- diffusivity ratios.

The model does not explicitly simulate molecular collision dynamics or ice crystal growth.

The diffusivity parameterization is therefore an effective representation of unresolved kinetic processes.

---

## 20. Rayleigh distillation is applied incrementally

The model assumes that the change in vapor isotope ratio at each temperature step can be represented by the Rayleigh relation:

\[
d\ln R = (\alpha - 1)d\ln f.
\]

Because the fractionation factors change with temperature, the calculation is performed incrementally.

The numerical solution therefore approximates a continuously varying trajectory using discrete steps.

---

## 21. The base model uses a fixed small temperature step

The publication and checked-in wrapper use an internal temperature increment of approximately:

```text
0.1 deg C
```

This is a numerical discretization assumption.

The eventual Python model may investigate convergence with smaller or larger steps, but the legacy step should be preserved during parity testing.

---

## 22. Euler-style numerical integration is considered adequate for the legacy model

The paper describes Euler numerics along the temperature trajectory.

The MATLAB implementation uses explicit sequential stepping.

This is both a numerical and practical modeling assumption.

The parity port should reproduce the legacy integration before replacing it with a higher-order method.

A higher-order method could be a later numerical improvement, but it would no longer be strict legacy parity.

---

## 23. Condensation temperature is not surface temperature

The reconstructed \(T_c\) represents an effective condensation-weighted atmospheric temperature.

It is not the physical surface temperature at the ice-core site.

For Antarctic applications, the model uses an empirical relationship between condensation and surface temperature.

The paper's base relationship is approximately:

\[
T_c = 0.69T_s - 8.2.
\]

The relationship itself carries uncertainty.

---

## 24. The condensation-to-surface relationship is an empirical parameterization

SWIM does not derive the local surface temperature directly from first principles of boundary-layer inversion physics.

Instead, the paper uses a relationship informed by observations and reanalysis/model results over Antarctica.

This conversion should remain conceptually separate from the core isotope distillation model.

---

## 25. The inverse problem is solved through a modeled isotope state space

The reconstruction method assumes that the forward SWIM state space contains a sufficiently informative and locally invertible mapping between isotope coordinates and:

```text
source temperature
condensation temperature
```

The model does not derive a closed-form inverse.

Instead, it interpolates the forward-model surface.

This introduces assumptions about:

- state-space coverage;
- interpolation behavior;
- uniqueness / local non-uniqueness;
- treatment of observations outside the modeled domain.

---

## 26. Logarithmic deuterium excess is preferred for reconstruction

The paper argues that \(d_{ln}\) is a more faithful coordinate for source conditions than the traditional linear \(d_{xs}\), particularly after strong distillation.

The preferred nonlinear inversion therefore uses:

```text
log-transformed delta18O
dln
```

rather than relying exclusively on the historical linear deuterium-excess definition.

This is a scientific design choice, not merely a plotting convention.

---

## 27. Atmospheric mixing is not explicitly resolved in the base forward path

The simplest SWIM framework treats individual pseudo-adiabatic trajectories as independent.

The paper separately investigates mixing among parcels and finds that mixing can broaden isotope distributions and introduce relatively modest changes in mean reconstructed relationships under the tested scenarios.

Thus, explicit stochastic parcel mixing is treated primarily as a sensitivity or uncertainty experiment rather than part of the base trajectory.

The exact MATLAB implementation of the paper's mixing experiments still needs to be mapped before any Python port of that functionality.

---

## 28. The model is intentionally simple and does not resolve all hydrological processes

SWIM is not an isotope-enabled general circulation model.

Important processes that are absent or highly parameterized include, among others:

- explicit atmospheric circulation;
- convection as a resolved dynamical process;
- detailed cloud microphysics;
- precipitation re-evaporation;
- post-depositional snow processes;
- spatially explicit source tagging within each trajectory;
- event-scale intermittency.

These omissions are part of the model hierarchy that makes SWIM interpretable and computationally lightweight.

---

## 29. Parameter uncertainty and structural uncertainty are distinct

The paper explores uncertainty both in parameter values and in model choices.

Examples include:

```text
supersaturation tuning
closure assumption
source RH parameterization
reanalysis product
evaporation fractionation parameters
condensation-to-surface conversion
mixing
```

The Python implementation should avoid representing all uncertainty as one generic parameter error.

Where practical, preserve the distinction between:

- uncertain numerical parameter;
- uncertain empirical relationship;
- alternative structural assumption.

---

## 30. Absolute reconstruction uncertainty is larger than relative variability uncertainty

The paper notes that many model perturbations shift reconstructed temperature surfaces while preserving similar patterns of variability.

Therefore uncertainty in absolute reconstructed \(T_0\) and \(T_c\) can be larger than uncertainty in their relative temporal changes.

This is an interpretation assumption relevant to later uncertainty workflows, not to the basic forward parity implementation.

---

## 31. Base assumptions versus legacy implementation quirks

The following should **not** automatically be treated as scientific assumptions:

- a suspicious Celsius/Kelvin conversion in a dated MATLAB file;
- one dated routine calling another dated routine unexpectedly;
- inconsistent function/file naming;
- stale absolute filesystem paths;
- branch-specific use of different dated evaporation routines;
- questionable comments in experimental code;
- MATLAB array orientation;
- MATLAB `griddata` implementation details.

Those belong in:

```text
docs/legacy/known_quirks.md
```

and/or:

```text
docs/porting/traceability.md
```

until their provenance is resolved.

---

## 32. Assumptions that must remain configurable in Python

The Python design should eventually expose, rather than bury, scientifically meaningful choices such as:

```text
closure assumption
reanalysis/source-condition dataset
hemisphere
season
source seawater isotope composition
supersaturation coefficients
transport pathway
fractionation parameterization
cloud phase parameterization
temperature step
```

During parity, defaults should match the frozen MATLAB baseline.

After parity, alternative scientifically justified parameterizations may be added explicitly rather than by editing internal constants.

---

## 33. Change policy

When changing one of these assumptions:

1. identify the assumption being changed;
2. state whether the change is scientific, numerical, or software-only;
3. identify which equations or model components are affected;
4. add a decision record under `docs/porting/decisions/`;
5. add or update tests;
6. preserve a way to reproduce the legacy result when useful.

A deliberate change in a scientific assumption should never be presented as a pure refactor.

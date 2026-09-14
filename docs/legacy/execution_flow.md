# Frozen MATLAB execution flow

The frozen `matlab-port-baseline-v1` follows this active MATLAB call chain:

```text
simple_water_isotope_model_2020
├── evaporation_2021
│   ├── T_RH_RHn_2020
│   └── d18Osw_to_dDsw
└── distillation_2020
    ├── fraction_il_brm_H10 (method "adj")
    ├── mixed_phased_supersaturation
    └── pseudo_adiabat_function
```

For each source/condensation grid pair, the wrapper first calculates the source
environment and initial vapor composition. It then constructs the inclusive
descending temperature path at 0.1 °C spacing. Mixed-phase saturation is
integrated first; pressure and remaining vapor are integrated in a second Euler
loop. Distillation applies effective fractionation step by step and stores the
last precipitation value in the state-space arrays.

`Tsite_Tsource_reconstruction_quick` transforms observations into
log-delta-18O and logarithmic deuterium-excess coordinates. MATLAB `griddata`
with the `natural` method returns condensation and source temperatures;
ordinary linear `griddata` returns saturated mixing ratio. `Ts_to_Tc_2020`
then converts condensation to surface temperature.

The Python equivalents are `swim.model.forward_trajectory`,
`swim.model.forward_state_space`, and
`swim.reconstruction.reconstruct_temperatures`. Detailed variable-level
mapping and preserved discrepancies are in `docs/porting/traceability.md`.

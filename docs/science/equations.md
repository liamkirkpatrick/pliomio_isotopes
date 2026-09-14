# SWIM equations and scientific formulation

## Purpose and scope

This document is a working scientific specification for the Simple Water Isotope Model (SWIM) described by:

> Markle, B. R., & Steig, E. J. (2022). Improving temperature reconstructions from ice-core water-isotope records. *Climate of the Past*, 18, 1321–1368. https://doi.org/10.5194/cp-18-1321-2022

It is intended to support a faithful MATLAB-to-Python port. It distinguishes:

- **Paper formulation** — equations and assumptions stated in the publication.
- **Checked-in MATLAB behavior** — equations and parameterizations observed in the GitHub repository.
- **Open questions** — places where the publication, dated MATLAB variants, or active call chain require further verification.

Repository snapshot inspected for this draft:

- `bradley-markle/simple_water_isotope_model`
- commit: `4db23b79aab8f2154254d9111f77e70a5ce8427d`

The current repository contains multiple dated implementations. Do not assume that the numerically latest filename is the implementation used for every published result.

---

## 1. Conceptual model

SWIM represents the isotopic evolution of water vapor that:

1. evaporates from the ocean under specified or climatologically inferred source conditions;
2. is transported along a prescribed thermodynamic cooling pathway from an initial source air temperature \(T_0\) to a final condensation temperature \(T_c\);
3. loses moisture through condensation and precipitation;
4. experiences equilibrium and kinetic isotope fractionation;
5. transitions between liquid and ice condensate as temperature falls; and
6. may become supersaturated with respect to ice at low temperature.

The model is organized in **temperature space**, not geographic space or explicit time. A single trajectory is a sequence

\[
T = T_0,\ T_0-\Delta T,\ \ldots,\ T_c.
\]

The paper's state-space calculations use \(\Delta T=0.1\,^\circ\mathrm{C}\) and a large ensemble of \((T_0,T_c)\) trajectories. The checked-in top-level MATLAB wrapper also hard-codes `dT = 0.1`.

The model does **not** attempt to model post-depositional processes.

---

## 2. Primary variables and notation

| Symbol | Meaning | Typical units / scale |
|---|---|---|
| \(T_0\) | initial source-region surface air / evaporation air temperature | °C |
| \(T_c\) | final condensation temperature | °C |
| \(T_s\) | surface air temperature at deposition site | °C |
| SST\(_0\) | source-region sea-surface temperature | °C |
| RH\(_0\) | source-region relative humidity | fraction or % depending on routine |
| RH\(_n\) | normalized source-region relative humidity | fraction |
| \(P\) | air pressure along pathway | kPa in MATLAB thermodynamic routines |
| \(e_s\) | saturation vapor pressure | kPa |
| \(r_s\) | saturated water-vapor mixing ratio | kg/kg in equations; MATLAB comments sometimes call it g/kg |
| \(f\) | fraction of initial vapor remaining | dimensionless |
| \(R\) | heavy/light isotope number ratio | dimensionless |
| \(\delta^{18}O,\delta D,\delta^{17}O\) | isotope delta values relative to VSMOW | ‰ |
| \(\alpha_{\rm eq}\) | equilibrium fractionation factor | dimensionless |
| \(\alpha_{\rm diff}\) | diffusive fractionation factor | dimensionless |
| \(\alpha_k\) | kinetic fractionation factor during condensation | dimensionless |
| \(\alpha_{\rm tot}\) | total fractionation factor | dimensionless |
| \(F_{\rm liq},F_{\rm ice}\) | liquid/ice condensate fractions | dimensionless |
| \(S_i\) | supersaturation with respect to ice | dimensionless |
| \(d_{xs}\) | linear deuterium excess | ‰ |
| \(d_{ln}\) | logarithmic deuterium excess | dimensionless in paper; per-mil-like transformed scale in MATLAB |

### Important temperature convention

Most public SWIM inputs are in degrees Celsius. Thermodynamic equations requiring absolute temperature use

\[
T_K = T + 273.15.
\]

This conversion is a parity-critical detail because one checked-in historical routine contains a likely sign error; see **Open questions and discrepancies**.

---

## 3. Isotope notation

### 3.1 Delta notation — paper Eq. (1)

For isotope ratio \(R_x\) and standard ratio \(R_{\rm std}\),

\[
\delta_x =
\frac{R_x-R_{\rm std}}{R_{\rm std}}.
\]

The paper reports \(\delta\) in per mil (‰). In code, conversion between ratio and per-mil delta is implemented as

\[
\delta_{x,\permil}
=
1000\left(\frac{R_x}{R_{\rm std}}-1\right),
\]

and

\[
R_x =
\left(1+\frac{\delta_{x,\permil}}{1000}\right)R_{\rm std}.
\]

Checked-in VSMOW ratios used in the core MATLAB routines are:

\[
R^{18}_{\rm VSMOW}=0.00200520,
\]

\[
R^{D}_{\rm VSMOW}=0.00015576,
\]

\[
R^{17}_{\rm VSMOW}=0.0003799.
\]

These constants should be preserved exactly during parity work.

### 3.2 Linear deuterium excess

\[
d_{xs} = \delta D - 8\delta^{18}O.
\]

The paper emphasizes that this historical linear definition develops substantial nonlinear bias during distillation, especially at cold Antarctic conditions.

### 3.3 Logarithmic isotope transform and \(d_{ln}\) — paper Eq. (4)

The paper defines

\[
\delta'_x = \ln(1+\delta_x),
\]

where \(\delta_x\) is expressed as a unitless ratio rather than the numerical per-mil value.

The logarithmic deuterium excess is

\[
d_{ln}
=
\delta'_D -
\left[
A(\delta'_{18})^2+B\delta'_{18}
\right],
\]

with

\[
A=-28.5,\qquad B=8.47.
\]

#### MATLAB scaling convention

The checked-in MATLAB routines use a transformed quantity scaled by 1000:

\[
\widetilde{\delta}'_x
=
1000\ln\left(1+\frac{\delta_{x,\permil}}{1000}\right)
=
1000\ln\left(\frac{R_x}{R_{\rm std}}\right).
\]

On this scale the same curve is implemented as

\[
d_{ln,\mathrm{MATLAB}}
=
\widetilde{\delta}'_D
-
\left[
-0.0285(\widetilde{\delta}'_{18})^2
+
8.47\widetilde{\delta}'_{18}
\right].
\]

This scaling distinction is easy to miss and must be tested explicitly in the Python port.

### 3.4 \(^{17}O\) excess

The code uses

\[
^{17}O_{xs}
=
10^6
\left[
\ln\left(1+\frac{\delta^{17}O}{1000}\right)
-
0.528
\ln\left(1+\frac{\delta^{18}O}{1000}\right)
\right].
\]

This is generally reported in per meg.

---

## 4. Source-region environmental conditions

### 4.1 Climatological mapping from \(T_0\)

When SST\(_0\) and RH\(_0\) are not supplied explicitly, SWIM obtains them from climatological relationships with source-region air temperature \(T_0\).

The paper describes three fit strategies and uses a smoothing cubic spline as the base method. The checked-in model calls `T_RH_RHn_2020(..., 'spline', ...)` from `evaporation_2021.m`.

The climatological fit data are stored under:

- `data/ncep_data/`
- `data/era_data/`

NCEP only supports annual means in the checked-in top-level wrapper; selecting NCEP forces `season = 'annual'`.

### 4.2 Normalized relative humidity — paper Eq. (A1)

\[
RH_n
=
RH
\frac{e_s(T_a)}{e_s(SST)}.
\]

This is the humidity relevant to kinetic fractionation during evaporation.

---

## 5. Saturation vapor pressure and mixing ratio

### 5.1 Saturated mixing ratio — paper Eq. (A2)

\[
r_s
=
\frac{R_d}{R_{wv}}
\frac{e_s}{P-e_s}.
\]

The MATLAB implementation uses

\[
\eta = \frac{R_d}{R_{wv}} = 0.622,
\]

so

\[
r_s
=
0.622\frac{e_s}{P-e_s}.
\]

### 5.2 Saturation vapor pressure over liquid water

The checked-in pseudo-adiabatic and evaporation routines use the Murphy and Koop formulation:

\[
e_{s,l}
=
10^{-3}\exp
\left[
54.842763
-\frac{6763.22}{T_K}
-4.21\ln T_K
+0.000367T_K
+\tanh\left(0.0415(T_K-218.8)\right)
\left(
53.878
-\frac{1331.22}{T_K}
-9.44523\ln T_K
+0.014025T_K
\right)
\right],
\]

with \(e_{s,l}\) in kPa.

### 5.3 Saturation vapor pressure over ice

\[
e_{s,i}
=
10^{-3}
\exp
\left[
9.550426
-\frac{5723.265}{T_K}
+3.53068\ln T_K
-0.00728332T_K
\right],
\]

again in kPa.

---

## 6. Cloud liquid/ice partitioning

The paper uses temperature-dependent cloud ice and liquid fractions based on Hu et al. (2010).

The active path in `distillation_2020.m` calls:

```text
fraction_il_brm_H10(T, 'adj')
```

For the active `"adj"` option,

\[
p(T)=
5.37
+0.4025T
+0.0847T^2
+0.007182T^3
+2.39\times10^{-4}T^4
+2.87\times10^{-6}T^5,
\]

\[
F_{\rm liq}
=
\frac{1}{1+\exp[-p(T)]},
\]

\[
F_{\rm ice}=1-F_{\rm liq}.
\]

This numerical polynomial is an **implementation detail**; the paper describes the satellite-based parameterization conceptually rather than listing these coefficients.

---

## 7. Pseudo-adiabatic transport

The paper describes pseudo-adiabatic transport following Bakhshaii and Stull (2013), with immediate removal of condensed moisture and mixed liquid/ice condensate.

The checked-in routine `pseudo_adiabat_function.m` integrates pressure as a function of a prescribed temperature grid.

### 7.1 Effective latent heat and heat-capacity terms

The code forms temperature-dependent mixed-phase effective quantities from \(F_{\rm liq}\) and \(F_{\rm ice}\).

For example,

\[
L_{\rm eff}
=
F_{\rm liq}L_v + F_{\rm ice}L_i.
\]

The implementation also forms an effective heat-capacity ratio from liquid-vapor and ice heat capacities.

### 7.2 Pseudo-adiabatic pressure gradient

To avoid confusing this coefficient with the supersaturation slope \(b\), define

\[
B_p
=
\frac{1+r_s/\eta}{1+r_s/c_{\rm eff}}.
\]

The MATLAB routine uses

\[
\frac{dP}{dT}
=
\frac{P}{B_p}
\frac{
C_{pd}
+
\frac{L_{\rm eff}^2 r_s\eta B_p}{R_dT_K^2}
}{
R_dT_K + L_{\rm eff}r_s
}.
\]

Euler stepping on the fixed temperature grid is then

\[
P_{i+1}
=
P_i
+
\left(\frac{dP}{dT}\right)_i
(T_{K,i+1}-T_{K,i}).
\]

The model initializes

\[
P_0=101.325\ {\rm kPa}.
\]

### 7.3 Fraction of water remaining

The paper defines, Eq. (A10),

\[
f
=
\frac{q}{q_0}
=
\frac{r_s}{r_{s,0}}.
\]

The checked-in pseudo-adiabatic routine returns exactly

\[
f_i=\frac{r_{s,i}}{r_{s,1}}.
\]

This quantity drives Rayleigh distillation.

---

## 8. Equilibrium fractionation

### 8.1 Definition — paper Eq. (A3)

For example, for \(^{18}O\) between liquid and vapor,

\[
^{18}\alpha_{l-v}
=
\frac{^{18}R_l}{^{18}R_v}.
\]

### 8.2 Active MATLAB equilibrium factors

The checked-in evaporation and distillation routines contain several historical alternatives, mostly commented out. The active formulas relevant to the example-model path include the following.

#### Deuterium, liquid-vapor

\[
\alpha^D_{eq,l}
=
\exp
\left[
\frac{
52.612
-76.248\times10^3/T_K
+24.844\times10^6/T_K^2
}{1000}
\right].
\]

#### Deuterium, ice-vapor

Using Lamb et al. (2017),

\[
\alpha^D_{eq,i}
=
\exp
\left(
-0.0559+\frac{13525}{T_K^2}
\right).
\]

#### \(^{18}O\), liquid-vapor

\[
\alpha^{18}_{eq,l}
=
\exp
\left[
\frac{
-2.0667
-0.4156\times10^3/T_K
+1.137\times10^6/T_K^2
}{1000}
\right].
\]

#### \(^{18}O\), ice-vapor

\[
\alpha^{18}_{eq,i}
=
\exp
\left[
\frac{
-28.224
+11.839\times10^3/T_K
}{1000}
\right].
\]

#### \(^{17}O\)

`evaporation_2021.m` uses

\[
\alpha^{17}_{eq,l}
=
(\alpha^{18}_{eq,l})^{0.529},
\qquad
\alpha^{17}_{eq,i}
=
(\alpha^{18}_{eq,i})^{0.529}.
\]

`distillation_2020.m`, however, uses

\[
\alpha^{17}_{eq,l}
=
(\alpha^{18}_{eq,l})^{0.529},
\qquad
\alpha^{17}_{eq,i}
=
(\alpha^{18}_{eq,i})^{0.531}.
\]

That difference should be preserved for parity until its provenance is resolved.

---

## 9. Kinetic fractionation during evaporation

### 9.1 Diffusive fractionation — paper Eq. (A4)

\[
\alpha_{\rm diff}
=
\left(\frac{D}{D^*}\right)^n,
\]

where \(D\) and \(D^*\) are light and heavy isotopologue diffusivities and \(n\) represents the turbulent/molecular diffusion regime.

### 9.2 Hydrogen/oxygen relationship — paper Eq. (A5)

\[
\phi_{\rm diff}
=
\frac{
{}^D\alpha_{\rm diff}-1
}{
{}^{18}\alpha_{\rm diff}-1
}.
\]

The paper discusses a preferred effective \(^{18}\alpha_{\rm diff}\approx1.009\) for initial evaporation and several alternative constraints.

The checked-in `evaporation_2021.m` contains both the fixed-\(\alpha\) formulation and a later Hellmann-and-Harvey-based formulation. In the inspected file, the active branch sets a diffusivity exponent \(N=0.302\) and calculates temperature-dependent diffusivity ratios. This is one reason the exact publication baseline should be pinned with reference output rather than inferred from filenames alone.

---

## 10. Ocean evaporation and closure assumptions

### 10.1 General evaporative fractionation — paper Eq. (A6)

The paper writes

\[
\alpha_{\rm evap}
=
\frac{R_o}{R_e}
=
\frac{
\alpha_{\rm eq}\alpha_{\rm diff}(1-RH_n)
}{
1-\alpha_{\rm eq}RH_n(R_v/R_o)
},
\]

where:

- \(R_o\): ocean isotope ratio,
- \(R_e\): net evaporate isotope ratio,
- \(R_v\): boundary-layer vapor isotope ratio.

### 10.2 Local closure — paper Eq. (A7)

Assuming the boundary-layer vapor is supplied only by local evaporate, \(R_v=R_e\):

\[
R_v
=
\frac{
R_o
}{
\alpha_{\rm eq}
\left[
\alpha_{\rm diff}
+
RH_n(1-\alpha_{\rm diff})
\right]
}.
\]

This same expression is active in the `closure == 'local'` branch of `evaporation_2021.m`.

### 10.3 Global closure — paper Eq. (A8)

The paper also defines a globally mixed end-member using a prescribed global \(\alpha_{\rm evap}=R_o/R_e\). Rearrangement gives

\[
R_v
=
R_o
\frac{
1-
\alpha_{\rm eq}\alpha_{\rm diff}(1-RH_n)/\alpha_{\rm evap}
}{
\alpha_{\rm eq}RH_n
}.
\]

The checked-in code uses:

\[
{}^{18}\alpha_{\rm evap}=1.0045,
\qquad
{}^{D}\alpha_{\rm evap}=1.0267
\]

for this branch.

### 10.4 Initial seawater composition

The active `evaporation_2021.m` default is

\[
\delta^{18}O_{\rm ocean}=-0.3\ \permil
\]

when no explicit seawater value is supplied.

It then calls `d18Osw_to_dDsw.m` to infer \(\delta D_{\rm ocean}\) from observed seawater relationships.

---

## 11. Rayleigh distillation during transport

### 11.1 Fundamental differential equation — paper Eq. (A9)

\[
\frac{d\ln R}{d\ln f}
=
\alpha-1.
\]

### 11.2 Total fractionation — paper Eq. (A11)

When kinetic fractionation is active,

\[
\alpha_{\rm tot}
=
\alpha_{\rm eq}\alpha_k,
\]

so

\[
d\ln R
=
(\alpha_{\rm tot}-1)d\ln f.
\]

The MATLAB code performs this explicitly at each temperature step:

\[
\Delta\ln R_i
=
(\alpha_i-1)
\left[
\ln(f_i)-\ln(f_{i-1})
\right],
\]

followed by

\[
R_{v,i}
=
\exp\left[\ln(R_{v,i-1})+\Delta\ln R_i\right].
\]

Precipitation is then calculated as

\[
R_{p,i}=\alpha_iR_{v,i}.
\]

If the change in \(f\) is zero, the MATLAB routine assigns the precipitation ratio `NaN`, because no new precipitation has formed.

---

## 12. Kinetic fractionation during condensation

### 12.1 Paper Eq. (A12)

For vapor supersaturated with respect to ice,

\[
\alpha_k
=
\frac{
S_i
}{
\alpha_{\rm eq}
\left(D/D^*\right)(S_i-1)+1
}.
\]

The checked-in distillation implementation evaluates this independently for D, \(^{18}O\), and \(^{17}O\).

### 12.2 Transport diffusivities

The paper discusses the commonly used fixed transport ratios

\[
D^{16}/D^{18}=1.0285,
\qquad
D^1/D^2=1.0251,
\]

and also the temperature-dependent diffusivity ratios of Hellmann and Harvey (2020).

The active code in both inspected `distillation_2020.m` and `distillation_2022.m` calculates the Hellmann–Harvey ratios:

\[
T_*=\frac{T_K}{100},
\]

\[
D_{r,\mathrm{HDO}}
=
0.98258-\frac{0.02546}{T_*}
+\frac{0.02421}{T_*^{5/2}},
\]

\[
D_{r,17}
=
0.98284+\frac{0.003517}{T_*^{1/2}}
-\frac{0.001996}{T_*^{5/2}},
\]

\[
D_{r,18}
=
0.96671+\frac{0.007406}{T_*^{1/2}}
-\frac{0.004861}{T_*^3},
\]

and then uses

\[
D_O=1/D_{r,18},\qquad
D_D=1/D_{r,\mathrm{HDO}},\qquad
D_{17}=1/D_{r,17}.
\]

Whether the published central state-space files were generated with this exact configuration should be verified from saved outputs / generation scripts.

---

## 13. Mixed-phase effective fractionation — paper Eq. (A13)

\[
\alpha_{\rm eff}
=
\alpha_{\rm tot(l-v)}F_{\rm liq}
+
\alpha_{\rm tot(i-v)}F_{\rm ice}.
\]

The distillation routine implements the same weighted combination separately for D, \(^{18}O\), and \(^{17}O\).

---

## 14. Supersaturation

### 14.1 Paper parameterization — Eq. (A14)

The simple linear parameterization is

\[
S_i=a-bT.
\]

For the paper's tuned local-closure configuration,

\[
a=1,\qquad
b=0.00525\ ^\circ\mathrm{C}^{-1}.
\]

The MATLAB thermodynamic routine generalizes this to

\[
S_i=a-bT-cT^2
\]

and clips values below one:

\[
S_i=\max(S_i,1).
\]

The standard example passes

\[
a=1,\quad b=0.00525,\quad c=0.
\]

### 14.2 Consistency between precipitation and isotope schemes

A scientifically important SWIM design choice is that the supersaturation controlling moisture removal is made consistent with the supersaturation entering kinetic fractionation. The paper emphasizes that inconsistent supersaturation assumptions between these parts of the model can create unphysical curvature in the isotope-temperature relationship.

This coupling must be preserved during the Python port.

---

## 15. Surface versus condensation temperature

The reconstructed \(T_c\) is not the surface temperature. It represents the condensation-weighted atmospheric temperature contributing net accumulation.

For Antarctica, the paper adopts the base relationship

\[
T_c=0.69T_s-8.2^\circ\mathrm{C},
\]

with an uncertainty of approximately \(\pm0.02\) in the slope.

The checked-in `Ts_to_Tc_2020.m` implements exactly:

```text
Tc = 0.69 * Ts - 8.2
```

and its inverse

\[
T_s=\frac{T_c+8.2}{0.69}.
\]

---

## 16. Forward isotope state space

The forward model defines a mapping

\[
(T_0,T_c)
\longrightarrow
\left(
\delta^{18}O,
\delta D,
d_{xs},
d_{ln},
\ldots
\right).
\]

Under the top-level MATLAB wrapper, each source temperature is paired with each site/condensation temperature. The endpoint isotope composition of precipitation is saved into 2-D state-space arrays.

The principal dimensions in the wrapper are:

```text
axis 1 / MATLAB i : source temperature T0
axis 2 / MATLAB j : condensation temperature Tc
axis 3 / MATLAB k : source RH, only when RH is explicitly swept
```

This orientation should be made explicit in Python rather than inferred from array shape.

---

## 17. Nonlinear temperature reconstruction

The publication's main reconstruction approach avoids a fixed linearization of isotope-temperature sensitivities.

Instead, SWIM's forward state space is inverted numerically:

\[
(\delta^{18}O,d_{ln})
\longrightarrow
(T_c,T_0).
\]

The checked-in `Tsite_Tsource_reconstruction_quick.m`:

1. transforms measured \(\delta^{18}O\) and \(\delta D\) to the MATLAB log-isotope scale;
2. calculates \(d_{ln}\);
3. loads a precomputed SWIM state space;
4. removes `NaN` regions;
5. uses MATLAB `griddata(..., 'natural')` over \((\widetilde{\delta}'_{18},d_{ln})\);
6. interpolates both \(T_c\) and \(T_0\).

The paper argues that \((\delta^{18}O,d_{ln})\) is preferable to \((\delta^{18}O,\delta D)\) or \((\delta^{18}O,d_{xs})\) because the \(T_0\) contours are more nearly orthogonal to the \(d_{ln}\) coordinate, reducing propagated reconstruction uncertainty.

---

## 18. Linear reconstruction equations retained for comparison

The traditional approximation in the paper is

\[
\Delta \delta^{18}O
=
\gamma_1\Delta T_{\rm site}
+
\gamma_2\Delta T_{\rm source},
\]

\[
\Delta d_{xs}
=
\beta_1\Delta T_{\rm site}
+
\beta_2\Delta T_{\rm source}.
\]

A central result of the paper is that the \(\beta\) and \(\gamma\) coefficients are functions of mean climate state; treating them as fixed over large temperature excursions can bias reconstructed temperatures, especially source-region temperature.

The Python implementation should therefore treat the nonlinear state-space inversion as first-class behavior, while preserving linear calculations only where needed for reproduction or comparison.

---

## 19. Base / example parameterization observed in the repository

| Quantity | Value / behavior | Source |
|---|---|---|
| source hemisphere | Southern (`SH=1`) | `run_SWIM_example.m` |
| closure | local | example and wrapper default |
| reanalysis | NCEP | example and wrapper default |
| season | annual | example; NCEP forces annual |
| pathway | pseudo-adiabatic | hard-coded in `simple_water_isotope_model_2020.m` |
| trajectory temperature step | 0.1 °C | wrapper |
| initial pressure | 101.325 kPa | evaporation / distillation thermodynamics |
| supersaturation \(a\) | 1 | example |
| supersaturation \(b\) | 0.00525 °C⁻¹ | example / paper tuned base |
| supersaturation \(c\) | 0 | example |
| default source \(\delta^{18}O_{sw}\) | −0.3 ‰ | `evaporation_2021.m` |
| cloud phase option | `adj` | `distillation_2020.m` |
| source climatology fit | spline | `evaporation_2021.m` → `T_RH_RHn_2020.m` |
| inverse interpolation | natural-neighbor-style `griddata` for method 1 | `Tsite_Tsource_reconstruction_quick.m` |

---

## 20. Open questions and paper/code discrepancies

These are not instructions to "fix" the MATLAB. They are items to resolve or lock down with reference runs before writing the parity implementation.

### OQ-1 — Which dated routines define the publication baseline?

`simple_water_isotope_model_2020.m` is the top-level function used by `run_SWIM_example.m`, but its active climatological branch calls:

```text
evaporation_2021.m
distillation_2020.m
```

The repository also contains `evaporation_2022.m` and `distillation_2022.m`.

**Action:** identify the exact functions and saved SWIM state-space file used for the paper's central results. Treat that combination as the publication baseline.

### OQ-2 — Explicit-RH branch uses a different evaporation version

Within `simple_water_isotope_model_2020.m`:

- climatological RH branch calls `evaporation_2021.m`;
- explicitly specified RH branch calls `evaporation_2020.m`.

This may be intentional or historical.

**Action:** preserve branch behavior until parity tests establish whether the distinction matters.

### OQ-3 — Likely Celsius-to-Kelvin error in `T_RH_RHn_2020.m`

The checked-in `T_RH_RHn_2020.m` contains:

```matlab
TK0 = T0 - 273.15;
SSTK0 = sst0 - 273.15;
```

inside its RHn calculation, despite the function documenting \(T_0\) and SST in °C.

`T_RH_RHn_2022.m` changes these lines to:

```matlab
TK0 = T0 + 273.15;
SSTK0 = sst0 + 273.15;
```

**Resolution:** the frozen Python parity path retains the 2020 source-condition
behavior when `evaporation_version="2021"` is selected. The Python default now
uses the corrected 2022 Kelvin conversion together with the other active
`evaporation_2022.m` diffusivity changes. A direct MATLAB evaporation reference
was captured before changing the default.

### OQ-4 — `T_RH_RHn_2022.m` declares the old function name

The file `T_RH_RHn_2022.m` begins with a function declaration named `T_RH_RHn_2020`.

This can matter in MATLAB because the primary function name and file name are expected to correspond.

**Action:** verify whether this file runs directly in the author's MATLAB environment and whether it was intended as a drop-in replacement rather than a separately callable function.

### OQ-5 — \(^{17}O\) ice equilibrium exponent differs by stage

The paper states a 0.529 equilibrium relationship in the evaporation discussion. `evaporation_2021.m` uses 0.529 for both liquid and ice, while `distillation_2020.m` uses 0.531 for ice.

**Action:** preserve observed behavior for parity and document the scientific choice before any cleanup.

### OQ-6 — Active transport diffusivity formulation

The paper discusses both fixed Jouzel–Merlivat diffusivity ratios and the Hellmann–Harvey temperature-dependent alternative. The checked-in `_2020` and `_2022` distillation files inspected here actively calculate Hellmann–Harvey ratios.

**Action:** determine which configuration generated the paper's central saved state space and reference outputs.

### OQ-7 — Stale reconstruction helpers

`Tsite_Tsource_reconstruction.m` contains calls to older/missing names such as `simple_water_istope_model_2018` and `Ts_to_Tc`, whereas `reconstruction_2020.m` actively uses `Tsite_Tsource_reconstruction_quick.m` and saved 2020 state-space files.

**Action:** treat `Tsite_Tsource_reconstruction_quick.m` plus the publication reconstruction script as the stronger candidate for the active inversion path.

---

## 21. Porting rule

For every equation above, distinguish three questions:

1. **What does the paper say the model should do?**
2. **What does the selected MATLAB baseline actually do?**
3. **What should the improved Python model eventually do?**

During the parity phase, question 2 controls numerical tests. Disagreements with question 1 should be documented, not silently repaired. Question 3 belongs to a later scientific-development phase.

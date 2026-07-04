# Water and brine thermodynamics

This document describes how the geochemistry solver computes the thermodynamic
properties of pure water and of saline solutions (brines), which ranges of
temperature, pressure and salinity are covered, and what the accuracy
limitations are.

The properties of water enter the solver in three places:

1. **Standard-state properties of aqueous species** (HKF model): the water
   density and dielectric constant, and their temperature and pressure
   derivatives, control the Born solvation terms and hence every
   temperature/pressure-dependent equilibrium constant.
2. **Activity coefficients**: the Debye-Hückel `A` and `B` parameters are
   functions of the water density and dielectric constant.
3. **Solution density**: reported as a solution property and available for
   flow coupling (buoyancy, mass/volume conversions).

## Pure water equation of state

Implemented in `opm/simulators/geochemistry/Thermo/water.{h,cpp}` (class
`water`). Given temperature `T` [K] and pressure `P` [Pa], `gibbsIAPWS(T, P)`
selects the appropriate formulation automatically:

| Conditions                                              | Formulation                       | `region_` |
|---------------------------------------------------------|-----------------------------------|-----------|
| 273.15–623.15 K, P<sub>sat</sub> ≤ P ≤ 100 MPa           | IAPWS-IF97 region 1 (liquid)      | 1         |
| 273.15–1073.15 K, P below P<sub>sat</sub>/B23, ≤ 100 MPa | IAPWS-IF97 region 2 (steam)       | 2         |
| 623.15–863.15 K, above the B23 line, ≤ 100 MPa           | IAPWS-IF97 region 3 (Helmholtz)   | 3         |
| 273.15–1073.15 K, 100 MPa < P ≤ 1000 MPa                 | IAPWS-95 (scientific formulation) | 95        |

The saturation line (IF97 region 4) is available as `water::PsatIAPWS(T)` for
273.15 K ≤ T ≤ 647.096 K, and the boundary between regions 2 and 3 as
`water::PB23IAPWS(T)`.

For every state the class provides the specific volume, internal energy,
enthalpy, entropy, heat capacities, speed of sound and Gibbs free energy, plus
the derivative properties needed by the HKF machinery: the isobaric thermal
expansivity `alpha_` [1/K], its temperature derivative `alpha_t_`, and the
isothermal compressibility `beta_` [1/Pa].

Implementation notes:

- Regions 1 and 2 are Gibbs-energy formulations evaluated directly from
  (T, P). Region 3 and IAPWS-95 are Helmholtz formulations in (density, T);
  the density is obtained from P with a bisection-safeguarded Newton
  iteration. Below the critical temperature the iteration is started on the
  correct side of the two-phase dome using the Wagner-Pruss auxiliary
  saturated-density correlations.
- Conditions outside the covered ranges raise `std::domain_error`
  (T < 273.15 K, T > 1073.15 K, P ≤ 0, P > 1000 MPa). IF97 region 5
  (T > 1073.15 K) is not implemented.
- **The melting line is not checked.** At low temperature and very high
  pressure (roughly P > 632 MPa at 273 K, the ice VI stability field) the
  returned properties refer to the metastable liquid.
- The IF97 and IAPWS-95 formulations are mutually consistent to within a few
  hundredths of a percent, so crossing the 100 MPa seam does not produce
  significant jumps in density or its derivatives.

Verification (see `tests/tests_using_catch2/test_thermo.cpp`): the
implementation reproduces the reference tables of the IAPWS releases — IF97
Table 5 (region 1), Table 15 (region 2), Table 33 (region 3), Table 35
(saturation line), and IAPWS-95 release Table 7 (700 MPa states at 300, 500
and 900 K) — as well as property continuity across all internal region
boundaries.

## Dielectric constant and Born functions

The static dielectric constant of water and the Born functions Z, Q, Y, X are
computed with the Johnson & Norton (1991) model
(`Thermo/eps_JN.{h,cpp}`), a polynomial in water density with
temperature-dependent coefficients. It is calibrated for 0–1000 °C and
1–5000 bar. Above 500 MPa the model extrapolates in density; the trend
remains smooth and physically reasonable up to 1000 MPa, but the Born
derivative functions carry increasing uncertainty there.

The HKF g-function correction for effective ionic radii (Shock et al., 1992,
`Thermo/ions.{h,cpp}`) is active only in the low-density window
155–355 °C between the saturation pressure and 1 kbar, and is identically
zero at liquid-like densities (ρ > 1 g/cm³), i.e. at high pressure.

## Solution (brine) density

`BasVec::solution_density()` (`Core/ChemBasVec.cpp`) computes the density of
the saline solution from the speciated composition. Per kg of solvent water,

```
mass   = 1 + Σ m_i M_i          [kg]
volume = 1/ρ_w + Σ m_i V°_i     [m³]
ρ_sol  = mass / volume          [kg/m³]
```

where the sums run over all aqueous species (basis species and complexes),
`m_i` are the molalities, `M_i` the molecular weights and `V°_i` the HKF
standard molal volumes at the current (T, P) as computed by `hkf::dGIons`.
The value is exposed as the `Solution_density` entry of the key solution
properties, next to `Water_density` (the pure-water density).

Design rationale:

- The sums use the **speciated** molalities, so ion pairs (e.g. NaCl°,
  CaCl⁺) contribute their own standard volumes. Since explicit complexation
  carries part of the solution non-ideality in this framework, this recovers
  part of the excess volume of mixing.
- The remaining approximation — standard (infinite-dilution) molal volumes
  with no Pitzer-type excess terms — is consistent with the extended
  Debye-Hückel activity model used by the solver: both degrade at comparable
  ionic strengths.
- Species without HKF volume parameters contribute zero volume; if the
  database provides no molecular weights, the pure-water density is returned.

Accuracy against measured NaCl solution densities at 25 °C / 1 bar
(`tests/tests_using_catch2/test_equilibrium_solver.cpp`):

| Salinity        | Measured [kg/m³] | Deviation      |
|-----------------|------------------|----------------|
| dilute limit    | 997.05           | exact          |
| 1 molal NaCl    | ~1036            | below 1%       |
| 4 molal NaCl    | ~1140            | below 3%       |

The deviation grows with molality because of the neglected excess volume;
above ~3 mol/kg ionic strength the density (like the activity model) should
be considered approximate at the percent level.

## Units

| Quantity                              | Storage/interface unit |
|---------------------------------------|------------------------|
| Temperature                           | K (input files: °C)    |
| Pressure                              | Pa (input files: Pa; tables: bar) |
| `ChemTable::mol_weight_`              | kg/mol                 |
| `ChemTable::mol_volume_` (aqueous)    | m³/mol                 |
| `StandardStateProperties::V`          | m³/mol                 |
| `Water_density`, `Solution_density`   | kg/m³                  |

Note (July 2026): two unit errors in the previously unused molal-volume
output of `hkf::dGIons`/`hkf::ionProperties` were fixed — a factor 1000 in
the overall conversion (values were L/mol instead of m³/mol) and a spurious
factor 10⁵ on the Born g-function pressure derivative inside its 155–355 °C
window. Thermo tables generated before this fix show aqueous-species volumes
that are wrong by these factors; Gibbs energies, entropies, heat capacities
and log K values were not affected. See
[hkf_molal_volume_unit_fixes.md](hkf_molal_volume_unit_fixes.md) for the
full analysis.

## References

- IAPWS, *Revised Release on the IAPWS Industrial Formulation 1997 for the
  Thermodynamic Properties of Water and Steam* (IAPWS-IF97), Lucerne, 2007.
- W. Wagner and A. Pruss, *The IAPWS Formulation 1995 for the Thermodynamic
  Properties of Ordinary Water Substance for General and Scientific Use*,
  J. Phys. Chem. Ref. Data **31** (2002) 387–535.
- J. W. Johnson and D. Norton, *Critical phenomena in hydrothermal systems*,
  Am. J. Sci. **291** (1991) 541–648.
- E. L. Shock, E. H. Oelkers, J. W. Johnson, D. A. Sverjensky and
  H. C. Helgeson, *Calculation of the thermodynamic properties of aqueous
  species at high pressures and temperatures*, J. Chem. Soc. Faraday Trans.
  **88** (1992) 803–826.
- J. C. Tanger and H. C. Helgeson, *Calculation of the thermodynamic and
  transport properties of aqueous species at high pressures and
  temperatures*, Am. J. Sci. **288** (1988) 19–98.

# HKF molal-volume unit fixes (July 2026)

Two unit errors in the standard molal volume of aqueous species were found
and fixed while implementing the salinity-dependent solution density. This
note documents what was wrong, what was affected, and — importantly — what
was **not** affected, as a reference for code review.

Both errors were in the volume expression shared by `hkf::dGIons` (the `MV`
output array) and `hkf::ionProperties` (the `StandardStateProperties::V`
field) in `opm/simulators/geochemistry/Thermo/hkf.cpp`.

## Bug 1: overall conversion factor 1000× too large

The HKF volume terms are accumulated in J/(mol·bar) and converted to SI with
the factor `Chat`:

```
before:  Chat = 41.84e-3 / cal2J   (= 1e-2)
after:   Chat = 41.84e-6 / cal2J   (= 1e-5, i.e. m^3*bar/J)
```

J/(mol·bar) = 1e-5 m³/mol, so the correct factor is 1e-5. With the old
factor the result was in **L/mol** while both the code documentation and the
consumers assumed m³/mol.

Empirical confirmation: with the fix, the computed conventional standard
molal volume of Na⁺ at 25 °C / 1 bar is −1.207 cm³/mol, in line with the
literature value of about −1.1 to −1.2 cm³/mol (the small spread reflects
the dielectric model; SUPCRT92 gives −1.11). Before the fix the same
quantity came out as −1206.7 "cm³/mol".

## Bug 2: Born g-function pressure derivative converted twice

The volume expression contains the pressure derivative of the effective Born
coefficient, `w_P`. The comment in the code assumed `w_P` was in J/(mol·Pa)
and multiplied by 1e5 to get to J/(mol·bar):

```
before:  ... - 1e5*omega*bornQ - (bornZ + 1)*w_P*1e5
after:   ... - 1e5*omega*bornQ - (bornZ + 1)*w_P
```

However, `w_P` is built from `born_g_P_`, which `ions::born_df` already
returns **per bar** (it converts the water compressibility with `beta_*1e5`
and differentiates the f-function with respect to pressure in bar). The
extra 1e5 therefore inflated this contribution by 100 000.

The `bornQ` term is different: `eps_JN` computes Q in 1/Pa, so its 1e5
factor is correct and is kept.

The g-function (Shock et al. 1992) is identically zero outside the window
155 °C < T < 355 °C, P<sub>sat</sub> < P < 1 kbar and at liquid-like
densities, so bug 2 was invisible at room temperature and only exploded
inside that window — e.g. V(Cl⁻) at 200 °C / 500 bar came out as
29 193 cm³/mol instead of ≈ 15 cm³/mol.

## What was affected

- **The V column of thermo tables** (`ThermoTable`, the recent thermo-table
  feature) for **aqueous species**: values 1000× too large everywhere, and
  additionally distorted by bug 2 inside the g-function window. Tables
  generated before the fix should be regenerated if the aqueous V column is
  used. Mineral and water rows were correct (they come from different code
  paths: the mineral molar volume is a database value, water uses the
  IAPWS equation of state directly).
- **The `MV` output of `hkf::dGIons`** (per-species molal volumes stored in
  `ChemTable::mol_volume_` for the aqueous tables): this output was dead
  code — computed but consumed by nothing — until the solution-density
  feature. No earlier result depended on it.

## What was NOT affected

**Equilibrium constants (logK) and all other standard-state properties were
never affected, at any time.** The reasons:

1. The Gibbs free energy in both `dGIons` and `ionProperties` is a separate
   expression that uses neither `Chat` nor `w_P`; its Born contribution
   (`(Wi - omega)*W + omega*Wref`) has always been in consistent units.
   logK values are derived from these Gibbs energies.
2. The entropy uses `w_T` and the heat capacity `w_TT` — temperature
   derivatives with no pressure-unit conversion involved.
3. The pressure dependence of the Gibbs energy is carried by the
   `a1..a4` terms and the Born term directly, not by integrating the molal
   volume.

Empirical confirmation: after the fix, all logK-based regression tests
(RunGeoChem transport suites, equilibrium reference cases, the thermo-table
pins for G, H, S, Cp of Na⁺ and the calcite logK) pass unchanged. The only
test expectation that had to be updated is the direct pin of the Na⁺ V
value in `tests/tests_using_catch2/test_thermo.cpp`, which had recorded the
pre-fix number.

## How the bugs were found

The solution-density feature sums `m_i * V°_i` over the speciated
composition (see [water_thermodynamics.md](water_thermodynamics.md)). The
resulting density of a 1 molal NaCl solution at 25 °C was ~15 kg/m³ instead
of ~1036 kg/m³ (bug 1), and a case at 200 °C / 500 bar remained absurd after
the first correction (bug 2). Densities now match measured NaCl data to
within 1% at 1 molal and 3% at 4 molal, and the individual ion volumes match
literature values, which pins down both conversion factors independently.

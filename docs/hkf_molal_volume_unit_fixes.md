# HKF molal-volume fixes (July 2026)

Two unit errors and one Born-coefficient error in the standard molal volume
of aqueous species were found and fixed while implementing and validating
the salinity-dependent solution density. This note documents what was wrong,
what was affected, and — importantly — what was **not** affected, as a
reference for code review.

All three errors were in the volume expression shared by `hkf::dGIons` (the
`MV` output array) and `hkf::ionProperties` (the
`StandardStateProperties::V` field) in
`opm/simulators/geochemistry/Thermo/hkf.cpp`.

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
before:  ... - 1e5*omega_ref*bornQ - (bornZ + 1)*w_P*1e5
after:   ... - 1e5*omega_ref*bornQ - (bornZ + 1)*w_P
```

However, `w_P` is built from `born_g_P_`, which `ions::born_df` already
returns **per bar** (it converts the water compressibility with `beta_*1e5`
and differentiates the f-function with respect to pressure in bar). The
extra 1e5 therefore inflated this contribution by 100 000 wherever the
g-function pressure derivative was nonzero.

The `bornQ` term is different: `eps_JN` computes Q in 1/Pa, so its 1e5
factor is correct and is kept.

The g-function (Shock et al. 1992) contains a density-dependent term and an
empirical difference function `f`. Only `f` is restricted to the window
155 °C < T < 355 °C and P<sub>sat</sub> < P < 1 kbar; the density term can
remain nonzero whenever ρ < 1 g/cm³. Bug 2 was negligible near room
conditions but grew rapidly as water density fell. For example, V(Cl⁻) at
200 °C / 500 bar came out as 29 193 cm³/mol instead of approximately
15 cm³/mol.

## Bug 3: the Born Q term used the reference omega

The Gibbs-energy Born term is

```text
(Wi - omega_ref)*W + omega_ref*Wref
```

where `Wi` is the effective Born coefficient at the current temperature and
pressure. Differentiating this expression with respect to pressure gives the
volume contribution

```text
-Wi*Q - (bornZ + 1)*w_P
```

but the volume code used `omega_ref` in the Q term:

```text
before:  ... - 1e5*omega_ref*bornQ - (bornZ + 1)*w_P
after:   ... - 1e5*Wi*bornQ        - (bornZ + 1)*w_P
```

The difference is negligible at the reference state, where `Wi` is
effectively `omega_ref`, but grows as the water density and effective ionic
radius change. At 400 °C and 3000 bar, for example, the Na⁺ volume changes
from 1.01110 to 1.00869 cm³/mol. The corrected value agrees with a centered
finite difference of the implemented Gibbs energy.

## Related fix: neutral species in the scalar path

The same consistency test showed that `hkf::ionProperties` initialized `Wi`
to zero for neutral species. In revised HKF, a neutral species uses its
pressure- and temperature-independent reference Born coefficient, so the
correct initialization is `Wi = omega_ref` with zero `w_T`, `w_TT`, and
`w_P`. The bulk `dGIons` Gibbs expression already represented this
convention.

This correction aligns scalar thermo-table rows with the bulk solver for
neutral aqueous complexes. Unlike the three volume-expression defects, it
can also change scalar-path G/H/S/Cp away from the reference state when the
neutral species has a nonzero `omega_ref`.

## What was affected

- **The V column of thermo tables** (`ThermoTable`, the recent thermo-table
  feature) for **aqueous species**: values 1000× too large everywhere, and
  additionally distorted by bugs 2 and 3 away from the reference state.
  Tables generated before the fixes should be regenerated if the aqueous V
  column is used. Mineral and water rows were correct (they come from
  different code paths: the mineral molar volume is a database value, water
  uses the IAPWS equation of state directly).
- **The `MV` output of `hkf::dGIons`** (per-species molal volumes stored in
  `ChemTable::mol_volume_` for the aqueous tables): this output was dead
  code — computed but consumed by nothing — until the solution-density
  feature. No earlier result depended on it.
- **Neutral aqueous rows in the scalar thermo-table path** could have
  inconsistent G/H/S/Cp/V values away from the reference state. The bulk
  equilibrium path already used the constant neutral-species Born
  coefficient.

## What was NOT affected

For charged species, equilibrium constants and the non-volume standard-state
properties were not affected by the three volume-expression defects:

1. The Gibbs free energy in both `dGIons` and `ionProperties` is a separate
   expression that uses neither `Chat` nor `w_P`; its Born contribution
   (`(Wi - omega)*W + omega*Wref`) has always been in consistent units.
   logK values are derived from these Gibbs energies.
2. The entropy uses `w_T` and the heat capacity `w_TT` — temperature
   derivatives with no pressure-unit conversion involved.
3. The pressure dependence of the Gibbs energy is carried by the
   `a1..a4` terms and the Born term directly, not by integrating the molal
   volume.

The related neutral-species scalar correction does not change the bulk
equilibrium path. All logK-based regression tests (RunGeoChem transport
suites, equilibrium reference cases, the thermo-table pins for G, H, S, Cp
of Na⁺ and the calcite logK) pass unchanged.

## How the bugs were found

The solution-density feature sums `m_i * V°_i` over the speciated
composition (see [water_thermodynamics.md](water_thermodynamics.md)). The
resulting density of a 1 molal NaCl solution at 25 °C was ~15 kg/m³ instead
of ~1036 kg/m³ (bug 1), and a case at 200 °C / 500 bar remained absurd after
the first correction (bug 2). The remaining coefficient mismatch was exposed
by checking the reported volume against a numerical pressure derivative of
the same Gibbs energy (bug 3). Comparing scalar and bulk HKF evaluations then
exposed the neutral-species initialization. Densities now match measured
NaCl data to within 1% at 1 molal and 3% at 4 molal, and the individual ion
volumes match literature values, which pins down both conversion factors
independently.

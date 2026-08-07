# Appelo et al. (2014) volume validation

These inputs exercise the water EOS and HKF standard molal volumes against
Figures 1-3 in:

C. A. J. Appelo, D. L. Parkhurst, and V. E. A. Post, "Equations for
calculating hydrogeochemical reactions of minerals and gases such as CO2 at
high pressures and temperatures", *Geochimica et Cosmochimica Acta* 125
(2014), 49-67, DOI: 10.1016/j.gca.2013.10.003.

Build the command-line solver first:

```sh
cmake --build build --target GeoChemX -j2
```

## What can be reproduced

| Paper figure | Validation in this directory |
| --- | --- |
| Figure 1 | Current SUPCRT-coefficient values at `I=0`; not the full HKFmoRR curves. |
| Figure 2 | The dashed SUPCRT intrinsic-volume curves and an independent SUPCRT92 reference. |
| Figure 3 | The explicitly documented worst-fit point in panel B. |

The current implementation evaluates HKF standard molal volumes at infinite
dilution. It does not implement the ionic-strength terms of the modified
Redlich-Rosenfeld equation (paper Eq. 6). Therefore, it cannot reproduce the
concentration-dependent slopes in Figure 1 or the full HKFmoRR density results
in Figure 3. Those require another model change, not different input data.

## Figure 1

Run:

```sh
build/GeoChemX THERMOTABLE \
  docs/validation/appelo_2014/fig1_infinite_dilution.in
```

At 25 C and 1 atm the current results are:

| Ion | `V0` [cm3/mol] |
| --- | ---: |
| H+ | 0.00000 |
| Na+ | -1.20668 |
| Ca+2 | -18.43403 |
| Cl- | 17.34629 |

Use additivity at infinite dilution:

```text
V0(HCl)   = V0(H+)   + V0(Cl-)     = 17.34629 cm3/mol
V0(NaCl)  = V0(Na+)  + V0(Cl-)     = 16.13961 cm3/mol
V0(CaCl2) = V0(Ca+2) + 2*V0(Cl-)   = 16.25855 cm3/mol
```

Paper Table 2 gives fitted HKFmoRR intrinsic volumes of approximately
`V0(Na+)=-1.52`, `V0(Ca+2)=-18.3`, and `V0(Cl-)=18.0 cm3/mol`. Those imply
Figure 1 intercepts of 18.0, 16.48, and 17.7 cm3/mol for HCl, NaCl, and
CaCl2, respectively. Thus, the values above are checks of the current
SUPCRT coefficient set, not exact reproductions of the paper's PHREEQC
intercepts. Running the solution solver over a concentration sweep also will
not produce the rising curves because the Eq. 6 excess-volume terms are
absent.

## Figure 2

`fig2_water_saturation_pressures.in` reports the IAPWS saturation-pressure
path. The ion input uses 1 atm through 95 C and `Psat + 0.001 bar` from
100-200 C to stay on the liquid side of the phase boundary. The first
temperature is 0.01 C to match the SUPCRT92 default saturation run.

Generate the table and a plot-friendly TSV:

```sh
build/GeoChemX THERMOTABLE \
  docs/validation/appelo_2014/fig2_supcrt_intrinsic_volumes.in \
  > /tmp/appelo_fig2.out

awk 'BEGIN {OFS="\t"; print "species","T_C","P_bar","V_cm3_mol"}
     $1 == "SPECIES" {species=$2; next}
     species != "" && $1 ~ /^[0-9]+$/ {print species,$2,$3,$9}' \
  /tmp/appelo_fig2.out > /tmp/appelo_fig2.tsv
```

Plot `V_cm3_mol` against `T_C`, grouped by `species`. These are the dashed
SUPCRT curves in Figure 2, not the solid PHREEQC HKFmoRR curves.
`fig2_supcrt92_reference.tsv` was generated independently with SUPCRT92 1.1
using the distribution database, the default liquid saturation path, and each
ion as a unit-stoichiometry reaction. The corresponding SUPCRT reaction input
is `fig2_supcrt92_reference.rxn`; the values are from its generic `.vxy`
volume output.

Compare all common samples numerically:

```sh
awk -f docs/validation/appelo_2014/compare_fig2.awk \
  docs/validation/appelo_2014/fig2_supcrt92_reference.tsv \
  /tmp/appelo_fig2.tsv
```

Overlay the reference points and GeoChemX curves with:

```sh
gnuplot docs/validation/appelo_2014/plot_fig2.gnuplot
```

This writes `/tmp/appelo_fig2.png`. At the common 25 C increments from
25-200 C, the maximum absolute GeoChemX-SUPCRT92 difference is
0.08933 cm3/mol (`SO4-2` at 200 C). Small residuals are expected because
GeoChemX uses the IAPWS water EOS and a slightly different pressure path.
The 25 C GeoChemX checkpoint is:

| Species | `V0` [cm3/mol] |
| --- | ---: |
| Na+ | -1.20668 |
| K+ | 9.00998 |
| Mg+2 | -22.00790 |
| Ca+2 | -18.43403 |
| Cl- | 17.34629 |
| HCO3- | 24.21563 |
| SO4-2 | 12.92941 |

## Figure 3

Appelo et al. identify the largest Figure 3 difference explicitly: a solution
with 0.93 NaCl and 1.36 CaCl2 at 5 C has measured density 1154.1 kg/m3.
The concentrations are used directly as mol/kg water in the input. This is
also the interpretation consistent with the paper's stated measured
mixed-salt volume:

```text
solute mass = 0.93 * 58.44277 + 1.36 * 110.984
            = 205.2900 g/kg water

mixed-salt amount = 0.93 + 2 * 1.36
                  = 3.65 mol/kg water

Vm(measured) = ((1000 + 205.2900) / 1.1541
                - 1000 / 0.999967) / 3.65
             = 12.14 cm3/mol
```

Run `fig3b_documented_point.dat` from a temporary directory because `EQSOLVER`
writes result files next to its input:

```sh
repo=$PWD
work=$(mktemp -d)
cp docs/validation/appelo_2014/fig3b_documented_point.dat "$work/"
(cd "$work" && "$repo/build/GeoChemX" EQSOLVER fig3b_documented_point.dat)
```

Current result:

```text
Water_density[kg/m3]    = 999.967
Solution_density[kg/m3] = 1166.75
Measured density        = 1154.1
Difference              = +12.65 kg/m3 (+1.10%)

Current mixed-salt Vm   = 9.04 cm3/mol
Measured mixed-salt Vm  = 12.14 cm3/mol
Paper HKFmoRR Vm        = 12.9 cm3/mol
```

The calculation converges, but its ionic strength is 4.4347 mol/kg. This is
well into the range where both the extended Debye-Huckel activity model and
the standard-volume density approximation are expected to degrade. The paper's
HKFmoRR result differs from the measurement by only 3.1 kg/m3, demonstrating
the effect of its concentration correction.

The full Figure 3 scatter plot additionally needs the composition tables from
the four experimental sources named in the paper caption. They are not
included in the Appelo et al. paper.

## EOS boundary checks

The two additional inputs isolate the EOS behavior from the paper comparison:

```sh
build/GeoChemX THERMOTABLE \
  docs/validation/appelo_2014/water_eos_states.in

build/GeoChemX THERMOTABLE \
  docs/validation/appelo_2014/hkf_volume_states.in
```

`water_eos_states.in` covers steam and liquid water at 150 C, region 3, the
999/1000 bar implementation boundary, and the high-pressure extension.
`hkf_volume_states.in` verifies finite charged-species volumes through
7000 bar, including the 1000 and 3000 bar points that previously exposed
derivative problems.

To check e.g. water density values from `water_eos_states.in` agains Python do
```
#pip install iapws
from iapws import IAPWS95, IAPWS97

Ts = [25, 150, 150, 200, 350, 400, 400, 400, 400, 400]
Ps = [1.01325, 1, 10, 500, 200, 250, 999, 1000, 3000, 7000]

for T_C, P_bar in zip(Ts, Ps):
  eos = IAPWS97 if P_bar <= 1000 else IAPWS95
  water = eos(P=P_bar / 10, T=T_C + 273.15)
  print(
    f"{T_C:5.0f} C {P_bar:8.3f} bar "
    f"{eos.__name__:7s} "
    f"{water.rho:10.4f} kg/m3 "
    f"{water.rho * 1e-3:.7f} g/cm3"
    )

```

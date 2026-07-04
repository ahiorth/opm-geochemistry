# Tutorials

Below are tutorials, and in some cases documentations, for the geochemistry solver itself and the coupling to OPM Flow.
The intention is to provide examples of the capabilities of the simulators and provide input files/simulation decks.

See also the reference documentation of the internal models:

- [Water and brine thermodynamics](../water_thermodynamics.md) — pure-water equation of state
  (IAPWS-IF97/IAPWS-95), validity ranges, dielectric constant, and salinity-dependent solution density.

## Documentation and tutorials for the geochemistry solver

The following tutorials showcase different aspects of the geochemistry solver, with some associated documentation of its
internal workings, and is associated with a corresponding keyword(s) for the standalone program `GeoChemX`:

- [Seawater speciation - `SOLUTION` keyword](./geochemx/solution/notebook/main_solution.ipynb)
- [Equilibrium with exchange sites - `IEXCHANGE` keyword](./geochemx/iexchange/notebook/main_iexchange.ipynb)
- [Equilibrium between seawater, CO2 and calcite - `EQUILIBRIUM_PHASES`
  keyword](./geochemx/equilibrium_phases/notebook/main_equilibrium_phases.ipynb)
- [Modeling gas phase at constant volume - `GAS_PHASE` keyword](./geochemx/gas_phase/notebook/main_gas_phase.ipynb)
- [Surface complexation - `COMPLEX` keyword](./geochemx/complex/notebook/main_complex.ipynb)
- [Core flooding - `RATE` and `INTERPOLATE` keywords](./geochemx/rate/notebook/main_rate.ipynb)
- [1D reactive transport](./geochemx/transport/notebook/main_transport.ipynb)
- [Advanced reactive transport example: Adsorption of sulphate and cation exchange in
  core](./geochemx/transportII/notebook/main_transportII.ipynb)

## Tutorials for OPM Flow coupled with geochemistry solver

The following tutorials showcase the coupled OPM Flow-geochemistry simulator, where the focus is on setting up and
running OPM Flow decks:

- [PHREEQC Example 11](./opm_flow/phreeqc_ex_11/PHREEQC_ex_11.md)
- [Mineral reactions](./opm_flow/mineral_interaction/magnesite_calcite_example.md)
- [Adsorption of sulphate and cation exchange in core](./opm_flow/adsorption/transportII.md)

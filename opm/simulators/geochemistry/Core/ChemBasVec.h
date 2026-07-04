/*
 * MIT License
 *
 * Copyright (C) 2025 Aksel Hiorth
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
*/
#ifndef CHEM_BAS_VEC_H
#define CHEM_BAS_VEC_H

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include <memory>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include <opm/simulators/geochemistry/Core/ChemInitChem.h>
#include <opm/simulators/geochemistry/Core/ChemTable.h>
#include <opm/simulators/geochemistry/Core/MineralKinetics.hpp>
#include <opm/simulators/geochemistry/Core/SurfaceChemistry.hpp>

#include <opm/simulators/geochemistry/Common/Constants.hpp>
#include <opm/simulators/geochemistry/Common/Enums.hpp>
#include <opm/simulators/geochemistry/Common/SerializeForTesting.hpp>
#include <opm/simulators/geochemistry/IO/ErrorMsg.hpp>
#include <opm/simulators/geochemistry/IO/FileFunctions.hpp>
#include <opm/simulators/geochemistry/IO/GeoLogger.hpp>
#include <opm/simulators/geochemistry/LinearAlgebra/LinearSolverRoutines.hpp>
#include <opm/simulators/geochemistry/Thermo/hkf.h>

// Forward declarations
class GCSolver;
class InitChem;

class BasVec {

  public:

    static BasVec* createFromInitChem(InitChem* ICS);
    ~BasVec();

    int dimension() const;

    void update_ctot_and_ctot_mineral(const std::vector<double>& dc, std::vector<double>& dc_min);

    // Helper functions which can be used by multiple output routines.
    void update_before_write_solution();
    void write_solution(const std::string& name, DiffuseLayerProperties* EDL_props=nullptr);
    void write_tot_conc(const std::string& name);
    void write_info(BasVecInfo* info);
    void write(const std::string& name);
    void write_solution_chemistry(const std::string& name,bool include_reactions=false);
    void write_gas_phase_solution(const std::string& name);
    void write_buffers(const std::string& name, bool include_reactions = false);

    // NB: For secondary species, must be called after update_before_write_solution().
    [[nodiscard]] double get_species_concentration(const std::string& species_name) const;
    [[nodiscard]] double get_delta_mineral(const std::string& mineral_name) const;

    void calc_ctot_aq();

    /** Updates concentration of the rock buffered species (if there are any). Uses the basis transformation matrix: log_10 a = beta^-1 *( log_10 a' + log_10 K). */
    void calc_rock_spec_conc();

    void calculate_gas_phase_pressure();

    void calc_complex_sparse(ChemTable& sm) const;
    double update_exchange_capacity_dl_model(double frac_DL);
    double update_ionic_strength();

    /* Solution (brine) density [kg/m^3] from the speciated composition and the
     * HKF standard molal volumes at the current (T, P); see the implementation
     * for details. */
    double solution_density() const;

    void add_new_buffer_mineral(int ICS_rock, int ICS_buffer);

    BasVec* convert_nlin_rate_equations_to_equilibrium(const double* cmin);

    [[nodiscard]] int index_of_most_supersaturated_mineral(const std::vector<int>& mineral_list,
                                                           int size,
                                                           const std::vector<int>& old_mineral_list,
                                                           int old_size) const;

    void update_mineral_concentrations_kinetics(double dt, double* c_min);
    void update_mineral_concentrations_equilibrium(const double* c0_min, const std::array<double, 3>& mass_phase, double* c_min);
    void update_ioxch_complex();

    void update_dG(hkf& HKF_EOS, double T, double P);
    void update_logK(hkf& HKF_EOS, double T);
    void update_molar_volume_gases_equilibrium_phases();

    void set_temperature_for_equilibrium_reactions(hkf& HKF_EOS);
    void set_temperature_for_mineral_kinetics();
    void update_activity_coefficients(double Io);

    void add_gas_phase_to_total_conc();

    [[nodiscard]] bool new_temperature() const;
    [[nodiscard]] bool new_pressure() const;

    // Note: When using these functions, it is assumed that the appropriate variables have already been updated!
    [[nodiscard]] double calc_pH() const;
    [[nodiscard]] double calc_exchanger_charge() const;  // eq/L
    [[nodiscard]] double calc_solution_charge() const;  // eq/L
    [[nodiscard]] double calc_surface_charge(bool coulombs_per_square_meter=false) const;  // If false, return surface charge in eq/L

  public:

    // Geochemical database
    BasVec* Vp_eq_{nullptr};  // All non-linear rate minerals replaced with equilibrium
    InitChem* ICS_{nullptr};

    bool equilibrate_{false};
    bool charge_balance_{false}; // Calculate pH from charge-balance?

    int surface_flag_{0}; // 0: No separate surface equation, 1: Extra surface equation (immediately after the mass balance equations)
    SurfaceChemistryOptions surface_options_{};

    // Compared to InitChem: No max_size_, size_aq_, size_io, or size_surf_.
    //
    int size_{0};  // Number of basis species
    int size_mass_{0};  // Size of mass balance
    int size_min_{0};  // Size of minerals and rock buffer
    int size_sup_min_{0};  // size of supersaturated minerals
    int size_rock_{0};  // Size of rock buffer
    int size_gas_{0};  // Size of gas buffers
    int size_gase_phase_{ 0 }; // number of components in a gas pahse
    int size_all_gas_buffers_{ 0 }; // size of all gas phases present independent on model

    // Compared to InitChem: No pos_. Additionally: pos_rel_sup_min_ and pos_rel_gas_.
    //
    std::vector<int> pos_mass_;  // Pos. to the species determined by mass balance
    std::vector<int> pos_buffer_;   // Pos. to rock buffer basis
    std::vector<int> pos_rock_;  // Pos. to rock buffered species
    std::vector<int> pos_min_;  // Pos. to total mineral buffer
    std::vector<int> pos_sup_min_;  // Pos. to supersaturated minerals
    std::vector<int> pos_min_bas_;  // Pos. to total basis mineral buffer
    std::vector<int> pos_sup_bas_;  // Pos. to basis species for supersaturated minerals
    std::vector<int> pos_gas_;  // Pos. to gas or oil buffer
    std::vector<int> pos_gas_phase_; // pos to gas phases

    std::vector<int> pos_all_gas_buffers_; // any gas buffer in the simulation, independent of model used

    std::vector<int> pos_rel_sup_min_;  // Pos. relative to input file
    std::vector<int> pos_rel_gas_;  // Pos. relative to pos_rock

    std::vector<int> phase_;  // = 0 water, 1=oil, 2=gas

    int pos_water_{-1};  // H2O
    int pos_pH_{-1};  // H+
    int pos_pe_{ -1 };  // Redox stuff
    int pos_exp_{-1};  // Pos. to surface potential basis specie, E

    bool has_ion_exchange_{false};
    std::vector<int> pos_exch_bas_; // Pos. to ion exchange basis species
    std::vector<int> pos_exch_comp_; // Pos. to ion exchange secondary complexes

    double WHTO_{1.0};  // Concentration of water (H2O)
    double Io_{0.0};  // Ionic strength

    double sch_{0.0};  // Surface charge
    double psi_{0.0};  // Surface potential
    double d_DL_{0.0};  // Size of diffuse layer (nm)
    double SA_{1.0};  // Specific surface area surface species (m^2/L)

    std::array<double, 3> Temp_{0.0, 0.0, 0.0};  // Phase-temperatures
    std::array<double, 3> Pres_{0.0, 0.0, 0.0};  // Phase-pressures

    std::vector<std::vector<double>> M_all_trans_;  // Transformed stoichiometric matrix for complexes

    double** beta_bas_inv_{nullptr};  // Inverse of basis switching matrix
    double** sm_buf_{nullptr};  // Reduced stoichiometric matrix of basis buffers
    double** sm_min_{nullptr};  // Reduced stoichiometric matrix of basis sup sat min and rock (eq) buffered

    // Total calculated concentrations of basis species
    std::array<std::vector<double>, GeochemicalComponentType::NUMBER_OF_NON_MINERAL_TYPES_> ctot_calc_;

    std::vector<double> ctot_;  // Total concentrations of basis species
    std::vector<double> ctot_ads_; // Total concentrations of adsorbed species
    std::vector<double> ctot_dl_excess_;  // Concentrations in the diffuse layer
    std::vector<double> delta_mineral_;  // Amounts of mineral precipitation
    std::vector<double> ctot_mineral_;  // Total concentrations of mineral phases (usually of a different size than delta_mineral)
    std::vector<double> log_m_;  // log10-concentrations for basis species
    std::vector<double> log_a_;  // log10-activities for basis species
    std::vector<double> log_g_;  // log10-activity coefficients
    std::vector<double> log_a_gas_;  // log10- activity for minerals

    std::vector<double> gas_phase_pressure_; // partial pressure gas phases
    std::vector<double> c_gas_phase_; // concentration of gas: (NB) mole gas molecules /volume of water 
    std::vector<double> mol_fraction_gas_phase_; 
    double gas_phase_frac_{ 1. };

    std::vector<MineralRateLaw> rate_laws_;  // Rate constants for each mineral phase

  private:

    std::vector<double> log_bas_;

  private:

    /** Helper function to find the basis surface complex corresponding to a
     * given secondary surface complex. We assume there is exactly one such complex. */
    [[nodiscard]] std::size_t indexOfBasisSurfaceComplex(std::size_t indexOfSecondarySurfaceComplex) const;

    void assign_basis2mineral();

    /** Computes basis transformation matrix and reduced stoichiometric matrix for basis buffers. */
    void update_basis_buffer_matrices(int old_size_rock=0);
    void update_sm_buf(int size_r, int* pos_b, int* pos_r, double** sm);
    void update_beta_inv();

    BasVec(InitChem* ICS);  // Note: Constructor is no longer called directly.

};

#endif

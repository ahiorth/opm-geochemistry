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
#ifndef WATER_H
#define WATER_H

#include <array>
#include <cmath>
#include <iostream>

/*
 * The International Association for the Properties of Water and Steam Lucerne, Switzerland
 * August 2007, Revised Release on the IAPWS Industrial Formulation 1997 (IAPWS-97).
 * (The revision only relates to the extension of region 5 to 50 MPa)
 *
 * Implemented regions:
 *
 *   - Region 1: compressed liquid, 273.15 K <= T <= 623.15 K, Psat(T) <= P <= 100 MPa
 *   - Region 2: steam/superheated vapour, 273.15 K <= T <= 1073.15 K,
 *               P below the saturation line (T <= 623.15 K) or below the B23
 *               boundary (623.15 K < T <= 863.15 K), P <= 100 MPa
 *   - Region 3: near-critical/supercritical, 623.15 K < T <= 863.15 K, above the
 *               B23 boundary, P <= 100 MPa (Helmholtz formulation, density is
 *               obtained with a safeguarded Newton iteration)
 *   - Region 4: the saturation line (PsatIAPWS)
 *
 * Region 5 (T > 1073.15 K) is not implemented. Conditions outside the covered
 * ranges raise std::domain_error.
 */
class water
{

public:

    water();

	void gibbsIAPWS(double T, double P);
	static double PsatIAPWS(double T);

    /* Pressure [Pa] on the boundary between regions 2 and 3 (IAPWS-97 Eq. 5),
     * valid for 623.15 K <= T <= 863.15 K. */
    static double PB23IAPWS(double T);

    void printProperties() const;

    // ****************************************************************************************
    //                                  PUBLIC VARIABLES
    // ****************************************************************************************
    //                              (exposed to, e.g., eps_JN)
    //
    double v_;  // specific volume [m^3/kg]
    double u_;  // specific internal energy [J/kg]
    double s_;  // specific entropy [J/kg]
    double h_;  // specific enthalpy [J/kg/K]
    double cp_;  // specific isobaric heat capacity [J/kg]
    double cv_;  // specific isochoric heat capacity [J/kg]
    double w_;  // speed of sound [m/s]
    double g_;  // specific Gibbs free energy [J/kg/K]

    double G_;  // specific Gibbs free energy [J/kg/K]
    double H_;  // specific enthalpy [J/mol/K]
    double denst_;  // 1/v_

    double alpha_;  // Isobaric thermal expansion [1/K]
    double alpha_t_;  //  d(alpha)/dt
    double beta_;  // Isothermal compressibility [1/Pa]
    double Psat_;  //  saturation pressure for given T (NaN above the critical temperature)

    int region_;  // IAPWS-97 region used in the last property evaluation


private:

    double P_;  // Pa
    double T_;  // K
    double R_;  // ideal gas constant kJ/kg/K
    double Mw_;  // mol weight H2O
    double Tcrit_;  // K
    double Pcrit_;  // Pa*
    double rho_crit_;  // kg/m^3

    /* Determines the IAPWS-97 region for (T, P) and calculates Gibbs free energy
    * and thermodynamic constants:
    *
    *   - specific volume [m^3/kg]
    *   - specific internal energy J/kg, specific entropy [J/kg]
    *   - specific enthalpy [J/kg/K]
    *   - specific isobaric heat capacity [J/kg]
    *   - specific isochoric heat capacity [J/kg/K]
    *   - speed of sound [m/s]
    *
    * Units for input pressure and temperature: [P]=Pa, [T]=K.
    * Throws std::domain_error outside the implemented regions.
    */
    void gibbsIAPWSlocal(double T, double P);

    /* Region-specific property evaluations. */
    void region1(double T, double P);
    void region2(double T, double P);
    void region3(double T, double P);

    /* Sets all member properties from the dimensionless Gibbs free energy g and
    * its (total) derivatives with respect to reduced pressure pi = P/pstar and
    * reduced temperature tau. Shared by regions 1 and 2 (Tables 3 and 12 of the
    * IAPWS-97 paper). */
    void setFromGibbs(double T, double P, double pstar, double tau,
                      double g, double g_p, double g_pp,
                      double g_t, double g_tt, double g_pt, double g_ptt);

    /* Dimensionless Helmholtz free energy phi(delta, tau) of region 3
    * (Eq. 28/Table 30 of the IAPWS-97 paper) and its derivatives.
    * delta = rho/rho_crit, tau = Tcrit/T. */
    static void phiRegion3(double delta, double tau,
                           double& phi, double& phi_d, double& phi_dd,
                           double& phi_t, double& phi_tt, double& phi_dt);

    /* Pressure [Pa] from the region-3 Helmholtz formulation, and its density
    * derivative dP/drho at constant T. */
    static double pressureRegion3(double rho, double T, double& dPdrho);

    /* Solves pressureRegion3(rho, T) = P for the density with a
    * bisection-safeguarded Newton iteration. */
    static double region3Density(double T, double P);

    /* Isobaric thermal expansion coefficient in region 3; used for the
    * finite-difference evaluation of alpha_t_. Does not touch member state. */
    static double region3Alpha(double T, double P);

    /* Auxiliary saturated liquid/vapour density correlations [kg/m^3]
    * (Wagner & Pruss auxiliary equations); used only for bracketing the
    * region-3 density iteration. */
    static double satLiquidDensity(double T);
    static double satVapourDensity(double T);

    /*
    * Calculates integer powers x^n, as well as the first and second derivatives.
    * Stores the results in variables passed in by reference.
    */
    static void mypow(double x, int n, double& xn, double& dxn, double& ddxn);

    template<typename real>
    static void nth_power(real x, int n, real& xn, real& dxn, real& ddxn)
    {
        xn = std::pow(x, n);
        const auto dbl_n = static_cast<real>(n);
        dxn = dbl_n*xn / x;
        ddxn = (dbl_n - 1.0)*dxn / x;
    }

};

#endif

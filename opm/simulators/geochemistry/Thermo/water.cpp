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
#include <opm/simulators/geochemistry/Thermo/water.h>

#include <cstdio>
#include <limits>
#include <stdexcept>
#include <string>

namespace {

// Fixed constants of the IAPWS-97 formulation, duplicated here so that the
// static (stateless) region-3 helper functions can use them.
constexpr double R_H2O = 0.461526e3;    // specific gas constant of water [J/kg/K]
constexpr double T_CRIT = 647.096;      // [K]
constexpr double RHO_CRIT = 322.0;      // [kg/m^3]

// IAPWS-95 uses a slightly different specific gas constant than IAPWS-97;
// the reducing constants T_CRIT and RHO_CRIT are the same.
constexpr double R_IAPWS95 = 461.51805;  // [J/kg/K]

std::string rangeError(const char* what, double T, double P)
{
    return std::string("water (IAPWS-97): ") + what
        + " for T=" + std::to_string(T) + " K, P=" + std::to_string(P) + " Pa";
}

} // anonymous namespace

water::water()
: v_(0.0)
, u_(0.0)
, s_(0.0)
, h_(0.0)
, cp_(0.0)
, cv_(0.0)
, w_(0.0)
, g_(0.0)
, G_(0.0)
, H_(0.0)
, denst_(0.0)
, alpha_(0.0)
, alpha_t_(0.0)
, beta_(0.0)
, Psat_(0.0)
, region_(0)
//
, P_(0.0)
, T_(-273.15)
, R_(R_H2O)
, Mw_(18.01528e-3)  // [kg/mol]
, Tcrit_(T_CRIT)
, Pcrit_(22.064e6)
, rho_crit_(RHO_CRIT)
{

}

/* Note: If T==T && P==P, the calculation is already done. Pressure is in Pascals, temperature in Kelvin. */
void water::gibbsIAPWS(double T, double P)
{
    if (T != T_ || P != P_){
        gibbsIAPWSlocal(T, P);
    }
}

void water::gibbsIAPWSlocal(double T, double P)
{
    if (!(T >= 273.15))
    {
        throw std::domain_error(rangeError("temperature below 273.15 K", T, P));
    }
    if (T > 1073.15)
    {
        throw std::domain_error(rangeError("temperature above 1073.15 K (region 5 is not implemented)", T, P));
    }
    if (!(P > 0.0))
    {
        throw std::domain_error(rangeError("non-positive pressure", T, P));
    }
    if (P > 1000.0e6)
    {
        throw std::domain_error(rangeError("pressure above 1000 MPa", T, P));
    }

    Psat_ = (T <= Tcrit_) ? PsatIAPWS(T) : std::numeric_limits<double>::quiet_NaN();

    if (P > 100.0e6)
    {
        // Outside the IAPWS-97 pressure range; use the IAPWS-95 scientific formulation
        regionIAPWS95(T, P);
    }
    else if (T <= 623.15)
    {
        if (P >= Psat_)
        {
            region1(T, P);
        }
        else
        {
            region2(T, P);
        }
    }
    else if (T <= 863.15 && P > PB23IAPWS(T))
    {
        region3(T, P);
    }
    else
    {
        region2(T, P);
    }

    // Cache only after a successful evaluation, so a failed call is retried.
    T_ = T;
    P_ = P;
}

/* Region 1: compressed liquid, 273.15 K <= T <= 623.15 K, Psat <= P <= 100 MPa. */
void water::region1(double T, double P)
{
    static constexpr std::array<int, 34> Ii = { 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2,
                                                   2, 2, 3, 3, 3, 4, 4, 4, 5, 8, 8, 21, 23, 29, 30, 31, 32 };

    static constexpr std::array<int, 34> Ji = { -2, -1, 0, 1, 2, 3, 4, 5, -9, -7, -1, 0, 1, 3, -3, 0, 1,
                                                 3, 17, -4, 0, 6, -5, -2, 10, -8, -11, -6, -29, -31, -38, -39, -40, -41 };

	static constexpr std::array<double, 34> ni =
    {
        1.46329712131670e-01, -8.45481871691140e-01, -3.75636036720400, 3.38551691683850, -9.57919633878720e-01,
        1.57720385132280e-01, -1.66164171995010e-02, 8.12146299835680e-04, 2.83190801238040e-04, -6.07063015658740e-04,
        -1.89900682184190e-02, -3.25297487705050e-02, -2.18417171754140e-02, -5.28383579699300e-05, -4.71843210732670e-04,
        -3.00017807930260e-04, 4.76613939069870e-05, -4.41418453308460e-06, -7.26949962975940e-16, -3.16796448450540e-05,
        -2.82707979853120e-06, -8.52051281201030e-10, -2.24252819080000e-06, -6.51712228956010e-07, -1.43417299379240e-13,
        -4.05169968601170e-07, -1.27343017416410e-09, -1.74248712306340e-10, -6.87621312955310e-19, 1.44783078285210e-20,
        2.63357816627950e-23, -1.19476226400710e-23, 1.82280945814040e-24, -9.35370872924580e-26
    };

	static constexpr double Ps = 16.53e6;  // 16.53 MPa
    static constexpr double Ts = 1386.0;  // K

    const double pi = P / Ps;
    const double pi_diff = 7.1 - pi;
    const double tau = Ts / T;
    const double tau_diff = tau - 1.222;

    double pi_pow = 0.0;
    double dpi_pow = 0.0;
    double ddpi_pow = 0.0;
    double tau_pow = 0.0;
    double dtau_pow = 0.0;
    double ddtau_pow = 0.0;

    double g = 0.0;  // dimensionless Gibbs free energy (in paper: lambda)
    // For the derivatives, we use shorter notation p~pi, and t~tau
    double dg_p = 0.0;
    double dg_pp = 0.0;
    double dg_t = 0.0;
    double dg_tt = 0.0;
    double dg_pt = 0.0;
    double dg_ptt = 0.0;

    // Equation (7) in IAPWS-97 paper ("the basic equation")
    for (std::size_t i=0; i < ni.size(); ++i)
    {
        nth_power(pi_diff, Ii[i], pi_pow, dpi_pow, ddpi_pow);
        nth_power(tau_diff, Ji[i], tau_pow, dtau_pow, ddtau_pow);

		g += ni[i] * pi_pow * tau_pow;
        dg_p += ni[i] * dpi_pow * tau_pow;
		dg_pp += ni[i] * ddpi_pow * tau_pow;
		dg_t += ni[i] * pi_pow * dtau_pow;
		dg_tt += ni[i] * pi_pow * ddtau_pow;
		dg_pt += ni[i] * dpi_pow * dtau_pow;
		dg_ptt += ni[i] * dpi_pow * ddtau_pow;
	}
    // Because the derivative of the kernel is -1 (chain rule, d/dpi)
    dg_p = -dg_p;
    dg_pt = -dg_pt;
    dg_ptt = -dg_ptt;

    setFromGibbs(T, P, Ps, tau, g, dg_p, dg_pp, dg_t, dg_tt, dg_pt, dg_ptt);
    region_ = 1;
}

/* Region 2: steam/superheated vapour (Eqs. 15-17, Tables 10-12 of the IAPWS-97 paper). */
void water::region2(double T, double P)
{
    // Ideal-gas part (Table 10)
    static constexpr std::array<int, 9> J0 = { 0, 1, -5, -4, -3, -2, -1, 2, 3 };

    static constexpr std::array<double, 9> n0 =
    {
        -0.96927686500217e1, 0.10086655968018e2, -0.56087911283020e-2, 0.71452738081455e-1, -0.40710498223928,
        0.14240819171444e1, -0.43839511319450e1, -0.28408632460772, 0.21268463753307e-1
    };

    // Residual part (Table 11)
    static constexpr std::array<int, 43> Ir = { 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 5, 6, 6, 6,
                                                7, 7, 7, 8, 8, 9, 10, 10, 10, 16, 16, 18, 20, 20, 20, 21, 22, 23, 24, 24, 24 };

    static constexpr std::array<int, 43> Jr = { 0, 1, 2, 3, 6, 1, 2, 4, 7, 36, 0, 1, 3, 6, 35, 1, 2, 3, 7, 3, 16, 35,
                                                0, 11, 25, 8, 36, 13, 4, 10, 14, 29, 50, 57, 20, 35, 48, 21, 53, 39, 26, 40, 58 };

    static constexpr std::array<double, 43> nr =
    {
        -0.17731742473213e-2, -0.17834862292358e-1, -0.45996013696365e-1, -0.57581259083432e-1, -0.50325278727930e-1,
        -0.33032641670203e-4, -0.18948987516315e-3, -0.39392777243355e-2, -0.43797295650573e-1, -0.26674547914087e-4,
        0.20481737692309e-7, 0.43870667284435e-6, -0.32277677238570e-4, -0.15033924542148e-2, -0.40668253562649e-1,
        -0.78847309559367e-9, 0.12790717852285e-7, 0.48225372718507e-6, 0.22922076337661e-5, -0.16714766451061e-10,
        -0.21171472321355e-2, -0.23895741934104e2, -0.59059564324270e-17, -0.12621808899101e-5, -0.38946842435739e-1,
        0.11256211360459e-10, -0.82311340897998e1, 0.19809712802088e-7, 0.10406965210174e-18, -0.10234747095929e-12,
        -0.10018179379511e-8, -0.80882908646985e-10, 0.10693031879409, -0.33662250574171, 0.89185845355421e-24,
        0.30629316876232e-12, -0.42002467698208e-5, -0.59056029685639e-25, 0.37826947613457e-5, -0.12768608934681e-14,
        0.73087610595061e-28, 0.55414715350778e-16, -0.94369707241210e-6
    };

    static constexpr double Ps = 1.0e6;  // 1 MPa
    static constexpr double Ts = 540.0;  // K

    const double pi = P / Ps;
    const double tau = Ts / T;
    const double tau_diff = tau - 0.5;

    // Ideal-gas part: gamma0 = ln(pi) + sum n0_i tau^J0_i (Eq. 16)
    double g = std::log(pi);
    double dg_p = 1.0 / pi;
    double dg_pp = -1.0 / (pi*pi);
    double dg_t = 0.0;
    double dg_tt = 0.0;
    double dg_pt = 0.0;
    double dg_ptt = 0.0;

    double pi_pow = 0.0;
    double dpi_pow = 0.0;
    double ddpi_pow = 0.0;
    double tau_pow = 0.0;
    double dtau_pow = 0.0;
    double ddtau_pow = 0.0;

    for (std::size_t i=0; i < n0.size(); ++i)
    {
        nth_power(tau, J0[i], tau_pow, dtau_pow, ddtau_pow);

        g += n0[i] * tau_pow;
        dg_t += n0[i] * dtau_pow;
        dg_tt += n0[i] * ddtau_pow;
    }

    // Residual part: gammar = sum nr_i pi^Ir_i (tau - 0.5)^Jr_i (Eq. 17)
    for (std::size_t i=0; i < nr.size(); ++i)
    {
        nth_power(pi, Ir[i], pi_pow, dpi_pow, ddpi_pow);
        nth_power(tau_diff, Jr[i], tau_pow, dtau_pow, ddtau_pow);

        g += nr[i] * pi_pow * tau_pow;
        dg_p += nr[i] * dpi_pow * tau_pow;
        dg_pp += nr[i] * ddpi_pow * tau_pow;
        dg_t += nr[i] * pi_pow * dtau_pow;
        dg_tt += nr[i] * pi_pow * ddtau_pow;
        dg_pt += nr[i] * dpi_pow * dtau_pow;
        dg_ptt += nr[i] * dpi_pow * ddtau_pow;
    }

    setFromGibbs(T, P, Ps, tau, g, dg_p, dg_pp, dg_t, dg_tt, dg_pt, dg_ptt);
    region_ = 2;
}

/* Region 3: near-critical/supercritical (Eq. 28, Tables 30-31 of the IAPWS-97 paper).
 * The formulation is in Helmholtz free energy phi(rho, T), so the density is first
 * obtained from P with a safeguarded Newton iteration. */
void water::region3(double T, double P)
{
    const double rho = region3Density(T, P);
    const double delta = rho / rho_crit_;
    const double tau = Tcrit_ / T;

    double phi = 0.0;
    double phi_d = 0.0;
    double phi_dd = 0.0;
    double phi_t = 0.0;
    double phi_tt = 0.0;
    double phi_dt = 0.0;
    phiRegion3(delta, tau, phi, phi_d, phi_dd, phi_t, phi_tt, phi_dt);

    setFromHelmholtz(T, rho, R_, phi, phi_d, phi_dd, phi_t, phi_tt, phi_dt);

    // d(alpha)/dT along the isobar by central differences; an analytical
    // expression would require third derivatives of phi. The probes are
    // bracketed around the converged density so that they stay on the same
    // branch of the two-phase dome as the base state.
    static constexpr double dT = 1.0e-2;
    alpha_t_ = (region3Alpha(T + dT, P, rho) - region3Alpha(T - dT, P, rho)) / (2.0*dT);

    region_ = 3;
}

/* IAPWS-95 (Wagner & Pruss 2002), used for 100 MPa < P <= 1000 MPa where the
 * industrial formulation does not apply. Helmholtz-based, like region 3. */
void water::regionIAPWS95(double T, double P)
{
    const double rho = densityIAPWS95(T, P);
    const double delta = rho / rho_crit_;
    const double tau = Tcrit_ / T;

    double phi = 0.0;
    double phi_d = 0.0;
    double phi_dd = 0.0;
    double phi_t = 0.0;
    double phi_tt = 0.0;
    double phi_dt = 0.0;
    phiIAPWS95(delta, tau, phi, phi_d, phi_dd, phi_t, phi_tt, phi_dt);

    setFromHelmholtz(T, rho, R_IAPWS95, phi, phi_d, phi_dd, phi_t, phi_tt, phi_dt);

    static constexpr double dT = 1.0e-2;
    alpha_t_ = (alphaIAPWS95(T + dT, P, rho) - alphaIAPWS95(T - dT, P, rho)) / (2.0*dT);

    region_ = 95;
}

/* Sets the member properties from the dimensionless Gibbs free energy and its
 * total derivatives with respect to pi = P/pstar and tau (regions 1 and 2). */
void water::setFromGibbs(double T, double P, double pstar, double tau,
                         double g, double g_p, double g_pp,
                         double g_t, double g_tt, double g_pt, double g_ptt)
{
    const double pi = P / pstar;

    // Table 3 of IAPWS-97 paper (identical structure for Table 12)
	v_ = pi*g_p*R_*T / P; // m^3/kg
	u_ = R_*T*(tau*g_t - pi*g_p);
	s_ = R_*(tau*g_t - g);
	h_ = R_*T*tau*g_t;
	cp_ = -tau*tau*g_tt*R_;
	const double cvi = (g_p - tau*g_pt);
	cv_ = cp_ + R_*cvi *cvi/ g_pp;
	w_ = R_*T*g_p*g_p / (cvi*cvi / tau / tau / g_tt - g_pp);
	w_ = sqrt(w_);

    g_ = g*R_*T;
    G_ = g_*Mw_;
    H_ = h_*Mw_;
    denst_ = 1.0 / v_;

    // beta_  = water isothermal compressibility = -1/v_(dv_/dp)_T
    // alpha_ = water isobaric thermal expansion = 1/v_(dv_/dT)_p
    // calculated analytically by replacing specific volume with v_=RT/pstar gamma_pi in
    // Table 3 of IAPWS-97 paper
	alpha_ = (g_p - tau*g_pt) / T / g_p;
	beta_ = -g_pp / pstar / g_p;
	alpha_t_ = (g_p - tau*g_pt);
	alpha_t_ *= alpha_t_;
	alpha_t_ = tau*tau*g_ptt*g_p - alpha_t_;
	alpha_t_ = alpha_t_ / (T*T*g_p*g_p);
}

/* Sets the member properties (except alpha_t_) from the dimensionless Helmholtz
 * free energy and its total derivatives with respect to delta and tau. Shared by
 * IAPWS-97 region 3 (Table 31) and IAPWS-95 (Table 3 of the 1995 release; the
 * formulas coincide when expressed in the total phi = phi0 + phir). */
void water::setFromHelmholtz(double T, double rho, double R,
                             double phi, double phi_d, double phi_dd,
                             double phi_t, double phi_tt, double phi_dt)
{
    const double delta = rho / rho_crit_;
    const double tau = Tcrit_ / T;

    v_ = 1.0 / rho;
    u_ = R*T*tau*phi_t;
    s_ = R*(tau*phi_t - phi);
    h_ = R*T*(tau*phi_t + delta*phi_d);
    cv_ = -R*tau*tau*phi_tt;
    const double num = delta*phi_d - delta*tau*phi_dt;
    const double den = 2.0*delta*phi_d + delta*delta*phi_dd;
    cp_ = cv_ + R*num*num / den;
    w_ = sqrt(R*T*(den - num*num / (tau*tau*phi_tt)));

    g_ = R*T*(phi + delta*phi_d);  // g = f + P/rho
    G_ = g_*Mw_;
    H_ = h_*Mw_;
    denst_ = rho;

    // alpha and beta from the P(rho, T) partial derivatives:
    //   beta  = 1/(rho (dP/drho)_T),  alpha = beta (dP/dT)_rho
    const double dPdT = rho*R*delta*(phi_d - tau*phi_dt);
    const double dPdrho = R*T*den;
    beta_ = 1.0 / (rho*dPdrho);
    alpha_ = beta_*dPdT;
}

/* Dimensionless Helmholtz free energy of region 3 (Eq. 28/Table 30) and its derivatives. */
void water::phiRegion3(double delta, double tau,
                       double& phi, double& phi_d, double& phi_dd,
                       double& phi_t, double& phi_tt, double& phi_dt)
{
    static constexpr double n1 = 0.10658070028513e1;

    static constexpr std::array<int, 39> Ii = { 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3,
                                                3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 8, 9, 9, 10, 10, 11 };

    static constexpr std::array<int, 39> Ji = { 0, 1, 2, 7, 10, 12, 23, 2, 6, 15, 17, 0, 2, 6, 7, 22, 26, 0, 2, 4,
                                                16, 26, 0, 2, 4, 26, 1, 3, 26, 0, 2, 26, 2, 26, 2, 26, 0, 1, 26 };

    static constexpr std::array<double, 39> ni =
    {
        -0.15732845290239e2, 0.20944396974307e2, -0.76867707878716e1, 0.26185947787954e1, -0.28080781148620e1,
        0.12053369696517e1, -0.84566812812502e-2, -0.12654315477714e1, -0.11524407806681e1, 0.88521043984318,
        -0.64207765181607, 0.38493460186671, -0.85214708824206, 0.48972281541877e1, -0.30502617256965e1,
        0.39420536879154e-1, 0.12558408424308, -0.27999329698710, 0.13899799569460e1, -0.20189915023570e1,
        -0.82147637173963e-2, -0.47596035734923, 0.43984074473500e-1, -0.44476435428739, 0.90572070719733,
        0.70522450087967, 0.10770512626332, -0.32913623258954, -0.50871062041158, -0.22175400873096e-1,
        0.94260751665092e-1, 0.16436278447961, -0.13503372241348e-1, -0.14834345352472e-1, 0.57922953628084e-3,
        0.32308904703711e-2, 0.80964802996215e-4, -0.16557679795037e-3, -0.44923899061815e-4
    };

    phi = n1*std::log(delta);
    phi_d = n1 / delta;
    phi_dd = -n1 / (delta*delta);
    phi_t = 0.0;
    phi_tt = 0.0;
    phi_dt = 0.0;

    double d_pow = 0.0;
    double dd_pow = 0.0;
    double ddd_pow = 0.0;
    double t_pow = 0.0;
    double dt_pow = 0.0;
    double ddt_pow = 0.0;

    for (std::size_t i=0; i < ni.size(); ++i)
    {
        nth_power(delta, Ii[i], d_pow, dd_pow, ddd_pow);
        nth_power(tau, Ji[i], t_pow, dt_pow, ddt_pow);

        phi += ni[i] * d_pow * t_pow;
        phi_d += ni[i] * dd_pow * t_pow;
        phi_dd += ni[i] * ddd_pow * t_pow;
        phi_t += ni[i] * d_pow * dt_pow;
        phi_tt += ni[i] * d_pow * ddt_pow;
        phi_dt += ni[i] * dd_pow * dt_pow;
    }
}

/* Pressure in region 3: P = rho R T delta phi_delta, plus dP/drho at constant T. */
double water::pressureRegion3(double rho, double T, double& dPdrho)
{
    const double delta = rho / RHO_CRIT;
    const double tau = T_CRIT / T;

    double phi = 0.0;
    double phi_d = 0.0;
    double phi_dd = 0.0;
    double phi_t = 0.0;
    double phi_tt = 0.0;
    double phi_dt = 0.0;
    phiRegion3(delta, tau, phi, phi_d, phi_dd, phi_t, phi_tt, phi_dt);

    dPdrho = R_H2O*T*(2.0*delta*phi_d + delta*delta*phi_dd);
    return rho*R_H2O*T*delta*phi_d;
}

/* Solves pressureRegion3(rho, T) = P for rho. Below the critical temperature the
 * initial bracket is placed on the correct side of the two-phase dome using the
 * auxiliary saturated density correlations. */
double water::region3Density(double T, double P)
{
    double lo = 1.0e-3;
    double hi = 800.0;

    if (T <= T_CRIT)
    {
        if (P >= PsatIAPWS(T))
        {
            lo = 0.98*satLiquidDensity(T);  // liquid branch
        }
        else
        {
            hi = 1.02*satVapourDensity(T);  // vapour branch
        }
    }

    return solveDensity(&water::pressureRegion3, T, P, lo, hi);
}

/* Bisection-safeguarded Newton iteration solving pfn(rho, T) = P for rho. */
double water::solveDensity(PressureFn pfn, double T, double P, double lo, double hi)
{
    double dPdrho = 0.0;
    double flo = pfn(lo, T, dPdrho) - P;
    double fhi = pfn(hi, T, dPdrho) - P;

    // The initial brackets (e.g. from the auxiliary saturated-density
    // correlations) can be slightly off; nudge them until the root is enclosed.
    for (int k=0; k < 400 && flo > 0.0; ++k)
    {
        lo *= 0.99;
        flo = pfn(lo, T, dPdrho) - P;
    }
    for (int k=0; k < 400 && fhi < 0.0; ++k)
    {
        hi *= 1.01;
        fhi = pfn(hi, T, dPdrho) - P;
    }
    if (flo > 0.0 || fhi < 0.0)
    {
        throw std::runtime_error(rangeError("could not bracket the density", T, P));
    }

    double rho = 0.5*(lo + hi);
    for (int it=0; it < 200; ++it)
    {
        const double f = pfn(rho, T, dPdrho) - P;
        if (std::fabs(f) <= 1.0e-10*P)
        {
            return rho;
        }

        if (f > 0.0)
        {
            hi = rho;
        }
        else
        {
            lo = rho;
        }

        // Newton step, bisection whenever the step is invalid or leaves the bracket
        double next = (dPdrho > 0.0) ? rho - f/dPdrho : 0.0;
        if (!(next > lo && next < hi))
        {
            next = 0.5*(lo + hi);
        }
        if (hi - lo < 1.0e-13*hi)
        {
            return rho;  // bracket exhausted (essentially flat isotherm)
        }
        rho = next;
    }

    throw std::runtime_error(rangeError("the density iteration did not converge", T, P));
}

/* Isobaric thermal expansion coefficient in region 3 (no member state is touched).
 * The bracket around rho_guess keeps the probe on the branch of the base state
 * and warm-starts the solve; solveDensity widens the bracket if needed. */
double water::region3Alpha(double T, double P, double rho_guess)
{
    const double rho = solveDensity(&water::pressureRegion3, T, P,
                                    0.9*rho_guess, 1.1*rho_guess);
    const double delta = rho / RHO_CRIT;
    const double tau = T_CRIT / T;

    double phi = 0.0;
    double phi_d = 0.0;
    double phi_dd = 0.0;
    double phi_t = 0.0;
    double phi_tt = 0.0;
    double phi_dt = 0.0;
    phiRegion3(delta, tau, phi, phi_d, phi_dd, phi_t, phi_tt, phi_dt);

    const double dPdT = rho*R_H2O*delta*(phi_d - tau*phi_dt);
    const double dPdrho = R_H2O*T*(2.0*delta*phi_d + delta*delta*phi_dd);

    return dPdT / (rho*dPdrho);
}

/* Dimensionless Helmholtz free energy of IAPWS-95, phi = phi0 + phir, and its
 * derivatives (Wagner & Pruss 2002, Tables 6.1/6.2 or the IAPWS-95 release,
 * Tables 1/2 with derivatives per Tables 4/5). */
void water::phiIAPWS95(double delta, double tau,
                       double& phi, double& phi_d, double& phi_dd,
                       double& phi_t, double& phi_tt, double& phi_dt)
{
    // ---------------------------------------------------------------- ideal part
    static constexpr std::array<double, 8> n0 =
    {
        -8.3204464837497, 6.6832105275932, 3.00632, 0.012436, 0.97315, 1.27950, 0.96956, 0.24873
    };
    static constexpr std::array<double, 5> gamma0 = { 1.28728967, 3.53734222, 7.74073708, 9.24437796, 27.5075105 };

    phi = std::log(delta) + n0[0] + n0[1]*tau + n0[2]*std::log(tau);
    phi_d = 1.0 / delta;
    phi_dd = -1.0 / (delta*delta);
    phi_t = n0[1] + n0[2] / tau;
    phi_tt = -n0[2] / (tau*tau);
    phi_dt = 0.0;

    for (std::size_t i=0; i < gamma0.size(); ++i)
    {
        const double e = std::exp(-gamma0[i]*tau);
        phi += n0[3+i]*std::log(1.0 - e);
        phi_t += n0[3+i]*gamma0[i]*(1.0/(1.0 - e) - 1.0);
        phi_tt -= n0[3+i]*gamma0[i]*gamma0[i]*e / ((1.0 - e)*(1.0 - e));
    }

    // ------------------------------------------------- residual: polynomial terms
    static constexpr std::array<int, 7> d1 = { 1, 1, 1, 2, 2, 3, 4 };
    static constexpr std::array<double, 7> t1 = { -0.5, 0.875, 1.0, 0.5, 0.75, 0.375, 1.0 };
    static constexpr std::array<double, 7> n1 =
    {
        0.12533547935523e-1, 0.78957634722828e1, -0.87803203303561e1, 0.31802509345418,
        -0.26145533859358, -0.78199751687981e-2, 0.88089493102134e-2
    };

    for (std::size_t i=0; i < n1.size(); ++i)
    {
        const double dd = static_cast<double>(d1[i]);
        const double dp = std::pow(delta, dd);
        const double tp = std::pow(tau, t1[i]);
        const double f = n1[i]*dp*tp;

        phi += f;
        phi_d += f*dd/delta;
        phi_dd += f*dd*(dd - 1.0)/(delta*delta);
        phi_t += f*t1[i]/tau;
        phi_tt += f*t1[i]*(t1[i] - 1.0)/(tau*tau);
        phi_dt += f*dd*t1[i]/(delta*tau);
    }

    // ------------------------------------------------ residual: exponential terms
    static constexpr std::array<int, 44> c2 = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                                2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                                                3, 3, 3, 3, 4, 6, 6, 6, 6 };
    static constexpr std::array<int, 44> d2 = { 1, 1, 1, 2, 2, 3, 4, 4, 5, 7, 9, 10, 11, 13, 15,
                                                1, 2, 2, 2, 3, 4, 4, 4, 5, 6, 6, 7, 9, 9, 9, 9, 9, 10, 10, 12,
                                                3, 4, 4, 5, 14, 3, 6, 6, 6 };
    static constexpr std::array<int, 44> t2 = { 4, 6, 12, 1, 5, 4, 2, 13, 9, 3, 4, 11, 4, 13, 1,
                                                7, 1, 9, 10, 10, 3, 7, 10, 10, 6, 10, 10, 1, 2, 3, 4, 8, 6, 9, 8,
                                                16, 22, 23, 23, 10, 50, 44, 46, 50 };
    static constexpr std::array<double, 44> n2 =
    {
        -0.66856572307965, 0.20433810950965, -0.66212605039687e-4, -0.19232721156002, -0.25709043003438,
        0.16074868486251, -0.40092828925807e-1, 0.39343422603254e-6, -0.75941377088144e-5, 0.56250979351888e-3,
        -0.15608652257135e-4, 0.11537996422951e-8, 0.36582165144204e-6, -0.13251180074668e-11, -0.62639586912454e-9,
        -0.10793600908932, 0.17611491008752e-1, 0.22132295167546, -0.40247669763528, 0.58083399985759,
        0.49969146990806e-2, -0.31358700712549e-1, -0.74315929710341, 0.47807329915480, 0.20527940895948e-1,
        -0.13636435110343, 0.14180634400617e-1, 0.83326504880713e-2, -0.29052336009585e-1, 0.38615085574206e-1,
        -0.20393486513704e-1, -0.16554050063734e-2, 0.19955571979541e-2, 0.15870308324157e-3, -0.16388568342530e-4,
        0.43613615723811e-1, 0.34994005463765e-1, -0.76788197844621e-1, 0.22446277332006e-1, -0.62689710414685e-4,
        -0.55711118565645e-9, -0.19905718354408, 0.31777497330738, -0.11841182425981
    };

    for (std::size_t i=0; i < n2.size(); ++i)
    {
        const double cc = static_cast<double>(c2[i]);
        const double dd = static_cast<double>(d2[i]);
        const double tt = static_cast<double>(t2[i]);
        const double dc = std::pow(delta, cc);
        const double f = n2[i]*std::pow(delta, dd)*std::pow(tau, tt)*std::exp(-dc);

        phi += f;
        phi_d += f*(dd - cc*dc)/delta;
        phi_dd += f*((dd - cc*dc)*(dd - 1.0 - cc*dc) - cc*cc*dc)/(delta*delta);
        phi_t += f*tt/tau;
        phi_tt += f*tt*(tt - 1.0)/(tau*tau);
        phi_dt += f*(dd - cc*dc)*tt/(delta*tau);
    }

    // --------------------------------------------------- residual: Gaussian terms
    static constexpr std::array<int, 3> d3 = { 3, 3, 3 };
    static constexpr std::array<int, 3> t3 = { 0, 1, 4 };
    static constexpr std::array<double, 3> n3 = { -0.31306260323435e2, 0.31546140237781e2, -0.25213154341695e4 };
    static constexpr std::array<double, 3> alpha3 = { 20.0, 20.0, 20.0 };
    static constexpr std::array<double, 3> beta3 = { 150.0, 150.0, 250.0 };
    static constexpr std::array<double, 3> gamma3 = { 1.21, 1.21, 1.25 };
    static constexpr std::array<double, 3> eps3 = { 1.0, 1.0, 1.0 };

    for (std::size_t i=0; i < n3.size(); ++i)
    {
        const double dd = static_cast<double>(d3[i]);
        const double tt = static_cast<double>(t3[i]);
        const double f = n3[i]*std::pow(delta, dd)*std::pow(tau, tt)
            *std::exp(-alpha3[i]*(delta - eps3[i])*(delta - eps3[i])
                      - beta3[i]*(tau - gamma3[i])*(tau - gamma3[i]));
        if (f == 0.0)
        {
            continue;
        }

        const double ad = dd/delta - 2.0*alpha3[i]*(delta - eps3[i]);
        const double at = tt/tau - 2.0*beta3[i]*(tau - gamma3[i]);

        phi += f;
        phi_d += f*ad;
        phi_dd += f*(ad*ad - dd/(delta*delta) - 2.0*alpha3[i]);
        phi_t += f*at;
        phi_tt += f*(at*at - tt/(tau*tau) - 2.0*beta3[i]);
        phi_dt += f*ad*at;
    }

    // ----------------------------------------------- residual: nonanalytic terms
    static constexpr std::array<double, 2> a4 = { 3.5, 3.5 };
    static constexpr std::array<double, 2> b4 = { 0.85, 0.95 };
    static constexpr std::array<double, 2> B4 = { 0.2, 0.2 };
    static constexpr std::array<double, 2> n4 = { -0.14874640856724, 0.31861088019884 };
    static constexpr std::array<double, 2> C4 = { 28.0, 32.0 };
    static constexpr std::array<double, 2> D4 = { 700.0, 800.0 };
    static constexpr std::array<double, 2> A4 = { 0.32, 0.32 };
    static constexpr std::array<double, 2> beta4 = { 0.3, 0.3 };

    // Nudge away from the removable singularity of the distance function at the
    // critical density.
    const double dm1 = (std::fabs(delta - 1.0) < 1.0e-10)
        ? ((delta >= 1.0) ? 1.0e-10 : -1.0e-10)
        : (delta - 1.0);
    const double dm1sq = dm1*dm1;

    for (std::size_t i=0; i < n4.size(); ++i)
    {
        const double psi = std::exp(-C4[i]*dm1sq - D4[i]*(tau - 1.0)*(tau - 1.0));
        if (psi == 0.0)
        {
            // The term (and all its derivatives) vanishes; skipping also avoids
            // overflow in the negative powers of (delta-1)^2 far from the
            // critical point.
            continue;
        }

        const double inv2beta = 1.0/(2.0*beta4[i]);
        const double theta = (1.0 - tau) + A4[i]*std::pow(dm1sq, inv2beta);
        const double Delta = theta*theta + B4[i]*std::pow(dm1sq, a4[i]);
        const double Db = std::pow(Delta, b4[i]);

        const double dDelta_dd = dm1*(A4[i]*theta*(2.0/beta4[i])*std::pow(dm1sq, inv2beta - 1.0)
                                      + 2.0*B4[i]*a4[i]*std::pow(dm1sq, a4[i] - 1.0));
        const double d2Delta_dd2 = dDelta_dd/dm1
            + dm1sq*(4.0*B4[i]*a4[i]*(a4[i] - 1.0)*std::pow(dm1sq, a4[i] - 2.0)
                     + 2.0*A4[i]*A4[i]/(beta4[i]*beta4[i])
                       *std::pow(dm1sq, inv2beta - 1.0)*std::pow(dm1sq, inv2beta - 1.0)
                     + A4[i]*theta*(4.0/beta4[i])*(inv2beta - 1.0)*std::pow(dm1sq, inv2beta - 2.0));

        const double dDb_dd = b4[i]*std::pow(Delta, b4[i] - 1.0)*dDelta_dd;
        const double d2Db_dd2 = b4[i]*(std::pow(Delta, b4[i] - 1.0)*d2Delta_dd2
                                       + (b4[i] - 1.0)*std::pow(Delta, b4[i] - 2.0)*dDelta_dd*dDelta_dd);
        const double dDb_dt = -2.0*theta*b4[i]*std::pow(Delta, b4[i] - 1.0);
        const double d2Db_dt2 = 2.0*b4[i]*std::pow(Delta, b4[i] - 1.0)
            + 4.0*theta*theta*b4[i]*(b4[i] - 1.0)*std::pow(Delta, b4[i] - 2.0);
        const double d2Db_ddt = -A4[i]*b4[i]*(2.0/beta4[i])*std::pow(Delta, b4[i] - 1.0)*dm1
                                  *std::pow(dm1sq, inv2beta - 1.0)
            - 2.0*theta*b4[i]*(b4[i] - 1.0)*std::pow(Delta, b4[i] - 2.0)*dDelta_dd;

        const double dpsi_dd = -2.0*C4[i]*dm1*psi;
        const double d2psi_dd2 = (2.0*C4[i]*dm1sq - 1.0)*2.0*C4[i]*psi;
        const double dpsi_dt = -2.0*D4[i]*(tau - 1.0)*psi;
        const double d2psi_dt2 = (2.0*D4[i]*(tau - 1.0)*(tau - 1.0) - 1.0)*2.0*D4[i]*psi;
        const double d2psi_ddt = 4.0*C4[i]*D4[i]*dm1*(tau - 1.0)*psi;

        phi += n4[i]*Db*delta*psi;
        phi_d += n4[i]*(Db*(psi + delta*dpsi_dd) + dDb_dd*delta*psi);
        phi_dd += n4[i]*(Db*(2.0*dpsi_dd + delta*d2psi_dd2)
                         + 2.0*dDb_dd*(psi + delta*dpsi_dd) + d2Db_dd2*delta*psi);
        phi_t += n4[i]*delta*(dDb_dt*psi + Db*dpsi_dt);
        phi_tt += n4[i]*delta*(d2Db_dt2*psi + 2.0*dDb_dt*dpsi_dt + Db*d2psi_dt2);
        phi_dt += n4[i]*(Db*(dpsi_dt + delta*d2psi_ddt) + delta*dDb_dd*dpsi_dt
                         + dDb_dt*(psi + delta*dpsi_dd) + d2Db_ddt*delta*psi);
    }
}

/* Pressure from IAPWS-95: P = rho R T delta phi_delta, plus dP/drho at constant T. */
double water::pressureIAPWS95(double rho, double T, double& dPdrho)
{
    const double delta = rho / RHO_CRIT;
    const double tau = T_CRIT / T;

    double phi = 0.0;
    double phi_d = 0.0;
    double phi_dd = 0.0;
    double phi_t = 0.0;
    double phi_tt = 0.0;
    double phi_dt = 0.0;
    phiIAPWS95(delta, tau, phi, phi_d, phi_dd, phi_t, phi_tt, phi_dt);

    dPdrho = R_IAPWS95*T*(2.0*delta*phi_d + delta*delta*phi_dd);
    return rho*R_IAPWS95*T*delta*phi_d;
}

/* Solves pressureIAPWS95(rho, T) = P for rho. Only used for P > 100 MPa, far
 * above the critical pressure, where the isotherms have a single (liquid-like
 * or supercritical) root. */
double water::densityIAPWS95(double T, double P)
{
    return solveDensity(&water::pressureIAPWS95, T, P, 1.0, 2000.0);
}

/* Isobaric thermal expansion coefficient from IAPWS-95 (no member state is touched).
 * The bracket around rho_guess warm-starts the solve; solveDensity widens the
 * bracket if needed. */
double water::alphaIAPWS95(double T, double P, double rho_guess)
{
    const double rho = solveDensity(&water::pressureIAPWS95, T, P,
                                    0.9*rho_guess, 1.1*rho_guess);
    const double delta = rho / RHO_CRIT;
    const double tau = T_CRIT / T;

    double phi = 0.0;
    double phi_d = 0.0;
    double phi_dd = 0.0;
    double phi_t = 0.0;
    double phi_tt = 0.0;
    double phi_dt = 0.0;
    phiIAPWS95(delta, tau, phi, phi_d, phi_dd, phi_t, phi_tt, phi_dt);

    const double dPdT = rho*R_IAPWS95*delta*(phi_d - tau*phi_dt);
    const double dPdrho = R_IAPWS95*T*(2.0*delta*phi_d + delta*delta*phi_dd);

    return dPdT / (rho*dPdrho);
}

/* Saturated liquid density [kg/m^3], auxiliary correlation of Wagner & Pruss (2002). */
double water::satLiquidDensity(double T)
{
    static constexpr std::array<double, 6> b = { 1.99274064, 1.09965342, -0.510839303,
                                                 -1.75493479, -45.5170352, -6.74694450e5 };

    const double th = std::max(1.0 - T/T_CRIT, 0.0);

    return RHO_CRIT*(1.0
        + b[0]*std::pow(th, 1.0/3.0)
        + b[1]*std::pow(th, 2.0/3.0)
        + b[2]*std::pow(th, 5.0/3.0)
        + b[3]*std::pow(th, 16.0/3.0)
        + b[4]*std::pow(th, 43.0/3.0)
        + b[5]*std::pow(th, 110.0/3.0));
}

/* Saturated vapour density [kg/m^3], auxiliary correlation of Wagner & Pruss (2002). */
double water::satVapourDensity(double T)
{
    static constexpr std::array<double, 6> c = { -2.03150240, -2.68302940, -5.38626492,
                                                 -17.2991605, -44.7586581, -63.9201063 };

    const double th = std::max(1.0 - T/T_CRIT, 0.0);

    return RHO_CRIT*std::exp(
          c[0]*std::pow(th, 2.0/6.0)
        + c[1]*std::pow(th, 4.0/6.0)
        + c[2]*std::pow(th, 8.0/6.0)
        + c[3]*std::pow(th, 18.0/6.0)
        + c[4]*std::pow(th, 37.0/6.0)
        + c[5]*std::pow(th, 71.0/6.0));
}

/* Returns the saturation pressure for 273.15 K <= T <= 647.096 K. */
double water::PsatIAPWS(double T)
{
    if (T < 273.15 || T > 647.096)
    {
        throw std::domain_error("water (IAPWS-97): temperature T=" + std::to_string(T)
            + " K outside the saturation-line range 273.15 K to 647.096 K");
    }

    static constexpr std::array<double, 10> n = { 0.11670521452767e4, -0.72421316703206e6, -0.17073846940092e2, 0.12020824702470e5,
                                             -0.32325550322333e7, 0.14915108613530e2, -0.48232657361591e4, 0.40511340542057e6,
                                             -0.23855557567849, 0.65017534844798e3 };
    static constexpr double Ts = 1.0;

    const double t = T / Ts;
    const double Th = t + n[8] / (t - n[9]);
    const double Th2 = Th*Th;
    const double A = Th2 + n[0] * Th + n[1];
    const double B = n[2] * Th2 + n[3] * Th + n[4];
    const double C = n[5] * Th2 + n[6] * Th + n[7];

    double p = -B + sqrt(B*B - 4 * A*C);
    p = 2.0*C / p;

    return 1.0e6*p*p*p*p;
}

/* Pressure on the region 2/3 boundary (Eq. 5 of the IAPWS-97 paper). */
double water::PB23IAPWS(double T)
{
    if (T < 623.15 || T > 863.15)
    {
        throw std::domain_error("water (IAPWS-97): temperature T=" + std::to_string(T)
            + " K outside the B23 boundary range 623.15 K to 863.15 K");
    }

    static constexpr double n1 = 0.34805185628969e3;
    static constexpr double n2 = -0.11671859879975e1;
    static constexpr double n3 = 0.10192970039326e-2;

    return 1.0e6*(n1 + n2*T + n3*T*T);
}

void water::printProperties() const
{
    printf("T[K]\tP[MPa]\tv[kg/m3]\th[kJ/kg]\tu[kJ/kg]\ts[kJ/kgK]\tcp[kJ/kgK]\tcv[kJ/kgK]\tw[m/s]\trho[kg/m3]\n");
    printf
    (
        //"{:4.8e}\t{:4.8e}\t{:4.8e}\t{:4.8e}\t{:4.8e}\t{:4.8e}\t{:4.8e}\t{:4.8e}\t{:4.8e}\t{:4.8e}\n",
        "%4.8e\t%4.8e\t%4.8e\t%4.8e\t%4.8e\t%4.8e\t%4.8e\t%4.8e\t%4.8e\t%4.8e\n",
        T_,
        1.0e-6*P_,
        v_,
        1.0e-3*h_,
        1.0e-3*u_,
        1.0e-3*s_,
        1.0e-3*cp_,
        1.0e-3*cv_,
        w_,
        denst_
    );
}

void water::mypow(double x, int n, double& xn, double& dxn, double& ddxn)
{
    if (n == 0)
    {
	xn = 1.0;
        dxn = 0.0;
        ddxn = 0.0;
    }
    else if (n == 1)
    {
        xn = x;
        dxn = 1.0;
        ddxn = 0.0;
    }
    else if (n == 2)
    {
        xn = x*x;
        dxn = 2.0*x;
        ddxn = 2.0;
    }
    else if (xn == 0.0)
    {
        xn = dxn = ddxn = 0.0;
    }
    else
    {
        const bool positive_power = (n>0);
        const int p = positive_power ? n : -n;
        double xp = 1.0;
        for (int i = 0; i < p; ++i)
        {
            xp *= x;
        }

        const auto dbl_n = static_cast<double>(n);
        xn = positive_power ? xp : 1.0 / xp;

        dxn = dbl_n*xn / x;
        ddxn = (dbl_n - 1.0)*dxn / x;
    }
}

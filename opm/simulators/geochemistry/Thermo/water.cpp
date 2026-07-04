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
    if (P > 100.0e6)
    {
        throw std::domain_error(rangeError("pressure above 100 MPa", T, P));
    }

    Psat_ = (T <= Tcrit_) ? PsatIAPWS(T) : std::numeric_limits<double>::quiet_NaN();

    if (T <= 623.15)
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

    // Table 31 of IAPWS-97 paper
    v_ = 1.0 / rho;
    u_ = R_*T*tau*phi_t;
    s_ = R_*(tau*phi_t - phi);
    h_ = R_*T*(tau*phi_t + delta*phi_d);
    cv_ = -R_*tau*tau*phi_tt;
    const double num = delta*phi_d - delta*tau*phi_dt;
    const double den = 2.0*delta*phi_d + delta*delta*phi_dd;
    cp_ = cv_ + R_*num*num / den;
    w_ = sqrt(R_*T*(den - num*num / (tau*tau*phi_tt)));

    g_ = R_*T*(phi + delta*phi_d);  // g = f + P/rho
    G_ = g_*Mw_;
    H_ = h_*Mw_;
    denst_ = rho;

    // alpha and beta from the P(rho, T) partial derivatives:
    //   beta  = 1/(rho (dP/drho)_T),  alpha = beta (dP/dT)_rho
    const double dPdT = rho*R_*delta*(phi_d - tau*phi_dt);
    const double dPdrho = R_*T*den;
    beta_ = 1.0 / (rho*dPdrho);
    alpha_ = beta_*dPdT;

    // d(alpha)/dT along the isobar by central differences; an analytical
    // expression would require third derivatives of phi.
    static constexpr double dT = 1.0e-2;
    alpha_t_ = (region3Alpha(T + dT, P) - region3Alpha(T - dT, P)) / (2.0*dT);

    region_ = 3;
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

    double dPdrho = 0.0;
    double flo = pressureRegion3(lo, T, dPdrho) - P;
    double fhi = pressureRegion3(hi, T, dPdrho) - P;

    // The brackets from the auxiliary correlations can be slightly off close to
    // the saturation line; nudge them until the root is enclosed.
    for (int k=0; k < 400 && flo > 0.0; ++k)
    {
        lo *= 0.99;
        flo = pressureRegion3(lo, T, dPdrho) - P;
    }
    for (int k=0; k < 400 && fhi < 0.0; ++k)
    {
        hi *= 1.01;
        fhi = pressureRegion3(hi, T, dPdrho) - P;
    }
    if (flo > 0.0 || fhi < 0.0)
    {
        throw std::runtime_error(rangeError("could not bracket the region-3 density", T, P));
    }

    double rho = 0.5*(lo + hi);
    for (int it=0; it < 200; ++it)
    {
        const double f = pressureRegion3(rho, T, dPdrho) - P;
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

    throw std::runtime_error(rangeError("region-3 density iteration did not converge", T, P));
}

/* Isobaric thermal expansion coefficient in region 3 (no member state is touched). */
double water::region3Alpha(double T, double P)
{
    const double rho = region3Density(T, P);
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

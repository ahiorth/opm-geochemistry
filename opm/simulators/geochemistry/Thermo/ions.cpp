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
#include <opm/simulators/geochemistry/Thermo/ions.h>

/* Based on data from:
 *
 *		Calculation of the thermodynamic properties of aqueous species at high pressures
 *		and temperatures.Effective electrostatic radii, dissociation constants and standard
 *		partial molal properties to 1000 �C and 5 kbar
 *
 *	Everett L. Shock,   Eric H. Oelkers,   James W. Johnson,   Dimitri A. Sverjensky and   Harold C. Helgeson
 *	J. Chem. Soc., Faraday Trans., 1992,88, 803-826
 *
 */

ions::ions(water* W)
: W_(W)
, eta_(1.66027e5 * 4.184)  //  J/mol
, born_f_(0.0)
, born_f_P_(0.0)
, born_f_T_(0.0)
, born_f_TT_(0.0)
, born_g_(0.0)
, born_g_P_(0.0)
, born_g_T_(0.0)
, born_g_TT_(0.0)
, born_gp2_(0.0)
, u_(0.0)
, u_T_(0.0)
, u_TT_(0.0)
{
}

/* T in Kelvin and P in Pascal.
 * The f correction is restricted to Psat < P < 1000 bar and 155 < T < 355 C.
 * The density-dependent part of g can remain nonzero outside that window. */
void ions::born_df(double T, double P)
{
    static constexpr double ag[] = { 0., -0.2037662e1, 0.5747e-2, -0.6557892e-5 };
    static constexpr double bg[] = { 0., 0.6107361e1, -0.1074377e-1, 0.1268348e-4 };
    static constexpr double af[] = { 0., 0.3666666e2, -0.1504956e-9, 0.5017997e-13 };
    
    double agg, bgg;
    
    born_f_ = born_f_P_ = born_f_T_ = born_f_TT_ = 0.;

    const double temperature_celsius = T - 273.15;
    const bool within_f_domain = temperature_celsius > 155.0
                                 && temperature_celsius < 355.0
                                 && P > W_->Psat_
                                 && P < 1000.0e5;
    double Tr = (temperature_celsius - 155.0) / 300.0;
    if (within_f_domain)
    {
        const double Pr = 1000.0 - P * 1.e-5;
        const double Tr_inv = 1. / Tr;
        const double Tr_p1 = pow(Tr, 4.8);
        const double Tr_p2 = pow(Tr, 16);
        const double P_p1 = Pr * Pr * Pr;
        const double P_p2 = P_p1 * Pr;
        born_f_ = Tr_p1 + af[1] * Tr_p2;
        born_f_P_ = born_f_;
        born_f_TT_ = born_f_T_ = af[2] * P_p1 + af[3] * P_p2;
        born_f_ *= born_f_T_;
        born_f_P_ *= -3. * af[2] * Pr * Pr - 4. * af[3] * P_p1;
        born_f_T_ *= 0.016 * Tr_p1 + 16. / 300. * af[1] * Tr_p2;
        born_f_T_ *= Tr_inv;
        born_f_TT_ *= 0.0608 / 300. * Tr_p1 + af[1] / 375. * Tr_p2;
        born_f_TT_ *= Tr_inv * Tr_inv;
    }
    
    
    Tr = temperature_celsius;
    agg = ag[1] + ag[2] * Tr + ag[3] * Tr * Tr;
    bgg = bg[1] + bg[2] * Tr + bg[3] * Tr * Tr;

    double agg_T = ag[2] + 2. * ag[3] * Tr;
    double agg_TT = 2. * ag[3];
    double bgg_T = bg[2] + 2. * bg[3] * Tr;
    double bgg_TT = 2. * bg[3];
    
    const double rho_ref = 1.e-3*W_->denst_;
    const double rho_1 = 1. - rho_ref;
    if (rho_1 > 0.)
    {
        const double rho_1_inv = 1. / rho_1;
        const double rho_p = pow(rho_1, bgg);
        const double log_rho_1 = log(rho_1);
        const double rho_ratio = rho_ref * rho_1_inv;

        // These identities use only rho_T = -alpha*rho and rho_P = beta*rho,
        // so they apply to every smooth EOS region that supplies alpha, beta,
        // and alpha_t = (d alpha/dT)_P.
        const double dlog_rho_p_dT =
            bgg_T * log_rho_1 + bgg * W_->alpha_ * rho_ratio;
        const double d2log_rho_p_dT2 =
            bgg_TT * log_rho_1
            + 2. * bgg_T * W_->alpha_ * rho_ratio
            + bgg * W_->alpha_t_ * rho_ratio
            - bgg * W_->alpha_ * W_->alpha_ * rho_ref
                * rho_1_inv * rho_1_inv;
        const double dlog_rho_p_dPbar =
            -bgg * W_->beta_ * 1.e5 * rho_ratio;

        const double rho_p_T = rho_p * dlog_rho_p_dT;
        const double rho_p_TT =
            rho_p * (dlog_rho_p_dT * dlog_rho_p_dT + d2log_rho_p_dT2);
        u_ = rho_p;
        u_T_ = rho_p_T;
        u_TT_ = rho_p_TT;

        born_g_ = agg * rho_p - born_f_;
        born_g_P_ = agg * rho_p * dlog_rho_p_dPbar - born_f_P_;
        born_g_T_ = agg_T * rho_p + agg * rho_p_T - born_f_T_;
        born_g_TT_ = agg_TT * rho_p + 2. * agg_T * rho_p_T + agg * rho_p_TT - born_f_TT_;
    }
    else {
        u_ = 0.;
        u_T_ = 0.;
        u_TT_ = 0.;
        born_g_ = 0.;
        born_g_P_ = 0.;
        born_g_T_ = 0.;
        born_g_TT_ = 0.;
    }
    
    
    
    
    born_gp2_ = born_g_ + 3.082;
    born_gp2_ *= born_gp2_;
}

/* Legacy temperature-only f/g evaluator. The pressure-aware path is born_df. */
void ions::born_f(double T)
{
    double ag[] = { 0., -0.2037662e1, 0.5747e-2, -0.6557892e-5 };
    double bg[] = { 0., 0.6107361e1, -0.1074377e-1, 0.1268348e-4 };
    double af[] = { 0., 0.3666666e2, -0.1504956e-9, 0.5017997e-13 };
    
    //double agg, bgg;
    //double Tr, Tr_p1, Tr_p2;
    const double temperature_celsius = T - 273.15;
    double Tr = (temperature_celsius - 155.0) / 300.0;
    if (temperature_celsius > 155.0 && temperature_celsius < 355.0)
    {
        const double Tr_p1 = pow(Tr, 4.8);
        const double Tr_p2 = pow(Tr, 16);
        born_f_ = Tr_p1 + af[1] * Tr_p2;
    }
    else{
        born_f_ = 0.;
    }
    
    Tr = temperature_celsius;
    const double agg = ag[1] + ag[2] * Tr + ag[3] * Tr * Tr;
    const double bgg = bg[1] + bg[2] * Tr + bg[3] * Tr * Tr;
    
    double rho_ref = W_->denst_ * 1.e-3;
    double rho_1 = 1. - rho_ref;
    if (rho_1 > 0.) {
        const double rho_p = pow(rho_1, bgg);
        born_g_ = agg * rho_p - born_f_;
    }
    else{
        born_g_ = 0.;
    }
    
    born_gp2_ = born_g_ + 3.082;
    born_gp2_ *= born_gp2_;
}

/*
 Calculation of the thermodynamic properties of aqueous species at high pressures
 and temperatures.Effective electrostatic radii, dissociation constants and standard
 partial molal properties to 1000 �C and 5 kbar
 
 Everett L. Shock,   Eric H. Oelkers,   James W. Johnson,   Dimitri A. Sverjensky and   Harold C. Helgeson
 J. Chem. Soc., Faraday Trans., 1992,88, 803-826
 
 double rho_ref = W_->denst_*1.e-3;
 if (rho_ref < 1.) born_df(T, P);
 else born_g_ = born_g_P_ = born_g_T_ = born_g_TT_ = 0.;
 */
/* assumes that born_f or born_df has been called before */
void ions::born(double Z, double& re_ref, double& w, double& w_T, double& w_TT, double& w_P)
{
    const double Zabs = fabs(Z);
    const double re = re_ref + Zabs * born_g_; /* effective ion radius @ P, T*/
    double ff = Zabs * Zabs / re;
    w = eta_ * (ff - Z / (3.082 + born_g_));
    ff *= Zabs / re;
    const double ff4 = ff * Zabs / re;
    ff -= Z / born_gp2_;
    w_P = -eta_ * ff * born_g_P_;
    w_T = -eta_ * ff * born_g_T_;
    w_TT = -eta_ * ff * born_g_TT_;
    w_TT += 2. * eta_ * (ff4 - Z / born_gp2_ / (3.082 + born_g_)) * born_g_T_ * born_g_T_;
}

void ions::born(double Z, double& re_ref, double& w)
{
    const double Zabs = fabs(Z);
    if (Zabs < 1e-8) {
        w = 0.;
    }
    else
    {
        const double re = re_ref + Zabs * born_g_; /* effective ion radius @ P, T*/
        w = eta_ * (Zabs * Zabs / re - Z / (3.082 + born_g_));
    }
}

void ions::compare_analytical_and_numerical_derivatives(double T, double P, double dT, double dP)
{
    W_->gibbsIAPWS(T, P);
    born_df(T, P);
    printf("Analytical:\n");
    printf("f\tdf/dP\tdf/dt\td2f/dT2\n");
    printf("%g\t%g\t%g\t%g\n", born_f_, born_f_P_, born_f_T_, born_f_TT_);
    printf("%g\t%g\t%g\t%g\n", born_g_, born_g_P_, born_g_T_, born_g_TT_);
    printf("%g\t%g\t%g\t%g\n", u_, 0., u_T_, u_TT_);
    const double f = born_f_;
    const double g = born_g_;
    const double u = u_;
    
    W_->gibbsIAPWS(T + dT, P);
    born_df(T + dT, P);
    const double df_T = born_f_;
    const double dg_T = born_g_;
    const double du_T = u_;
    
    W_->gibbsIAPWS(T - dT, P);
    born_df(T - dT, P);
    const double df_T2 = born_f_;
    const double dg_T2 = born_g_;
    const double du_T2 = u_;
    
    W_->gibbsIAPWS(T, P - dP);
    born_df(T, P - dP);
    const double df_P2 = born_f_;
    const double dg_P2 = born_g_;
    
    W_->gibbsIAPWS(T, P + dP);
    born_df(T, P + dP);
    
    const double df_P = born_f_;
    const double dg_P = born_g_;
    printf("Numerical:\n");
    printf("f\tdf/dP\tdf/dt\td2f/dT2\n");
    printf("%g\t%g\t%g\t%g\n", f, (df_P - df_P2) / (2. * dP * 1e-5), (df_T - df_T2) / (2. * dT), (df_T + df_T2 - 2. * f) / dT / dT);
    printf("%g\t%g\t%g\t%g\n", g, (dg_P - dg_P2) / (2. * dP * 1e-5), (dg_T - dg_T2) / (2. * dT), (dg_T + dg_T2 - 2. * g) / dT / dT);
    printf("%g\t%g\t%g\t%g\n", u, 0., (du_T - du_T2) / (2. * dT), (du_T + du_T2 - 2. * u) / dT / dT);
    
}

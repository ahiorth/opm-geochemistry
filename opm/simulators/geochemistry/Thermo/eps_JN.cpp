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
#include <opm/simulators/geochemistry/Thermo/eps_JN.h>

#include <stdexcept>
#include <string>

/* Water Psat_ is NaN above the critical temperature, so the comparison is
 * false there and supercritical states pass. */
void eps_JN::requireLiquidOrSupercritical(double T, double P) const
{
    if (P < W_->Psat_)
    {
        throw std::domain_error(
            "eps_JN (Johnson-Norton): water is vapour at T=" + std::to_string(T)
            + " K, P=" + std::to_string(P) + " Pa (below the saturation pressure "
            + std::to_string(W_->Psat_) + " Pa). The dielectric constant and the"
            " HKF aqueous-species properties require liquid or supercritical water.");
    }
}

eps_JN::eps_JN(water* W)
: permittivity_(0.0)
, permittivity_P_(0.0)
, permittivity_T_(0.0)
, permittivity_TT_(0.0)
, bornX_(0.0)
, bornY_(0.0)
, bornZ_(0.0)
, bornQ_(0.0)
, W_(W)
, Tref_(298.15)
, Tref_inv_(1.0/298.15)
, a_()
, k_(5, 0.0)
, k_t_(5, 0.0)
, k_tt_(5, 0.0)
{
    a_.push_back(0.);
    a_.push_back(0.1470333593e2);
    a_.push_back(0.2128462733e3);
    a_.push_back(-0.1154445173e3);
    a_.push_back(0.1955210915e2);
    a_.push_back(-0.8330347980e2);
    a_.push_back(0.3213240048e2);
    a_.push_back(-0.6694098645e1);
    a_.push_back(-0.3786202045e2);
    a_.push_back(0.6887359646e2);
    a_.push_back(-0.2729401652e2);
}

void eps_JN::ke(double T)
{
    const double Ts = T / Tref_;
    const double Ts_inv = 1. / Ts;
    
    k_[0] = 1.0;
    k_[1] = a_[1] * Ts_inv;
    k_[2] = a_[2] * Ts_inv + a_[3] + a_[4] * Ts;
    k_[3] = a_[5] * Ts_inv + a_[6] * Ts + a_[7] * Ts * Ts;
    k_[4] = a_[8] * Ts_inv * Ts_inv + a_[9] * Ts_inv + a_[10];
}


void eps_JN::ke_t(double T)
{
    const double Ts = T / Tref_;
    const double Ts_inv = 1. / Ts;
    
    k_t_[0] = 0.0;
    k_t_[1] = -a_[1] * Tref_inv_ * Ts_inv * Ts_inv;
    k_t_[2] = (-a_[2] * Ts_inv * Ts_inv + a_[4]) * Tref_inv_;
    k_t_[3] = (-a_[5] * Ts_inv * Ts_inv + a_[6] + 2. * a_[7] * Ts) * Tref_inv_;
    k_t_[4] = (-2. * a_[8] * Ts_inv * Ts_inv * Ts_inv - a_[9] * Ts_inv * Ts_inv) * Tref_inv_;
}

void eps_JN::ke_tt(double T)
{
    const double Ts = T / Tref_;
    const double Ts_inv = 1. / Ts;
    
    k_tt_[0] = 0.0;
    k_tt_[1] = 2. * a_[1] * Tref_inv_ * Tref_inv_ * Ts_inv * Ts_inv * Ts_inv;
    k_tt_[2] = 2. * a_[2] * Tref_inv_ * Tref_inv_ * Ts_inv * Ts_inv * Ts_inv;
    k_tt_[3] = (2. * a_[5] * Ts_inv * Ts_inv * Ts_inv + 2. * a_[7]) * Tref_inv_ * Tref_inv_;
    k_tt_[4] = (6. * a_[8] * Ts_inv * Ts_inv * Ts_inv * Ts_inv + 2. * a_[9] * Ts_inv * Ts_inv * Ts_inv) * Tref_inv_ * Tref_inv_;
}

/* permittivity of water*/
/* calls water properties - if unknown*/
void eps_JN::permittivity(double T, double P)
{
    ke(T);
    W_->gibbsIAPWS(T, P);
    requireLiquidOrSupercritical(T, P);

    const double rho_1 = 1.0e-3 / W_->v_;
    const double rho_2 = rho_1 * rho_1;
    const double rho_3 = rho_1 * rho_2;
    const double rho_4 = rho_1 * rho_3;
    const double beta = W_->beta_;
    
    permittivity_ = k_[0] + k_[1]*rho_1 + k_[2]*rho_2 + k_[3]*rho_3 + k_[4]*rho_4;
    permittivity_P_ = beta*(k_[1]*rho_1 + 2.0*k_[2]*rho_2 + 3.0*k_[3]*rho_3 + 4.0*k_[4]*rho_4);
    bornZ_ = -1.0 / permittivity_;
    bornQ_ = bornZ_*bornZ_*permittivity_P_;
}

/* permittivity of water and its derivatives */
void eps_JN::permittivity_TP(double T, double P)
{
    
    ke(T);
    ke_t(T);
    ke_tt(T);
    W_->gibbsIAPWS(T, P);
    requireLiquidOrSupercritical(T, P);

    const double alpha = W_->alpha_;
    const double beta = W_->beta_;
    const double alpha_t = W_->alpha_t_;
    
    const double rho_1 = 1.0e-3 / W_->v_;
    const double rho_2 = rho_1 * rho_1;
    const double rho_3 = rho_1 * rho_2;
    const double rho_4 = rho_1 * rho_3;
    
    permittivity_ = k_[0] + k_[1] * rho_1 + k_[2] * rho_2 + k_[3] * rho_3 + k_[4] * rho_4;
    permittivity_P_ = beta*(k_[1] * rho_1 + 2. * k_[2] * rho_2 + 3. * k_[3] * rho_3 + 4. * k_[4] * rho_4);
    
    permittivity_T_ = k_t_[0] + (k_t_[1] - alpha * k_[1]) * rho_1
    + (k_t_[2] - 2. * alpha * k_[2]) * rho_2
    + (k_t_[3] - 3. * alpha * k_[3]) * rho_3
    + (k_t_[4] - 4. * alpha * k_[4]) * rho_4;
    
    permittivity_TT_ = k_tt_[0]
    + (k_tt_[1] - 2. * alpha * k_t_[1] - alpha_t * k_[1] + alpha * alpha * k_[1]) * rho_1
    + (k_tt_[2] - 4. * alpha * k_t_[2] - 2. * alpha_t * k_[2] + 4. * alpha * alpha * k_[2]) * rho_2
    + (k_tt_[3] - 6. * alpha * k_t_[3] - 3. * alpha_t * k_[3] + 9. * alpha * alpha * k_[3]) * rho_3
    + (k_tt_[4] - 8. * alpha * k_t_[4] - 4. * alpha_t * k_[4] + 16. * alpha * alpha * k_[4]) * rho_4;
    
    bornZ_ = -1.0 / permittivity_;
    
    const double z2 = bornZ_ * bornZ_;
    bornX_ = permittivity_TT_ + 2.0 * bornZ_ * permittivity_T_ * permittivity_T_;
    bornX_ *= z2;
    
    bornY_ = z2 * permittivity_T_;
    bornQ_ = z2 * permittivity_P_;
}

void eps_JN::print_permittivity(double T, double P)
{
    permittivity_TP(T, P);
    printf("T=%g\tP=%g\n", T, P);
    printf("permittivity\tperm_P\tperm_T\tperm_TT\n");
    printf("%g\t%g\t%g\t%g\n", permittivity_, permittivity_P_, permittivity_T_, permittivity_TT_);
}

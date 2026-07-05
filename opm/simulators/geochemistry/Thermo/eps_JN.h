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
#ifndef EPS_JN_H
#define EPS_JN_H

#include <iostream>
#include <stdlib.h>
#include <vector>

#include <opm/simulators/geochemistry/Thermo/water.h>

//Johnson, J.W.and Norton, D. (1991) Critical phenomena in hydrothermal systems : 
//state, thermodynamic, electrostatic, and transport properties of H2O in the critical region.
// Am.J.Sci. 291, 541--648.
// http ://www.ajsonline.org/cgi/content/abstract/291/6/541

class eps_JN {
    
public:
    
    eps_JN(water* W);
    
    void permittivity(double T, double P);
    void permittivity_TP(double T, double P);
    void print_permittivity(double T=298.15, double P=1.0e5);
    
    double permittivity_;  // dimensionless
    double permittivity_P_;  // 1/Pa
    double permittivity_T_;  // 1/K
    double permittivity_TT_;  // 1/K^2
    double bornX_;  // 1/K^2
    double bornY_;  // 1/K
    double bornZ_;  // dimensionless
    double bornQ_;  // 1/Pa
    
private:
    
    water* W_;
    
    double Tref_;
    double Tref_inv_;
    
    std::vector<double> a_;
    std::vector<double> k_;
    std::vector<double> k_t_;
    std::vector<double> k_tt_;
    
    void ke(double T);
    void ke_t(double T);
    void ke_tt(double T);

    /* The Johnson-Norton correlation (and the HKF machinery built on it)
     * applies to liquid and supercritical water only; throws
     * std::domain_error for sub-saturation (vapour) states, which the water
     * EOS itself computes normally. Assumes W_->gibbsIAPWS(T, P) has been
     * called. */
    void requireLiquidOrSupercritical(double T, double P) const;
    
};

#endif

#if __has_include(<catch2/catch.hpp>)
#include <catch2/catch.hpp>
#else
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#endif

#include <array>
#include <cmath>
#include <map>
#include <stdexcept>
#include <string>
#include <fstream>
#include <sstream>

#include <opm/simulators/geochemistry/StandaloneSolvers.hpp>
#include <opm/simulators/geochemistry/Thermo/hkf.h>
#include <opm/simulators/geochemistry/Thermo/eps_JN.h>
#include <opm/simulators/geochemistry/Thermo/ions.h>
#include <opm/simulators/geochemistry/Thermo/ThermoTable.h>
#include <opm/simulators/geochemistry/Thermo/thermodata.h>
#include <opm/simulators/geochemistry/Thermo/water.h>

static constexpr double abs_tolerance_ = 1.0e-7;
static constexpr double abs_tolerance2_ = 1.0e-4;  // for most water props.

/*
* Stores reference cases 1, 2, and 3 from Table 5 of the IAPWS-97
* industrial thermodynamic formulation.
*
* See also:
*
*   - https://github.com/jjgomera/iapws.git
*   - http://twt.mpei.ac.ru/mcs/worksheets/iapws/IAPWS-IF97-Region1.xmcd
*   - http://www.casprod.eu/~leon/research/iapws-if97/listing.html
*/
struct WaterPropertiesIAPWS_97
{
    WaterPropertiesIAPWS_97(int reference_case)
    {
        switch (reference_case) {
            case 1:
            {
                pressure_ = 3.0e6;
                temperature_ = 300.0;
                v_ = 0.100215168e-2;
                h_ = 0.115331273e3;
                u_ = 0.112324818e3;
                s_ = 0.392294792;
                cp_ = 0.417301218e1;
                w_ = 0.150773921e4;
                break;
            }
            case 2:
            {
                pressure_ = 80.0e6;
                temperature_ = 300.0;
                v_ = 0.971180894e-3;
                h_ = 0.184142828e3;
                u_ = 0.106448356e3;
                s_ = 0.368563852;
                cp_ = 0.401008987e1;
                w_ = 0.163469054e4;
                break;
            }
            case 3:
            {
                pressure_ = 3.0e6;
                temperature_ = 500.0;
                v_ = 0.120241800e-2;
                h_ = 0.975542239e3;
                u_ = 0.971934985e3;
                s_ = 0.258041912e1;
                cp_ = 0.465580682e1;
                w_ = 0.124071337e4;
                break;
            }
        }
    }

    double pressure_;
    double temperature_;

    double v_;  // specific volume (m3/kg)
    double h_;  // specific enthalpy (kJ/kg)
    double u_;  // specific internal energy (kJ/kg)
    double s_;  // specific entropy (kJ/kg/K)
    double cp_;  // specific isobaric heat capacity (kJ/kg/K)
    double w_;  // speed of sound (m/s)

    //double cv_;  // specific isochoric heat capacity (kJ/kg/K)
};


TEST_CASE("Test calculation of thermodynamic properties for water")
{

    water water_props;

    static const std::array<WaterPropertiesIAPWS_97, 3> reference_cases =
    {
        WaterPropertiesIAPWS_97(1),
        WaterPropertiesIAPWS_97(2),
        WaterPropertiesIAPWS_97(3),
    };

    for(const auto& ref_props: reference_cases)
    {
        const double T = ref_props.temperature_;
        const double P = ref_props.pressure_;

        water_props.gibbsIAPWS(T, P);

        CHECK_THAT(water_props.v_, Catch::Matchers::WithinAbs(ref_props.v_, abs_tolerance_));

        // Note: Reference results are in kJ/kg or kJ/kg/K, we use Joules...
        CHECK_THAT(1.0e-3*water_props.h_, Catch::Matchers::WithinAbs(ref_props.h_, abs_tolerance2_));
        CHECK_THAT(1.0e-3*water_props.u_, Catch::Matchers::WithinAbs(ref_props.u_, abs_tolerance2_));
        CHECK_THAT(1.0e-3*water_props.s_, Catch::Matchers::WithinAbs(ref_props.s_, abs_tolerance2_));
        CHECK_THAT(1.0e-3*water_props.cp_, Catch::Matchers::WithinAbs(ref_props.cp_, abs_tolerance2_));
        CHECK_THAT(water_props.w_, Catch::Matchers::WithinAbs(ref_props.w_, abs_tolerance2_));

        water_props.printProperties();
    }

}

TEST_CASE("Test calculation of saturation pressure for water")
{
    water water_props;

    // Table 35 of IAPWS97 paper
    static constexpr std::array<std::pair<double, double>, 3> table35_tests =
    {
        std::make_pair<double, double>(300.0, 0.353658941e-2),
        std::make_pair<double, double>(500.0, 0.263889776e1),
        std::make_pair<double, double>(600.0, 0.123443146e2)
    };

    for(const auto& [T, Psat] : table35_tests)
    {
        const double satPress = 1.0e-6*water_props.PsatIAPWS(T);
        CHECK_THAT(satPress, Catch::Matchers::WithinAbs(Psat, abs_tolerance_));
    }
}

/* Reference values for regions 2 and 3 of IAPWS-97, Tables 15 and 33.
 * For region 3 the paper tabulates (T, rho) inputs; the pressures below are the
 * corresponding table values, so the returned density must reproduce rho. */
struct WaterPropertiesHighTP
{
    double pressure_;     // Pa
    double temperature_;  // K
    int region_;
    double v_;   // specific volume (m3/kg)
    double h_;   // specific enthalpy (kJ/kg)
    double u_;   // specific internal energy (kJ/kg)
    double s_;   // specific entropy (kJ/kg/K)
    double cp_;  // specific isobaric heat capacity (kJ/kg/K)
    double w_;   // speed of sound (m/s)
};

TEST_CASE("Test calculation of water properties in IAPWS-97 region 2 (steam)")
{
    water water_props;

    // Table 15 of IAPWS-97 paper
    static constexpr std::array<WaterPropertiesHighTP, 3> reference_cases =
    {{
        { 0.0035e6, 300.0, 2, 0.394913866e2, 0.254991145e4, 0.241169160e4, 0.852238967e1, 0.191300162e1, 0.427920172e3 },
        { 0.0035e6, 700.0, 2, 0.923015898e2, 0.333568375e4, 0.301262819e4, 0.101749996e2, 0.208141274e1, 0.644289068e3 },
        { 30.0e6,   700.0, 2, 0.542946619e-2, 0.263149474e4, 0.246861076e4, 0.517540298e1, 0.103505092e2, 0.480386523e3 },
    }};

    for(const auto& ref_props: reference_cases)
    {
        water_props.gibbsIAPWS(ref_props.temperature_, ref_props.pressure_);

        CHECK(water_props.region_ == ref_props.region_);
        CHECK_THAT(water_props.v_, Catch::Matchers::WithinRel(ref_props.v_, 1.0e-8));
        CHECK_THAT(1.0e-3*water_props.h_, Catch::Matchers::WithinAbs(ref_props.h_, abs_tolerance2_));
        CHECK_THAT(1.0e-3*water_props.u_, Catch::Matchers::WithinAbs(ref_props.u_, abs_tolerance2_));
        CHECK_THAT(1.0e-3*water_props.s_, Catch::Matchers::WithinAbs(ref_props.s_, abs_tolerance2_));
        CHECK_THAT(1.0e-3*water_props.cp_, Catch::Matchers::WithinAbs(ref_props.cp_, abs_tolerance2_));
        CHECK_THAT(water_props.w_, Catch::Matchers::WithinAbs(ref_props.w_, abs_tolerance2_));
    }
}

TEST_CASE("Test calculation of water properties in IAPWS-97 region 3 (supercritical)")
{
    water water_props;

    // Table 33 of IAPWS-97 paper: inputs there are (T, rho) = (650, 500), (650, 200)
    // and (750, 500); the pressures below are the tabulated P(T, rho), so the
    // density iteration has to recover rho.
    static constexpr std::array<WaterPropertiesHighTP, 3> reference_cases =
    {{
        { 0.255837018e8, 650.0, 3, 1.0/500.0, 0.186343019e4, 0.181226279e4, 0.405427273e1, 0.138935717e2, 0.502005554e3 },
        { 0.222930643e8, 650.0, 3, 1.0/200.0, 0.237512401e4, 0.226365868e4, 0.485438792e1, 0.446579342e2, 0.383444594e3 },
        { 0.783095639e8, 750.0, 3, 1.0/500.0, 0.225868845e4, 0.210206932e4, 0.446971906e1, 0.634165359e1, 0.760696041e3 },
    }};

    for(const auto& ref_props: reference_cases)
    {
        water_props.gibbsIAPWS(ref_props.temperature_, ref_props.pressure_);

        CHECK(water_props.region_ == ref_props.region_);
        CHECK_THAT(water_props.v_, Catch::Matchers::WithinRel(ref_props.v_, 1.0e-6));
        CHECK_THAT(1.0e-3*water_props.h_, Catch::Matchers::WithinAbs(ref_props.h_, abs_tolerance2_));
        CHECK_THAT(1.0e-3*water_props.u_, Catch::Matchers::WithinAbs(ref_props.u_, abs_tolerance2_));
        CHECK_THAT(1.0e-3*water_props.s_, Catch::Matchers::WithinAbs(ref_props.s_, abs_tolerance2_));
        CHECK_THAT(1.0e-3*water_props.cp_, Catch::Matchers::WithinAbs(ref_props.cp_, 1.0e-3));
        CHECK_THAT(water_props.w_, Catch::Matchers::WithinAbs(ref_props.w_, 1.0e-3));
    }
}

TEST_CASE("Test IAPWS-97 region boundaries and dispatch")
{
    water water_props;

    // B23 boundary verification pair from the IAPWS-97 paper (below Eq. 6):
    // T = 623.15 K <-> P = 16.5291643 MPa
    CHECK_THAT(1.0e-6*water::PB23IAPWS(623.15), Catch::Matchers::WithinAbs(16.5291643, 1.0e-6));

    // Region 1 (liquid) just below saturation temperature, steam just below Psat
    water_props.gibbsIAPWS(600.0, 15.0e6);
    CHECK(water_props.region_ == 1);
    water_props.gibbsIAPWS(600.0, 10.0e6);  // Psat(600) = 12.34 MPa
    CHECK(water_props.region_ == 2);
    water_props.gibbsIAPWS(900.0, 50.0e6);  // high-T steam, above B23 range
    CHECK(water_props.region_ == 2);

    // Property continuity across the region 1/3 boundary (T = 623.15 K isotherm).
    // IAPWS-97 guarantees consistency to within ~0.05% in v at the boundaries;
    // test with a slightly looser tolerance.
    water water_r1;
    water water_r3;
    water_r1.gibbsIAPWS(623.15, 50.0e6);
    CHECK(water_r1.region_ == 1);
    water_r3.gibbsIAPWS(623.16, 50.0e6);
    CHECK(water_r3.region_ == 3);
    CHECK_THAT(water_r3.v_, Catch::Matchers::WithinRel(water_r1.v_, 5.0e-3));
    CHECK_THAT(water_r3.h_, Catch::Matchers::WithinRel(water_r1.h_, 5.0e-3));
    CHECK_THAT(water_r3.s_, Catch::Matchers::WithinRel(water_r1.s_, 5.0e-3));
    CHECK_THAT(water_r3.alpha_, Catch::Matchers::WithinRel(water_r1.alpha_, 2.0e-2));
    CHECK_THAT(water_r3.beta_, Catch::Matchers::WithinRel(water_r1.beta_, 2.0e-2));

    // Property continuity across the region 2/3 (B23) boundary at 700 K
    const double PB23_700 = water::PB23IAPWS(700.0);
    water water_r2;
    water_r2.gibbsIAPWS(700.0, 0.999*PB23_700);
    CHECK(water_r2.region_ == 2);
    water_r3.gibbsIAPWS(700.0, 1.001*PB23_700);
    CHECK(water_r3.region_ == 3);
    CHECK_THAT(water_r3.v_, Catch::Matchers::WithinRel(water_r2.v_, 2.0e-2));
    CHECK_THAT(water_r3.h_, Catch::Matchers::WithinRel(water_r2.h_, 2.0e-2));
    CHECK_THAT(water_r3.s_, Catch::Matchers::WithinRel(water_r2.s_, 2.0e-2));
}

TEST_CASE("Test water properties above 100 MPa with IAPWS-95")
{
    water water_props;

    // Single-phase reference states from Table 7 of the IAPWS-95 release
    // (Wagner & Pruss 2002, Table 6.6): inputs there are (T, rho); the
    // pressures below are the tabulated P(T, rho), so the density iteration
    // has to recover rho. Only the 700 MPa states exceed the IAPWS-97 range
    // and are reachable through the (T, P) interface.
    struct Iapws95Reference
    {
        double pressure_;     // Pa
        double temperature_;  // K
        double rho_;  // density (kg/m3)
        double cv_;   // specific isochoric heat capacity (kJ/kg/K)
        double w_;    // speed of sound (m/s)
        double s_;    // specific entropy (kJ/kg/K)
    };

    static constexpr std::array<Iapws95Reference, 3> reference_cases =
    {{
        { 0.700004704e9, 300.0, 0.1188202e4, 0.346135580e1, 0.244357992e4, 0.132609616 },
        { 0.700000405e9, 500.0, 0.1084564e4, 0.307437693e1, 0.241200877e4, 0.203237509e1 },
        { 0.700000006e9, 900.0, 0.8707690e3, 0.266422350e1, 0.201933608e4, 0.417223802e1 },
    }};

    for(const auto& ref_props: reference_cases)
    {
        water_props.gibbsIAPWS(ref_props.temperature_, ref_props.pressure_);

        CHECK(water_props.region_ == 95);
        CHECK_THAT(water_props.denst_, Catch::Matchers::WithinRel(ref_props.rho_, 1.0e-6));
        CHECK_THAT(1.0e-3*water_props.cv_, Catch::Matchers::WithinAbs(ref_props.cv_, abs_tolerance2_));
        CHECK_THAT(water_props.w_, Catch::Matchers::WithinAbs(ref_props.w_, 1.0e-3));
        CHECK_THAT(1.0e-3*water_props.s_, Catch::Matchers::WithinAbs(ref_props.s_, 1.0e-6));
    }

    // Continuity across the 100 MPa seam between IAPWS-97 and IAPWS-95
    // (IAPWS-97 was fitted to IAPWS-95; consistency is a few 0.01%).
    water water_97;
    water water_95;

    water_97.gibbsIAPWS(500.0, 99.99e6);  // region 1 side
    CHECK(water_97.region_ == 1);
    water_95.gibbsIAPWS(500.0, 100.01e6);
    CHECK(water_95.region_ == 95);
    CHECK_THAT(water_95.v_, Catch::Matchers::WithinRel(water_97.v_, 2.0e-3));
    CHECK_THAT(water_95.h_, Catch::Matchers::WithinRel(water_97.h_, 2.0e-3));
    CHECK_THAT(water_95.s_, Catch::Matchers::WithinRel(water_97.s_, 2.0e-3));
    CHECK_THAT(water_95.alpha_, Catch::Matchers::WithinRel(water_97.alpha_, 1.0e-2));
    CHECK_THAT(water_95.beta_, Catch::Matchers::WithinRel(water_97.beta_, 1.0e-2));

    water_97.gibbsIAPWS(650.0, 99.99e6);  // region 3 side (exercises the near-critical terms)
    CHECK(water_97.region_ == 3);
    water_95.gibbsIAPWS(650.0, 100.01e6);
    CHECK(water_95.region_ == 95);
    CHECK_THAT(water_95.v_, Catch::Matchers::WithinRel(water_97.v_, 2.0e-3));
    CHECK_THAT(water_95.h_, Catch::Matchers::WithinRel(water_97.h_, 2.0e-3));
    CHECK_THAT(water_95.s_, Catch::Matchers::WithinRel(water_97.s_, 2.0e-3));

    water_97.gibbsIAPWS(1000.0, 99.99e6);  // region 2 side
    CHECK(water_97.region_ == 2);
    water_95.gibbsIAPWS(1000.0, 100.01e6);
    CHECK(water_95.region_ == 95);
    CHECK_THAT(water_95.v_, Catch::Matchers::WithinRel(water_97.v_, 2.0e-3));
    CHECK_THAT(water_95.h_, Catch::Matchers::WithinRel(water_97.h_, 2.0e-3));
    CHECK_THAT(water_95.s_, Catch::Matchers::WithinRel(water_97.s_, 2.0e-3));
}

TEST_CASE("HKF machinery rejects sub-saturation (steam) water states")
{
    hkf hkf_props;

    // 150 C at 1 bar is steam (Psat = 4.76 bar): the water EOS computes it,
    // but the dielectric/HKF chain requires liquid or supercritical water.
    CHECK_THROWS_AS(hkf_props.WaterProp(423.15, 1.0e5), std::domain_error);

    // Same temperature above the saturation pressure works.
    hkf_props.WaterProp(423.15, 1.0e6);
    CHECK(hkf_props.rhow_ > 900.0);

    // The pure-water EOS itself still handles steam.
    water water_props;
    water_props.gibbsIAPWS(423.15, 1.0e5);
    CHECK(water_props.region_ == 2);
}

TEST_CASE("Region 3 alpha_t is consistent and stays on one branch near saturation")
{
    static constexpr double T = 630.0;
    const double Psat = water::PsatIAPWS(T);

    // Validate alpha_t_ against a coarse finite difference of alpha_ at a
    // pressure safely above the saturation line (all probes liquid).
    water water_props;
    water_props.gibbsIAPWS(T, Psat + 5.0e5);
    CHECK(water_props.region_ == 3);
    const double alpha_t_ref = water_props.alpha_t_;

    water water_plus;
    water water_minus;
    water_plus.gibbsIAPWS(T + 0.5, Psat + 5.0e5);
    water_minus.gibbsIAPWS(T - 0.5, Psat + 5.0e5);
    const double alpha_t_fd = (water_plus.alpha_ - water_minus.alpha_) / 1.0;
    CHECK_THAT(alpha_t_ref, Catch::Matchers::WithinRel(alpha_t_fd, 0.15));

    // Within one finite-difference step of the saturation line
    // (dPsat/dT ~ 0.35 MPa/K, probe dT = 0.01 K): the internal probes must
    // stay on the liquid branch instead of flipping to the vapour root,
    // which would produce a wildly wrong derivative.
    water_props.gibbsIAPWS(T, Psat + 2.0e3);
    CHECK(water_props.denst_ > 500.0);  // liquid branch
    CHECK(std::isfinite(water_props.alpha_t_));
    CHECK(water_props.alpha_t_*alpha_t_ref > 0.0);  // same sign
    CHECK(std::fabs(water_props.alpha_t_) < 20.0*std::fabs(alpha_t_ref));
}

TEST_CASE("Test IAPWS-97 out-of-range conditions raise exceptions")
{
    water water_props;

    CHECK_THROWS_AS(water_props.gibbsIAPWS(200.0, 1.0e6), std::domain_error);    // too cold
    CHECK_THROWS_AS(water_props.gibbsIAPWS(1200.0, 1.0e6), std::domain_error);   // region 5
    CHECK_THROWS_AS(water_props.gibbsIAPWS(400.0, 1500.0e6), std::domain_error); // above 1000 MPa
    CHECK_THROWS_AS(water_props.gibbsIAPWS(400.0, -1.0e6), std::domain_error);   // negative P
    CHECK_THROWS_AS(water::PsatIAPWS(700.0), std::domain_error);                 // above Tcrit

    // A failed call must not poison the (T, P) cache: a subsequent valid call
    // at the previously failing temperature has to be evaluated properly.
    water_props.gibbsIAPWS(400.0, 1.0e6);
    CHECK(water_props.region_ == 1);  // liquid: P > Psat(400 K) = 0.246 MPa
    CHECK(water_props.v_ > 0.0);
}

TEST_CASE("Test HKF standard-state water properties")
{
    hkf hkf_props;
    hkf_props.WaterProp(298.15, 1.0e5);

    CHECK_THAT(hkf_props.rhow_, Catch::Matchers::WithinAbs(997.047435408, 1.0e-6));
    CHECK_THAT(hkf_props.epsw_, Catch::Matchers::WithinAbs(78.2438908377, 1.0e-6));
    CHECK_THAT(hkf_props.G_, Catch::Matchers::WithinAbs(-237181.719094, 1.0e-6));
    CHECK_THAT(hkf_props.H_, Catch::Matchers::WithinAbs(-285830.819484, 1.0e-6));
    CHECK_THAT(hkf_props.S_, Catch::Matchers::WithinAbs(69.9280637281, 1.0e-6));
    CHECK_THAT(hkf_props.V_, Catch::Matchers::WithinAbs(1.80686287936e-05, 1.0e-12));
    CHECK_THAT(hkf_props.Cp_, Catch::Matchers::WithinAbs(75.3381005398, 1.0e-6));
}

TEST_CASE("Thermo table calculator matches water table-style output")
{
    ThermoTableCalculator calculator;
    const auto rows = calculator.evaluate("H2O", {25.0}, {1.0});

    REQUIRE(rows.size() == 1);
    const auto& row = rows.front();

    CHECK_THAT(row.T, Catch::Matchers::WithinAbs(25.0, abs_tolerance_));
    CHECK_THAT(row.P, Catch::Matchers::WithinAbs(1.0, abs_tolerance_));
    CHECK_THAT(row.rho, Catch::Matchers::WithinAbs(0.9970614, 2.0e-5));
    CHECK_THAT(row.logK, Catch::Matchers::WithinAbs(41.55238, 1.0e-3));
    CHECK_THAT(row.G, Catch::Matchers::WithinAbs(-237181.4, 1.0e1));
    CHECK_THAT(row.H, Catch::Matchers::WithinAbs(-285837.3, 1.0e1));
    CHECK_THAT(row.S, Catch::Matchers::WithinAbs(69.92418, 1.0e-2));
    CHECK_THAT(row.V, Catch::Matchers::WithinAbs(18.06830, 1.0e-2));
    CHECK_THAT(row.Cp, Catch::Matchers::WithinAbs(75.36053, 5.0e-2));
}

TEST_CASE("Thermo table calculator broadcasts vector inputs")
{
    ThermoTableCalculator calculator;
    const auto rows = calculator.evaluate("H2O", {25.0, 50.0}, {1.0});

    REQUIRE(rows.size() == 2);
    CHECK_THAT(rows[0].T, Catch::Matchers::WithinAbs(25.0, abs_tolerance_));
    CHECK_THAT(rows[1].T, Catch::Matchers::WithinAbs(50.0, abs_tolerance_));
    CHECK_THAT(rows[0].P, Catch::Matchers::WithinAbs(1.0, abs_tolerance_));
    CHECK_THAT(rows[1].P, Catch::Matchers::WithinAbs(1.0, abs_tolerance_));
    CHECK(rows[1].rho < rows[0].rho);
}

TEST_CASE("Thermo table calculator formats table output")
{
    ThermoTableCalculator calculator;
    const auto table = calculator.evaluateFormatted("H2O", {25.0, 50.0}, {1.0});

    CHECK(table.find("T") != std::string::npos);
    CHECK(table.find("P") != std::string::npos);
    CHECK(table.find("rho") != std::string::npos);
    CHECK(table.find("Cp") != std::string::npos);
    CHECK(table.find("   1   25.00") != std::string::npos);
    CHECK(table.find("   2   50.00") != std::string::npos);
    CHECK(table.find("75.33810") != std::string::npos);
    CHECK(table.find("75.29589") != std::string::npos);
}

TEST_CASE("Thermo table calculator exposes direct IAPWS water saturation pressure")
{
    ThermoTableCalculator calculator;
    const auto rows = calculator.evaluateWaterSaturationPressure({26.85, 226.85});

    REQUIRE(rows.size() == 2);
    CHECK_THAT(rows[0].T, Catch::Matchers::WithinAbs(26.85, abs_tolerance_));
    CHECK_THAT(rows[0].Psat, Catch::Matchers::WithinAbs(0.0353658941, 1.0e-8));
    CHECK_THAT(rows[1].Psat, Catch::Matchers::WithinAbs(26.3889776, 1.0e-6));

    const auto formatted = calculator.evaluateWaterSaturationPressureFormatted({25.0});
    CHECK(formatted.find("Psat") != std::string::npos);
}

TEST_CASE("Thermo table calculator reproduces HKF basis species reference data")
{
    ThermoTableCalculator calculator;
    const auto rows = calculator.evaluate("Na+", {25.0}, {1.0});

    REQUIRE(rows.size() == 1);
    const auto& row = rows.front();

    CHECK_THAT(row.G, Catch::Matchers::WithinAbs(-62591.0*UnitConversionFactors::cal2J_, abs_tolerance2_));
    CHECK_THAT(row.H, Catch::Matchers::WithinAbs(-57433.0*UnitConversionFactors::cal2J_, abs_tolerance2_));
    CHECK_THAT(row.S, Catch::Matchers::WithinAbs(13.96*UnitConversionFactors::cal2J_, abs_tolerance2_));
    // V from the HKF a-parameters and the Born Q term, in cm^3/mol (approx.
    // -1.1 in SUPCRT92; the difference comes from the Johnson-Norton
    // dielectric model used here). The pre-2025-07 reference value -1206.7
    // reflected a unit error (factor 1000) in hkf::ionProperties.
    CHECK_THAT(row.V, Catch::Matchers::WithinAbs(-1.2067361105, 1.0e-6));
    CHECK_THAT(row.Cp, Catch::Matchers::WithinAbs(38.1192278706, 1.0e-4));
}

TEST_CASE("Thermo table calculator reproduces HKF mineral reference data")
{
    ThermoTableCalculator calculator;
    const auto rows = calculator.evaluate("CALCITE", {25.0}, {1.0});

    REQUIRE(rows.size() == 1);
    const auto& row = rows.front();

    CHECK_THAT(row.G, Catch::Matchers::WithinAbs(-269880.0*UnitConversionFactors::cal2J_, abs_tolerance2_));
    CHECK_THAT(row.H, Catch::Matchers::WithinAbs(-288552.0*UnitConversionFactors::cal2J_, abs_tolerance2_));
    CHECK_THAT(row.S, Catch::Matchers::WithinAbs(22.15*UnitConversionFactors::cal2J_, abs_tolerance2_));
    CHECK_THAT(row.V, Catch::Matchers::WithinAbs(36.934, abs_tolerance2_));
    CHECK_THAT(row.Cp, Catch::Matchers::WithinAbs(81.876, 2.0e-1));
    CHECK_THAT(row.logK, Catch::Matchers::WithinAbs(1.84865, 1.0e-4));
}

TEST_CASE("Thermo table calculator reproduces HKF gas reference data")
{
    ThermoTableCalculator calculator;
    const auto rows = calculator.evaluate("CO2,g", {25.0}, {1.0});

    REQUIRE(rows.size() == 1);
    const auto& row = rows.front();

    CHECK_THAT(row.G, Catch::Matchers::WithinAbs(-94254.0*UnitConversionFactors::cal2J_, abs_tolerance2_));
    CHECK_THAT(row.H, Catch::Matchers::WithinAbs(-94051.0*UnitConversionFactors::cal2J_, abs_tolerance2_));
    CHECK_THAT(row.S, Catch::Matchers::WithinAbs(51.085*UnitConversionFactors::cal2J_, abs_tolerance2_));
    CHECK_THAT(row.Cp, Catch::Matchers::WithinAbs(37.155, 2.0e-1));
    CHECK_THAT(row.logK, Catch::Matchers::WithinAbs(-7.81373, 1.0e-4));
}

TEST_CASE("Thermo table calculator derives H2O,g logK from IAPWS saturation pressure")
{
    ThermoTableCalculator calculator;
    const auto rows = calculator.evaluate("H2O,g", {25.0}, {1.0});

    REQUIRE(rows.size() == 1);
    const auto& row = rows.front();

    const double expected_logK = std::log10(PhysicalConstants::standard_gas_pressure / water::PsatIAPWS(298.15));
    CHECK_THAT(row.logK, Catch::Matchers::WithinAbs(expected_logK, 1.0e-10));
    CHECK_THAT(row.logK, Catch::Matchers::WithinAbs(1.49897, 1.0e-4));
}

TEST_CASE("Thermo table standalone solver reads debug input")
{
    const std::string input = R"(SPECIESLIST
H2O
CALCITE
/ end
TEMPS
25.0 50.0
/ end
PRESSURES
1.0e5
/ end
TempUnit C
PresUnit Pa
IncludeIndex 1
)";

    std::istringstream input_stream(input);
    ThermoTableSolver solver;
    const auto output = solver.solve("thermotable_debug", input_stream);

    CHECK(output.find("SPECIES H2O") != std::string::npos);
    CHECK(output.find("SPECIES CALCITE") != std::string::npos);
    CHECK(output.find("IAPWS97_PSAT_H2O") != std::string::npos);
    CHECK(output.find("   1   25.00") != std::string::npos);
    CHECK(output.find("   2   50.00") != std::string::npos);
    CHECK(output.find("18.06863") != std::string::npos);
    CHECK(output.find("36.93400") != std::string::npos);
}

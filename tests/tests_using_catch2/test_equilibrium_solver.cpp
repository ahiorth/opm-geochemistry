#if __has_include(<catch2/catch.hpp>)
#include <catch2/catch.hpp>
#else
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#endif

#include <fstream>
#include <string>
#include <filesystem>
#include <opm/simulators/geochemistry/GeoChemIF.h>

#include "reference_results_eqsolver/EquilibriumTestCases/SetupAndRunEqsolver.hpp"
#include "reference_results_eqsolver/EquilibriumTestCases/RegisterEquilibriumTests.hpp"


/* Reads reference simulation data from a binary file. */
BasVecInfo getReferenceResults(GeoChemTestCase* testcase)
{
    const std::string archive_file = testcase->name() + ".cereal";

    BasVecInfo info_backup;
    std::ifstream cereal_infile(archive_file, std::ios::binary);
    cereal::BinaryInputArchive archive_infile(cereal_infile);
    archive_infile(info_backup);

    return info_backup;
}


TEST_CASE( "Test equilibrium solver") {

    static constexpr double abs_tolerance = 1.0e-5;
    // static constexpr double rel_tolerance = 1.0e-3;

    auto testFactory = EquilibriumTestCaseFactory();
    const auto& all_tests = testFactory.getAllEquilibriumTestCases();

    for(std::size_t i=0; i < all_tests.size(); ++i)
    {
        GeoChemTestCase* current_test = all_tests[i].get();

        std::cout << current_test->name() << "\n";
        std::filesystem::path cwd = std::filesystem::current_path();
        std::cout << "Current working directory: " << cwd << '\n';
        const auto reference_results = getReferenceResults(current_test);
        const auto eq_sol_results = gen_results_from_equilibrium_calculation(current_test);

        const auto reference_solution_props = reference_results.key_solution_properties_;
        for(auto const& [key, value]: eq_sol_results.key_solution_properties_){
            auto it = reference_solution_props.find(key);
            if(it != reference_solution_props.end()){
                const auto expected_value = it->second;
                CHECK_THAT(value, Catch::Matchers::WithinAbs(expected_value, abs_tolerance));
                // CHECK_THAT(value, Catch::Matchers::WithinRel(expected_value, rel_tolerance));
            }
        }

        const auto reference_species_data = reference_results.species_properties_;
        for(auto const& [key, value]: eq_sol_results.species_properties_){
            auto it = reference_species_data.find(key);
            if(it != reference_species_data.end()){
                const auto expected_value = it->second;
                if(std::fabs(value-expected_value) > abs_tolerance){
                    std::cout << "Key=(" << key.first << ", " << key.second << ")";
                    std::cout << ", Value=" << value << ", expected=" << expected_value << ".\n";
                }
                CHECK_THAT(value, Catch::Matchers::WithinAbs(expected_value, abs_tolerance));
                // CHECK_THAT(value, Catch::Matchers::WithinRel(expected_value, rel_tolerance));
            }
        }
    }

}

/* NaCl solution for the density tests below. */
struct NaClDensityCase: public GeoChemTestCase
{
    explicit NaClDensityCase(double molality, double T_celsius=25.0, double P_pascal=1.0e5)
    : molality_(molality)
    , T_celsius_(T_celsius)
    , P_pascal_(P_pascal)
    {}

    std::string name() const override
    {
        return "nacl_density";
    }

    std::string getInputFileAsString() const override
    {
        std::ostringstream input;
        input << "Temp " << T_celsius_ << "\n";
        input << "Pres " << P_pascal_ << "\n";
        input << "equilibrate 1\n\ngeochem\n\nsolution 0\npH 7\n";
        input << "Na " << molality_ << "\n";
        input << "Cl " << molality_ << "\n";
        input << "/end\n\n/end";
        return input.str();
    }

    double molality_;
    double T_celsius_;
    double P_pascal_;
};

double solutionDensity(GeoChemTestCase* testcase)
{
    const auto results = gen_results_from_equilibrium_calculation(testcase);
    return results.key_solution_properties_.at("Solution_density");
}

TEST_CASE("Solution density at different salinities")
{
    // Dilute limit: the solution density reduces to the pure-water density
    {
        NaClDensityCase dilute(1.0e-6);
        const auto results = gen_results_from_equilibrium_calculation(&dilute);
        const double rho_w = results.key_solution_properties_.at("Water_density");
        const double rho = results.key_solution_properties_.at("Solution_density");

        CHECK_THAT(rho_w, Catch::Matchers::WithinRel(997.047, 1.0e-3));
        CHECK_THAT(rho, Catch::Matchers::WithinRel(rho_w, 1.0e-4));
    }

    // Measured NaCl solution densities at 25 C / 1 bar (CRC handbook). The
    // HKF standard molal volumes ignore the excess volume of mixing, so the
    // deviation grows with molality; the tolerances reflect that.
    {
        NaClDensityCase molal_1(1.0);
        CHECK_THAT(solutionDensity(&molal_1), Catch::Matchers::WithinRel(1036.0, 1.0e-2));
    }
    {
        NaClDensityCase molal_4(4.0);
        CHECK_THAT(solutionDensity(&molal_4), Catch::Matchers::WithinRel(1139.5, 3.0e-2));
    }

    // Elevated temperature and pressure (1 molal, 200 C, 500 bar): water
    // expands much more than the salt contribution, so the brine stays
    // denser than pure water while both drop well below the 25 C value.
    {
        NaClDensityCase hot(1.0, 200.0, 5.0e7);
        const auto results = gen_results_from_equilibrium_calculation(&hot);
        const double rho_w = results.key_solution_properties_.at("Water_density");
        const double rho = results.key_solution_properties_.at("Solution_density");

        CHECK(rho > rho_w);
        CHECK(rho < 1000.0);
        CHECK(rho > 900.0);
    }
}

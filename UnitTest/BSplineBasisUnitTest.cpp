#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>

#include <libnurbs/Basis/BSplineBasis.hpp>

#include <stdexcept>

using namespace Catch;
using namespace libnurbs;
using namespace std;

TEST_CASE("Core/BSplineBasis", "[libnurbs_BSplineBasis]")
{
    SECTION("Evaluate (p=2)")
    {
        int degree = 2;
        KnotVector U{{0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0}};
        VecX result_0_0 = BSplineBasis::Evaluate(degree, U, 0.0);
        REQUIRE(result_0_0.size() == 3);
        REQUIRE(result_0_0(0) == Approx(1.0));
        REQUIRE(result_0_0(1) == Approx(0.0));
        REQUIRE(result_0_0(2) == Approx(0.0));
        VecX result_1_0 = BSplineBasis::Evaluate(degree, U, 1.0);
        REQUIRE(result_1_0.size() == 3);
        REQUIRE(result_1_0(0) == Approx(0.0));
        REQUIRE(result_1_0(1) == Approx(0.0));
        REQUIRE(result_1_0(2) == Approx(1.0));
        VecX result_0_5 = BSplineBasis::Evaluate(degree, U, 0.5);
        REQUIRE(result_0_5.size() == 3);
        REQUIRE(result_0_5(0) == Approx(0.5));
        REQUIRE(result_0_5(1) == Approx(0.5));
        REQUIRE(result_0_5(2) == Approx(0.0));
    }
}
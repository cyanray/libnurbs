#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>

#include <libnurbs/Basis/BSplineBasis.hpp>
#include <libnurbs/Core/KnotVector.hpp>

#include <stdexcept>

using namespace Catch;
using namespace libnurbs;
using namespace std;

TEST_CASE("Basis/BSplineBasis", "[libnurbs_BSplineBasis]")
{
    SECTION("Evaluate (p=2)")
    {
        int degree = 2;
        KnotVector U{{0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0}};
        {
            Numeric x = 0.0;
            VecX result = BSplineBasis::Evaluate(degree, U.Values(), U.FindSpanIndex(degree, x), x);
            REQUIRE(result.size() == 3);
            REQUIRE(result(0) == Approx(1.0));
            REQUIRE(result(1) == Approx(0.0));
            REQUIRE(result(2) == Approx(0.0));
        }
        {
            Numeric x = 1.0;
            VecX result = BSplineBasis::Evaluate(degree, U.Values(), U.FindSpanIndex(degree, x), x);
            REQUIRE(result.size() == 3);
            REQUIRE(result(0) == Approx(0.0));
            REQUIRE(result(1) == Approx(0.0));
            REQUIRE(result(2) == Approx(1.0));
        }
        {
            Numeric x = 0.5;
            VecX result = BSplineBasis::Evaluate(degree, U.Values(), U.FindSpanIndex(degree, x), x);
            REQUIRE(result.size() == 3);
            REQUIRE(result(0) == Approx(0.5));
            REQUIRE(result(1) == Approx(0.5));
            REQUIRE(result(2) == Approx(0.0));
        }
    }
}
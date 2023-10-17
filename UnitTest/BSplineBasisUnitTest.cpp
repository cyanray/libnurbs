#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>

#include <libnurbs/Basis/BSplineBasis.hpp>
#include <libnurbs/Core/KnotVector.hpp>

#include <stdexcept>
#include <sstream>

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
            VecX result = BSplineBasis::Evaluate(degree, U, x);
            REQUIRE(result.size() == 3);
            REQUIRE(result(0) == Approx(1.0));
            REQUIRE(result(1) == Approx(0.0));
            REQUIRE(result(2) == Approx(0.0));
        }
        {
            Numeric x = 1.0;
            VecX result = BSplineBasis::Evaluate(degree, U, x);
            REQUIRE(result.size() == 3);
            REQUIRE(result(0) == Approx(0.0));
            REQUIRE(result(1) == Approx(0.0));
            REQUIRE(result(2) == Approx(1.0));
        }
        {
            Numeric x = 0.5;
            VecX result = BSplineBasis::Evaluate(degree, U, x);
            REQUIRE(result.size() == 3);
            REQUIRE(result(0) == Approx(0.5));
            REQUIRE(result(1) == Approx(0.5));
            REQUIRE(result(2) == Approx(0.0));
        }
    }

    SECTION("derivative")
    {
        Numeric x1 = -1.0 / std::sqrt(3), x2 = -x1;
        auto conv = [](Numeric x) { return 0.25 * x + 0.25; };
        Numeric k1 = conv(x1), k2 = conv(x2);
        KnotVector U{{0, 0, 0, 0.5, 1, 1, 1}};
        VecX dN1 = BSplineBasis::EvaluateDerivative(2, U, k1, 1).row(1);
        VecX dN2 = BSplineBasis::EvaluateDerivative(2, U, k2, 1).row(1);
        MatX tN1 = dN1 * dN1.transpose();
        MatX tN2 = dN2 * dN2.transpose();
        MatX result = (tN1 + tN2) / 8.0;
        REQUIRE(result(0, 0) == Approx(1.333333));
        REQUIRE(result(0, 1) == Approx(-1));
        REQUIRE(result(0, 2) == Approx(-0.333333));
        REQUIRE(result(1, 0) == Approx(-1));
        REQUIRE(result(1, 1) == Approx(1));
        REQUIRE(result(1, 2) == Approx(0));
        REQUIRE(result(2, 0) == Approx(-0.333333));
        REQUIRE(result(2, 1) == Approx(0));
        REQUIRE(result(2, 2) == Approx(0.333333));
    }

}
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


TEST_CASE("BSplineBasis/Evaluate (p=2)", "[basis, evaluate]")
{
    int degree = 2;
    KnotVector U{{0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0}};
    SECTION("x = 0.0")
    {
        Numeric x = 0.0;
        VecX result = BSplineBasis::Evaluate(degree, U, x);
        REQUIRE(result.size() == 3);
        REQUIRE(result(0) == Approx(1.0));
        REQUIRE(result(1) == Approx(0.0));
        REQUIRE(result(2) == Approx(0.0));
    }

    SECTION("x = 1.0")
    {
        Numeric x = 1.0;
        VecX result = BSplineBasis::Evaluate(degree, U, x);
        REQUIRE(result.size() == 3);
        REQUIRE(result(0) == Approx(0.0));
        REQUIRE(result(1) == Approx(0.0));
        REQUIRE(result(2) == Approx(1.0));
    }

    SECTION("x = 0.5")
    {
        Numeric x = 0.5;
        VecX result = BSplineBasis::Evaluate(degree, U, x);
        REQUIRE(result.size() == 3);
        REQUIRE(result(0) == Approx(0.5));
        REQUIRE(result(1) == Approx(0.5));
        REQUIRE(result(2) == Approx(0.0));
    }
}


TEST_CASE("BSplineBasis/Derivative (p=2)", "[basis, derivative]")
{
    Numeric x1 = -1.0 / std::sqrt(3), x2 = -x1;
    auto conv = [](Numeric x) { return 0.25 * x + 0.25; };
    Numeric k1 = conv(x1), k2 = conv(x2);
    KnotVector U{{0, 0, 0, 0.5, 1, 1, 1}};
    VecX dN1 = BSplineBasis::EvaluateDerivative(2, U, k1, 1);
    VecX dN2 = BSplineBasis::EvaluateDerivative(2, U, k2, 1);
    MatX tN1 = dN1 * dN1.transpose();
    MatX tN2 = dN2 * dN2.transpose();
    MatX result = (tN1 + tN2) / 8.0;

    SECTION("u = 0 & v = 0, 1, 2")
    {
        REQUIRE(result(0, 0) == Approx(1.333333));
        REQUIRE(result(0, 1) == Approx(-1));
        REQUIRE(result(0, 2) == Approx(-0.333333));
    }

    SECTION("u = 1 & v = 0, 1, 2")
    {
        REQUIRE(result(1, 0) == Approx(-1));
        REQUIRE(result(1, 1) == Approx(1));
        REQUIRE(result(1, 2) == Approx(0));
    }

    SECTION("u = 2 & v = 0, 1, 2")
    {
        REQUIRE(result(2, 0) == Approx(-0.333333));
        REQUIRE(result(2, 1) == Approx(0));
        REQUIRE(result(2, 2) == Approx(0.333333));
    }
}

TEST_CASE("BSplineBasis/EvaluateAll above degree", "[basis][jet]")
{
    KnotVector knots{{0, 0, 0, 0.5, 1, 1, 1}};
    const Numeric x = 0.37;
    const MatX derivatives = BSplineBasis::EvaluateAll(2, knots, x, 4);

    REQUIRE(derivatives.rows() == 5);
    REQUIRE(derivatives.cols() == 3);
    for (int order = 0; order <= 2; ++order)
    {
        const VecX expected = order == 0
            ? BSplineBasis::Evaluate(2, knots, x)
            : BSplineBasis::EvaluateDerivative(2, knots, x, order);
        REQUIRE((derivatives.row(order).transpose() - expected).norm() == Approx(0.0).margin(1e-12));
    }
    REQUIRE(derivatives.row(3).norm() == Approx(0.0).margin(1e-12));
    REQUIRE(derivatives.row(4).norm() == Approx(0.0).margin(1e-12));
}

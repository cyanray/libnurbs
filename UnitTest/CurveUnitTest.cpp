#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>

#include <libnurbs/Curve/Curve.hpp>

#include <stdexcept>

using namespace Catch;
using namespace libnurbs;
using namespace std;

TEST_CASE("Core/Curve", "[libnurbs_BSplineCurve]")
{
    SECTION("Evaluate (p=2)")
    {
        int degree = 2;
        KnotVector U{{0.0, 0.0, 0.0, 1.0, 1.0, 1.0}};
        vector<Vec4> controlPoints{{0.0, 0.0, 0.0, 1.0},
                                   {1.0, 0.0, 0.0, 1.0},
                                   {1.0, 1.0, 0.0, 1.0}};
        Curve curve;
        curve.Degree = degree;
        curve.Knots = U;
        curve.ControlPoints = controlPoints;
        {
            Numeric x = 0.0;
            Vec3 result = curve.Evaluate(x);
            REQUIRE(result(0) == Approx(0.0));
            REQUIRE(result(1) == Approx(0.0));
            REQUIRE(result(2) == Approx(0.0));
        }
        {
            Numeric x = 1.0;
            Vec3 result = curve.Evaluate(x);
            REQUIRE(result(0) == Approx(1.0));
            REQUIRE(result(1) == Approx(1.0));
            REQUIRE(result(2) == Approx(0.0));
        }
        {
            Numeric x = 0.5;
            Vec3 result = curve.Evaluate(x);
            REQUIRE(result(0) == Approx(0.75));
            REQUIRE(result(1) == Approx(0.25));
            REQUIRE(result(2) == Approx(0.0));
        }
        {
            Numeric x = 0.3;
            Vec3 result = curve.Evaluate(x);
            REQUIRE(result(0) == Approx(0.51));
            REQUIRE(result(1) == Approx(0.09));
            REQUIRE(result(2) == Approx(0.0));
        }
        {
            Numeric x = 0.3;
            Vec3 result = curve.Evaluate(x);
            REQUIRE(result(0) == Approx(0.51));
            REQUIRE(result(1) == Approx(0.09));
            REQUIRE(result(2) == Approx(0.0));
        }
        {
            Numeric x = 0.678;
            Vec3 result = curve.Evaluate(x);
            REQUIRE(result(0) == Approx(0.896316));
            REQUIRE(result(1) == Approx(0.459684));
            REQUIRE(result(2) == Approx(0.0));
        }

    }
}
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>

#include <libnurbs/Surface/Surface.hpp>

#include <stdexcept>

using namespace Catch;
using namespace libnurbs;
using namespace std;


TEST_CASE("Surface/Evaluate", "[surface, evaluate]")
{
    ControlPointGrid grid;
    grid.UCount = 3;
    grid.VCount = 2;
    grid.Values.emplace_back(1, 0, 0, 1);
    grid.Values.emplace_back(1, 1, 0, 1);
    grid.Values.emplace_back(0, 1, 0, 2);
    grid.Values.emplace_back(1, 0, 1, 1);
    grid.Values.emplace_back(1, 1, 1, 1);
    grid.Values.emplace_back(0, 1, 1, 2);

    Surface surface;
    surface.DegreeU = 2;
    surface.DegreeV = 1;
    surface.KnotsU = KnotVector{{0.0, 0.0, 0.0, 1.0, 1.0, 1.0}};
    surface.KnotsV = KnotVector{{0.0, 0.0, 1.0, 1.0}};
    surface.ControlPoints = grid;

    SECTION("u = 0.5 & v = 0.5")
    {
        auto value = surface.Evaluate(0.5, 0.5);
        INFO("value:" << value.transpose());
        REQUIRE(value.x() == Approx(0.6));
        REQUIRE(value.y() == Approx(0.8));
        REQUIRE(value.z() == Approx(0.5));
    }
}

TEST_CASE("Surface/EvaluateDerivative", "[surface, derivative]")
{
    ControlPointGrid grid;
    grid.UCount = 3;
    grid.VCount = 2;
    grid.Values.emplace_back(1, 0, 0, 1);
    grid.Values.emplace_back(1, 1, 0, 1);
    grid.Values.emplace_back(0, 1, 0, 2);
    grid.Values.emplace_back(1, 0, 1, 1);
    grid.Values.emplace_back(1, 1, 1, 1);
    grid.Values.emplace_back(0, 1, 1, 2);

    Surface surface;
    surface.DegreeU = 2;
    surface.DegreeV = 1;
    surface.KnotsU = KnotVector{{0.0, 0.0, 0.0, 1.0, 1.0, 1.0}};
    surface.KnotsV = KnotVector{{0.0, 0.0, 1.0, 1.0}};
    surface.ControlPoints = grid;

    SECTION("du for (0.0, 0.0)")
    {
        auto value = surface.EvaluateDerivative(0.0, 0.0, 1, 0);
        INFO("value:" << value.transpose());
        REQUIRE(value.x() == Approx(0.0));
        REQUIRE(value.y() == Approx(2.0));
        REQUIRE(value.z() == Approx(0.0));
    }

    SECTION("dv for (0.0, 0.0)")
    {
        auto value = surface.EvaluateDerivative(0.0, 0.0, 0, 1);
        INFO("value:" << value.transpose());
        REQUIRE(value.x() == Approx(0.0));
        REQUIRE(value.y() == Approx(0.0));
        REQUIRE(value.z() == Approx(1.0));
    }

    SECTION("dudv for (0.0, 0.0)")
    {
        auto value = surface.EvaluateDerivative(0.0, 0.0, 1, 1);
        INFO("value:" << value.transpose());
        REQUIRE(value.x() == Approx(0.0));
        REQUIRE(value.y() == Approx(0.0));
        REQUIRE(value.z() == Approx(0.0));
    }

    SECTION("du for (0.5, 0.5)")
    {
        auto value = surface.EvaluateDerivative(0.5, 0.5, 1, 0);
        INFO("value:" << value.transpose());
        REQUIRE(value.x() == Approx(-1.0));
        REQUIRE(value.y() == Approx(1.0));
        REQUIRE(value.z() == Approx(0.0));
    }

    SECTION("dv for (0.5, 0.5)")
    {
        auto value = surface.EvaluateDerivative(0.5, 0.5, 0, 1);
        INFO("value:" << value.transpose());
        REQUIRE(value.x() == Approx(0.0));
        REQUIRE(value.y() == Approx(0.0));
        REQUIRE(value.z() == Approx(1.0));
    }

    SECTION("du for (1.0, 1.0)")
    {
        auto value = surface.EvaluateDerivative(1.0, 1.0, 1, 0);
        INFO("value:" << value.transpose());
        REQUIRE(value.x() == Approx(-2.0));
        REQUIRE(value.y() == Approx(0.0));
        REQUIRE(value.z() == Approx(0.0));
    }

    SECTION("dv for (1.0, 1.0)")
    {
        auto value = surface.EvaluateDerivative(1.0, 1.0, 0, 1);
        INFO("value:" << value.transpose());
        REQUIRE(value.x() == Approx(0.0));
        REQUIRE(value.y() == Approx(0.0));
        REQUIRE(value.z() == Approx(1.0));
    }


}

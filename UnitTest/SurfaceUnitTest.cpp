#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>

#include <libnurbs/Surface/Surface.hpp>

#include <stdexcept>

#include "libnurbs/Geometry/GeomRect.hpp"

using namespace Catch;
using namespace libnurbs;
using namespace std;


TEST_CASE("Surface/Evaluate", "[surface][evaluate]")
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

TEST_CASE("Surface/EvaluateDerivative", "[surface][evaluate_derivative]")
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


TEST_CASE("Surface/SearchParameter", "[surface][search_parameter]")
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

    SECTION("u = 0.0 & v = 0.0")
    {
        auto [u, v] = surface.SearchParameter({1.0, 0.0, 0.0});
        REQUIRE(std::abs(u) == Approx(0.0));
        REQUIRE(std::abs(v) == Approx(0.0));
    }

    SECTION("u = 1.0 & v = 1.0")
    {
        auto [u, v] = surface.SearchParameter({0.0, 1.0, 1.0});
        REQUIRE(u == Approx(1.0));
        REQUIRE(v == Approx(1.0));
    }

    SECTION("u = 0.5 & v = 0.5")
    {
        auto [u, v] = surface.SearchParameter({0.6, 0.8, 0.5});
        REQUIRE(u == Approx(0.5));
        REQUIRE(v == Approx(0.5));
    }

    SECTION("u = 0.00123 & v = 0.123")
    {
        auto [u, v] = surface.SearchParameter(surface.Evaluate(0.00123, 0.123));
        REQUIRE(u == Approx(0.00123));
        REQUIRE(v == Approx(0.123));
    }

    SECTION("u = 0.00123 & v = 0.999")
    {
        auto [u, v] = surface.SearchParameter(surface.Evaluate(0.00123, 0.999));
        REQUIRE(u == Approx(0.00123));
        REQUIRE(v == Approx(0.999));
    }
}


TEST_CASE("Surface/SearchParameter 2", "[surface][search_parameter]")
{
    auto geom_rect1 = GeomRect::Make({0, 0, 0}, {1, 0, 0},
                                     {0, 2, 0}, {1, 2, 0});
    geom_rect1.DegreeU = 3;
    geom_rect1.DegreeV = 3;
    geom_rect1.ControlPointCountU = 5;
    geom_rect1.ControlPointCountV = 5;

    auto geom_rect2 = GeomRect::Make({1, 0, 0}, {2.0, 0, 0},
                                     {1, 2, 0}, {2, 2, 0});
    geom_rect2.DegreeU = 3;
    geom_rect2.DegreeV = 3;
    geom_rect2.ControlPointCountU = 5;
    geom_rect2.ControlPointCountV = 5;

    auto rect1 = geom_rect1.GetSurface();
    auto rect2 = geom_rect2.GetSurface();

    const Numeric Zero = 1e-18;

    SECTION("0.0123")
    {
        Numeric v = 0.0123;
        auto point = rect1.Evaluate(1, v);
        auto [xi, eta] = rect2.SearchParameter(point);
        INFO("[xi, eta]: "<< xi << ", " << eta);
        REQUIRE(std::abs(xi - Zero) < 1e-16);
        REQUIRE(eta == Approx(v));
    }

    SECTION("0.33499526089999998")
    {
        Numeric v = 0.33499526089999998;
        auto point = rect1.Evaluate(1, v);
        auto [xi, eta] = rect2.SearchParameter(point);
        INFO("[xi, eta]: "<< xi << ", " << eta);
        REQUIRE(std::abs(xi - Zero) < 1e-16);
        REQUIRE(eta == Approx(v));
    }
}


TEST_CASE("Surface/SearchParameter 3", "[surface][search_parameter]")
{
    auto geom_rect1 = GeomRect::Make({0, 0, 0}, {1, 0, 0},
                                     {0, 2, 0}, {1, 2, 0});
    geom_rect1.DegreeU = 4;
    geom_rect1.DegreeV = 4;
    geom_rect1.ControlPointCountU = 10;
    geom_rect1.ControlPointCountV = 10;

    auto geom_rect2 = GeomRect::Make({1, 0, 0}, {2.0, 0, 0},
                                     {1, 2, 0}, {2, 2, 0});
    geom_rect2.DegreeU = 4;
    geom_rect2.DegreeV = 4;
    geom_rect2.ControlPointCountU = 10;
    geom_rect2.ControlPointCountV = 10;

    auto rect1 = geom_rect1.GetSurface();
    auto rect2 = geom_rect2.GetSurface();

    const Numeric Zero = 1e-18;

    SECTION("0.11166508696666666")
    {
        Numeric v = 0.11166508696666666;
        auto point = rect1.Evaluate(1, v);
        auto [xi, eta] = rect2.SearchParameter(point, 0, 0.11);
        INFO("[xi, eta]: "<< xi << ", " << eta);
        REQUIRE(std::abs(xi - Zero) < 1e-16);
        REQUIRE(eta == Approx(v));
    }

    SECTION("0.33499526089999998")
    {
        Numeric v = 0.33499526089999998;
        auto point = rect1.Evaluate(1, v);
        auto [xi, eta] = rect2.SearchParameter(point, 0, 0.33);
        INFO("[xi, eta]: "<< xi << ", " << eta);
        REQUIRE(std::abs(xi - Zero) < 1e-16);
        REQUIRE(eta == Approx(v));
    }
}


TEST_CASE("Surface/SearchParameter 4", "[surface][search_parameter]")
{
    auto geom_rect1 = GeomRect::Make({0, 0, 0}, {1, 0, 0},
                                     {0, 2, 0}, {1, 2, 0});
    geom_rect1.DegreeU = 4;
    geom_rect1.DegreeV = 4;
    geom_rect1.ControlPointCountU = 8;
    geom_rect1.ControlPointCountV = 9;

    auto geom_rect2 = GeomRect::Make({1, 0, 0}, {2.0, 0, 0},
                                     {1, 2, 0}, {2, 2, 0});
    geom_rect2.DegreeU = 4;
    geom_rect2.DegreeV = 4;
    geom_rect2.ControlPointCountU = 9;
    geom_rect2.ControlPointCountV = 10;

    auto rect1 = geom_rect1.GetSurface();
    auto rect2 = geom_rect2.GetSurface();

    const Numeric Zero = 1e-18;

    SECTION("1")
    {
        Numeric v = 0.58611363116000004;
        auto point1 = rect1.Evaluate(1, v);
        auto [xi, eta] = rect2.SearchParameter(point1, 0, v, 1E-12);
        auto point2 = rect2.Evaluate(xi, eta);
        INFO("[xi, eta]: "<< xi << ", " << eta);
        INFO("point 1: " << std::setprecision(10) << point1.transpose());
        INFO("point 2: " << std::setprecision(10) << point2.transpose());
        REQUIRE(std::abs(xi - Zero) < 1e-16);
        REQUIRE(point1.x() == Approx(point2.x()).epsilon(1e-10));
        REQUIRE(point1.y() == Approx(point2.y()).epsilon(1e-10));
        REQUIRE(point1.z() == Approx(point2.z()).epsilon(1e-10));
    }
}

TEST_CASE("Surface/InsertKnotU 1", "[surface][insert_knot]")
{
    auto geom_rect1 = GeomRect::Make({0, 0, 0}, {1, 0, 0},
                                     {0, 2, 0}, {1, 2, 0});
    geom_rect1.DegreeU = 3;
    geom_rect1.DegreeV = 3;
    geom_rect1.ControlPointCountU = 5;
    geom_rect1.ControlPointCountV = 5;

    auto surface = geom_rect1.GetSurface();

    for (Numeric knot: vector{0.1, 0.2, 0.345, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9})
    {
        auto new_surface = surface.InsertKnotU(knot);
        for (Numeric u: vector{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 1.0})
        {
            for (Numeric v: vector{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 1.0})
            {
                Vec3 value = surface.Evaluate(u, v);
                Vec3 new_value = new_surface.Evaluate(u, v);
                INFO("knot: "<< knot);
                INFO("u: " << u << ", v: " << v);
                INFO("value: " << value.transpose());
                INFO("new_value: " << new_value.transpose());
                REQUIRE(value.x() == Approx(new_value.x()));
                REQUIRE(value.y() == Approx(new_value.y()));
                REQUIRE(value.z() == Approx(new_value.z()));
            }
        }
    }
}

TEST_CASE("Surface/InsertKnotU 2", "[surface][insert_knot]")
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

    for (Numeric knot: vector{0.1, 0.2, 0.345, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9})
    {
        auto new_surface = surface.InsertKnotU(knot);
        for (Numeric u: vector{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 1.0})
        {
            for (Numeric v: vector{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 1.0})
            {
                Vec3 value = surface.Evaluate(u, v);
                Vec3 new_value = new_surface.Evaluate(u, v);
                INFO("knot: "<< knot);
                INFO("u: " << u << ", v: " << v);
                INFO("value: " << value.transpose());
                INFO("new_value: " << new_value.transpose());
                REQUIRE(value.x() == Approx(new_value.x()));
                REQUIRE(value.y() == Approx(new_value.y()));
                REQUIRE(value.z() == Approx(new_value.z()));
            }
        }
    }
}

TEST_CASE("Surface/InsertKnotU 3", "[surface][insert_knot]")
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

    for (Numeric knot1: vector{0.1, 0.2, 0.345, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9})
    {
        for (Numeric knot2: vector{0.1, 0.2, 0.345, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9})
        {
            auto new_surface = surface.InsertKnotU(knot1).InsertKnotU(knot2);
            for (Numeric u: vector{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 1.0})
            {
                for (Numeric v: vector{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 1.0})
                {
                    Vec3 value = surface.Evaluate(u, v);
                    Vec3 new_value = new_surface.Evaluate(u, v);
                    INFO("knot 1: "<< knot1);
                    INFO("knot 2: "<< knot2);
                    INFO("u: " << u << ", v: " << v);
                    INFO("value: " << value.transpose());
                    INFO("new_value: " << new_value.transpose());
                    REQUIRE(value.x() == Approx(new_value.x()));
                    REQUIRE(value.y() == Approx(new_value.y()));
                    REQUIRE(value.z() == Approx(new_value.z()));
                }
            }
        }
    }
}


TEST_CASE("Surface/InsertKnotV 1", "[surface][insert_knot]")
{
    auto geom_rect1 = GeomRect::Make({0, 0, 0}, {1, 0, 0},
                                     {0, 2, 0}, {1, 2, 0});
    geom_rect1.DegreeU = 3;
    geom_rect1.DegreeV = 3;
    geom_rect1.ControlPointCountU = 5;
    geom_rect1.ControlPointCountV = 5;

    auto surface = geom_rect1.GetSurface();

    for (Numeric knot: vector{0.1, 0.2, 0.345, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9})
    {
        auto new_surface = surface.InsertKnotV(knot);
        for (Numeric u: vector{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 1.0})
        {
            for (Numeric v: vector{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 1.0})
            {
                Vec3 value = surface.Evaluate(u, v);
                Vec3 new_value = new_surface.Evaluate(u, v);
                INFO("knot: "<< knot);
                INFO("u: " << u << ", v: " << v);
                INFO("value: " << value.transpose());
                INFO("new_value: " << new_value.transpose());
                REQUIRE(value.x() == Approx(new_value.x()));
                REQUIRE(value.y() == Approx(new_value.y()));
                REQUIRE(value.z() == Approx(new_value.z()));
            }
        }
    }
}

TEST_CASE("Surface/InsertKnotV 2", "[surface][insert_knot]")
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

    for (Numeric knot: vector{0.1, 0.2, 0.345, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9})
    {
        auto new_surface = surface.InsertKnotV(knot);
        for (Numeric u: vector{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 1.0})
        {
            for (Numeric v: vector{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 1.0})
            {
                Vec3 value = surface.Evaluate(u, v);
                Vec3 new_value = new_surface.Evaluate(u, v);
                INFO("knot: "<< knot);
                INFO("u: " << u << ", v: " << v);
                INFO("value: " << value.transpose());
                INFO("new_value: " << new_value.transpose());
                REQUIRE(value.x() == Approx(new_value.x()));
                REQUIRE(value.y() == Approx(new_value.y()));
                REQUIRE(value.z() == Approx(new_value.z()));
            }
        }
    }
}

TEST_CASE("Surface/InsertKnotV 3", "[surface][insert_knot]")
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

    for (Numeric knot1: vector{0.1, 0.2, 0.345, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9})
    {
        for (Numeric knot2: vector{0.1, 0.2, 0.345, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9})
        {
            auto new_surface = surface.InsertKnotV(knot1).InsertKnotV(knot2);
            for (Numeric u: vector{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 1.0})
            {
                for (Numeric v: vector{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 1.0})
                {
                    Vec3 value = surface.Evaluate(u, v);
                    Vec3 new_value = new_surface.Evaluate(u, v);
                    INFO("knot 1: "<< knot1);
                    INFO("knot 2: "<< knot2);
                    INFO("u: " << u << ", v: " << v);
                    INFO("value: " << value.transpose());
                    INFO("new_value: " << new_value.transpose());
                    REQUIRE(value.x() == Approx(new_value.x()));
                    REQUIRE(value.y() == Approx(new_value.y()));
                    REQUIRE(value.z() == Approx(new_value.z()));
                }
            }
        }
    }
}


TEST_CASE("Surface/InsertKnotUV", "[surface][insert_knot]")
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

    SECTION("UV")
    {
        for (Numeric knot1: vector{0.1, 0.2, 0.345, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9})
        {
            for (Numeric knot2: vector{0.1, 0.2, 0.345, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9})
            {
                auto new_surface = surface.InsertKnotU(knot1).InsertKnotV(knot2);
                for (Numeric u: vector{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 1.0})
                {
                    for (Numeric v: vector{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 1.0})
                    {
                        Vec3 value = surface.Evaluate(u, v);
                        Vec3 new_value = new_surface.Evaluate(u, v);
                        INFO("knot 1: "<< knot1);
                        INFO("knot 2: "<< knot2);
                        INFO("u: " << u << ", v: " << v);
                        INFO("value: " << value.transpose());
                        INFO("new_value: " << new_value.transpose());
                        REQUIRE(value.x() == Approx(new_value.x()));
                        REQUIRE(value.y() == Approx(new_value.y()));
                        REQUIRE(value.z() == Approx(new_value.z()));
                    }
                }
            }
        }
    }

    SECTION("VU")
    {
        for (Numeric knot1: vector{0.1, 0.2, 0.345, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9})
        {
            for (Numeric knot2: vector{0.1, 0.2, 0.345, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9})
            {
                auto new_surface = surface.InsertKnotV(knot1).InsertKnotU(knot2);
                for (Numeric u: vector{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 1.0})
                {
                    for (Numeric v: vector{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 1.0})
                    {
                        Vec3 value = surface.Evaluate(u, v);
                        Vec3 new_value = new_surface.Evaluate(u, v);
                        INFO("knot 1: "<< knot1);
                        INFO("knot 2: "<< knot2);
                        INFO("u: " << u << ", v: " << v);
                        INFO("value: " << value.transpose());
                        INFO("new_value: " << new_value.transpose());
                        REQUIRE(value.x() == Approx(new_value.x()));
                        REQUIRE(value.y() == Approx(new_value.y()));
                        REQUIRE(value.z() == Approx(new_value.z()));
                    }
                }
            }
        }
    }
}

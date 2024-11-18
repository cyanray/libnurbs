#include <stdexcept>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <libnurbs/Geometry/GeomRect.hpp>
#include <libnurbs/Surface/Surface.hpp>

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
        REQUIRE(value.x() == Approx(-1.28));
        REQUIRE(value.y() == Approx(0.96));
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
        REQUIRE(value.x() == Approx(-1.0));
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
    auto geom_rect1 = GeomRect::Make({0, 0, 0},
                                     {1, 0, 0},
                                     {0, 2, 0},
                                     {1, 2, 0});
    geom_rect1.DegreeU = 3;
    geom_rect1.DegreeV = 3;
    geom_rect1.ControlPointCountU = 5;
    geom_rect1.ControlPointCountV = 5;

    auto geom_rect2 = GeomRect::Make({1, 0, 0},
                                     {2.0, 0, 0},
                                     {1, 2, 0},
                                     {2, 2, 0});
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
    auto geom_rect1 = GeomRect::Make({0, 0, 0},
                                     {1, 0, 0},
                                     {0, 2, 0},
                                     {1, 2, 0});
    geom_rect1.DegreeU = 4;
    geom_rect1.DegreeV = 4;
    geom_rect1.ControlPointCountU = 10;
    geom_rect1.ControlPointCountV = 10;

    auto geom_rect2 = GeomRect::Make({1, 0, 0},
                                     {2.0, 0, 0},
                                     {1, 2, 0},
                                     {2, 2, 0});
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
    auto geom_rect1 = GeomRect::Make({0, 0, 0},
                                     {1, 0, 0},
                                     {0, 2, 0},
                                     {1, 2, 0});
    geom_rect1.DegreeU = 4;
    geom_rect1.DegreeV = 4;
    geom_rect1.ControlPointCountU = 8;
    geom_rect1.ControlPointCountV = 9;

    auto geom_rect2 = GeomRect::Make({1, 0, 0},
                                     {2.0, 0, 0},
                                     {1, 2, 0},
                                     {2, 2, 0});
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

TEST_CASE("Surface/SearchParameterOn", "[surface][search_parameter]")
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

    SECTION("u = 0.0 & v = 0.0 (a)")
    {
        auto [u, v] = surface.SearchParameterOn({1.0, 0.0, 0.0}, 0, 0);
        REQUIRE(std::abs(u) == Approx(0.0));
        REQUIRE(std::abs(v) == Approx(0.0));
    }

    SECTION("u = 0.0 & v = 0.0 (b)")
    {
        auto [u, v] = surface.SearchParameterOn({1.0, 0.0, 0.0}, 1, 0);
        REQUIRE(std::abs(u) == Approx(0.0));
        REQUIRE(std::abs(v) == Approx(0.0));
    }


    SECTION("u = 0.5 & v = 0.5")
    {
        auto [u, v] = surface.SearchParameterOn({0.6, 0.8, 0.5}, 0, 0.5);
        REQUIRE(u == Approx(0.5));
        REQUIRE(v == Approx(0.5));
    }

    SECTION("u = 0 & v = 0.123")
    {
        auto [u, v] = surface.SearchParameterOn(surface.Evaluate(0, 0.123), 0, 0);
        REQUIRE(u == Approx(0));
        REQUIRE(v == Approx(0.123));
    }

    SECTION("u = 0.123 & v = 0")
    {
        auto [u, v] = surface.SearchParameterOn(surface.Evaluate(0.123, 0), 1, 0);
        REQUIRE(u == Approx(0.123));
        REQUIRE(v == Approx(0));
    }
}

TEST_CASE("Surface/BinarySearchParameterOn", "[surface][search_parameter]")
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

    SECTION("u = 0.0 & v = 0.0 (a)")
    {
        auto [u, v] = surface.BinarySearchParameterOn({1.0, 0.0, 0.0}, 0, 0);
        REQUIRE(std::abs(u) < 1e-20);
        REQUIRE(std::abs(v) < 1e-20);
    }

    SECTION("u = 0.0 & v = 0.0 (b)")
    {
        auto [u, v] = surface.BinarySearchParameterOn({1.0, 0.0, 0.0}, 1, 0);
        REQUIRE(std::abs(u) < 1e-20);
        REQUIRE(std::abs(v) < 1e-20);
    }


    SECTION("u = 0.5 & v = 0.5")
    {
        auto [u, v] = surface.BinarySearchParameterOn({0.6, 0.8, 0.5}, 0, 0.5);
        REQUIRE(u == Approx(0.5));
        REQUIRE(v == Approx(0.5));
    }

    SECTION("u = 0 & v = 0.123")
    {
        auto [u, v] = surface.BinarySearchParameterOn(surface.Evaluate(0, 0.123), 0, 0);
        REQUIRE(u == Approx(0));
        REQUIRE(v == Approx(0.123));
    }

    SECTION("u = 0.123 & v = 0")
    {
        auto [u, v] = surface.BinarySearchParameterOn(surface.Evaluate(0.123, 0), 1, 0);
        REQUIRE(u == Approx(0.123));
        REQUIRE(v == Approx(0));
    }

    SECTION("u = 0.123 & v = 1")
    {
        auto [u, v] = surface.BinarySearchParameterOn(surface.Evaluate(0.123, 1), 1, 1);
        REQUIRE(u == Approx(0.123));
        REQUIRE(v == Approx(1));
    }
}

TEST_CASE("Surface/BinarySearchParameterOn 2", "[surface][search_parameter]")
{
    auto get_surf = []()
    {
        ControlPointGrid grid;
        grid.UCount = 8;
        grid.VCount = 8;
        grid.Values.emplace_back(-0.0294261835515e-3, 0.1506999969482, 0.0941875, 1);
        grid.Values.emplace_back(4.9402571349272e-6, 0.1521587224039, 0.0753500315527, 1);
        grid.Values.emplace_back(0.0576485078669e-3, 0.1517236898364, 0.0376750417598, 1);
        grid.Values.emplace_back(0.0481794104219e-3, 0.1514366598289, -0.0188374776592, 1);
        grid.Values.emplace_back(0.0190505260041e-3, 0.1322965328062, -0.075349982688, 1);
        grid.Values.emplace_back(0.2702707795319e-6, 0.0802878958362, -0.131862495488, 1);
        grid.Values.emplace_back(-0.0189456996064e-3, 0.0311588038952, -0.1695374986589, 1);
        grid.Values.emplace_back(-0.0292324908078e-3, 0, -0.188375, 1);
        grid.Values.emplace_back(-0.0294261835515e-3, 0.1306066640218, 0.0941875, 1);
        grid.Values.emplace_back(4.9402571349265e-6, 0.13187089275, 0.0753500315527, 1);
        grid.Values.emplace_back(0.0576485078669e-3, 0.1314938645249, 0.0376750417598, 1);
        grid.Values.emplace_back(0.0481794104219e-3, 0.1312451051851, -0.0188374776592, 1);
        grid.Values.emplace_back(0.0190505260041e-3, 0.1146569950987, -0.075349982688, 1);
        grid.Values.emplace_back(0.2702707795315e-6, 0.0695828430581, -0.131862495488, 1);
        grid.Values.emplace_back(-0.0189456996064e-3, 0.0270042967092, -0.1695374986589, 1);
        grid.Values.emplace_back(-0.0292324908078e-3, 0, -0.188375, 1);
        grid.Values.emplace_back(-0.0294261835515e-3, 0.0904199981689, 0.0941875, 1);
        grid.Values.emplace_back(4.9402571349249e-6, 0.0912952334423, 0.0753500315527, 1);
        grid.Values.emplace_back(0.0576485078669e-3, 0.0910342139019, 0.0376750417598, 1);
        grid.Values.emplace_back(0.0481794104219e-3, 0.0908619958974, -0.0188374776592, 1);
        grid.Values.emplace_back(0.0190505260041e-3, 0.0793779196837, -0.075349982688, 1);
        grid.Values.emplace_back(0.2702707795307e-6, 0.0481727375017, -0.131862495488, 1);
        grid.Values.emplace_back(-0.0189456996064e-3, 0.0186952823371, -0.1695374986589, 1);
        grid.Values.emplace_back(-0.0292324908078e-3, 0, -0.188375, 1);
        grid.Values.emplace_back(-0.0294261835515e-3, 0.0301399993896, 0.0941875, 1);
        grid.Values.emplace_back(4.9402571349226e-6, 0.0304317444808, 0.0753500315527, 1);
        grid.Values.emplace_back(0.0576485078669e-3, 0.0303447379673, 0.0376750417598, 1);
        grid.Values.emplace_back(0.0481794104219e-3, 0.0302873319658, -0.0188374776592, 1);
        grid.Values.emplace_back(0.0190505260041e-3, 0.0264593065612, -0.075349982688, 1);
        grid.Values.emplace_back(0.2702707795295e-6, 0.0160575791672, -0.131862495488, 1);
        grid.Values.emplace_back(-0.0189456996064e-3, 6.2317607790403e-3, -0.1695374986589, 1);
        grid.Values.emplace_back(-0.0292324908078e-3, 0, -0.188375, 1);
        grid.Values.emplace_back(-0.0294261835515e-3, -0.0301399993896, 0.0941875, 1);
        grid.Values.emplace_back(4.9402571349203e-6, -0.0304317444808, 0.0753500315527, 1);
        grid.Values.emplace_back(0.0576485078669e-3, -0.0303447379673, 0.0376750417598, 1);
        grid.Values.emplace_back(0.0481794104219e-3, -0.0302873319658, -0.0188374776592, 1);
        grid.Values.emplace_back(0.0190505260041e-3, -0.0264593065612, -0.075349982688, 1);
        grid.Values.emplace_back(0.2702707795283e-6, -0.0160575791672, -0.131862495488, 1);
        grid.Values.emplace_back(-0.0189456996064e-3, -6.2317607790402e-3, -0.1695374986589, 1);
        grid.Values.emplace_back(-0.0292324908078e-3, 0, -0.188375, 1);
        grid.Values.emplace_back(-0.0294261835515e-3, -0.0904199981689, 0.0941875, 1);
        grid.Values.emplace_back(4.9402571349179e-6, -0.0912952334423, 0.0753500315527, 1);
        grid.Values.emplace_back(0.0576485078669e-3, -0.0910342139019, 0.0376750417598, 1);
        grid.Values.emplace_back(0.0481794104219e-3, -0.0908619958974, -0.0188374776592, 1);
        grid.Values.emplace_back(0.0190505260041e-3, -0.0793779196837, -0.075349982688, 1);
        grid.Values.emplace_back(0.270270779527e-6, -0.0481727375017, -0.131862495488, 1);
        grid.Values.emplace_back(-0.0189456996064e-3, -0.0186952823371, -0.1695374986589, 1);
        grid.Values.emplace_back(-0.0292324908078e-3, 0, -0.188375, 1);
        grid.Values.emplace_back(-0.0294261835515e-3, -0.1306066640218, 0.0941875, 1);
        grid.Values.emplace_back(4.9402571349164e-6, -0.13187089275, 0.0753500315527, 1);
        grid.Values.emplace_back(0.0576485078669e-3, -0.1314938645249, 0.0376750417598, 1);
        grid.Values.emplace_back(0.0481794104219e-3, -0.1312451051851, -0.0188374776592, 1);
        grid.Values.emplace_back(0.0190505260041e-3, -0.1146569950987, -0.075349982688, 1);
        grid.Values.emplace_back(0.2702707795262e-6, -0.0695828430581, -0.131862495488, 1);
        grid.Values.emplace_back(-0.0189456996064e-3, -0.0270042967092, -0.1695374986589, 1);
        grid.Values.emplace_back(-0.0292324908078e-3, 0, -0.188375, 1);
        grid.Values.emplace_back(-0.0294261835515e-3, -0.1506999969482, 0.0941875, 1);
        grid.Values.emplace_back(4.9402571349156e-6, -0.1521587224039, 0.0753500315527, 1);
        grid.Values.emplace_back(0.0576485078669e-3, -0.1517236898364, 0.0376750417598, 1);
        grid.Values.emplace_back(0.0481794104219e-3, -0.1514366598289, -0.0188374776592, 1);
        grid.Values.emplace_back(0.0190505260041e-3, -0.1322965328062, -0.075349982688, 1);
        grid.Values.emplace_back(0.2702707795258e-6, -0.0802878958362, -0.131862495488, 1);
        grid.Values.emplace_back(-0.0189456996064e-3, -0.0311588038952, -0.1695374986589, 1);
        grid.Values.emplace_back(-0.0292324908078e-3, 0, -0.188375, 1);
        Surface surface;
        surface.DegreeU = 3;
        surface.DegreeV = 3;
        surface.KnotsU = KnotVector{{0, 0, 0, 0, 0.2, 0.4, 0.6, 0.8, 1, 1, 1, 1}};
        surface.KnotsV = KnotVector{
            {0, 0, 0, 0, 0.19999999999999998, 0.39999999999999997, 0.6, 0.7999999999999999, 1, 1, 1, 1}
        };
        surface.ControlPoints = grid;
        return surface;
    };

    Surface surface = get_surf();

    SECTION("u = 0.0 & v = 0.0 (a)")
    {
        auto point = surface.Evaluate(0, 0);
        auto [u, v] = surface.BinarySearchParameterOn(point, 0, 0, 1e-10);
        auto point_searched = surface.Evaluate(u, v);
        INFO("point: " << point.transpose());
        INFO("point_search: " << point_searched.transpose());
        REQUIRE(std::abs((point - point_searched).norm()) < 1e-10);
    }

    SECTION("u = 0.0 & v = 0.0 (b)")
    {
        auto point = surface.Evaluate(0, 0);
        auto [u, v] = surface.BinarySearchParameterOn(point, 1, 0, 1e-10);
        auto point_searched = surface.Evaluate(u, v);
        INFO("point: " << point.transpose());
        INFO("point_search: " << point_searched.transpose());
        REQUIRE(std::abs((point - point_searched).norm()) < 1e-10);
    }

    SECTION("u = 1.0 & v = 0.0 (a)")
    {
        auto point = surface.Evaluate(1, 0);
        auto [u, v] = surface.BinarySearchParameterOn(point, 0, 1, 1e-10);
        auto point_searched = surface.Evaluate(u, v);
        INFO("point: " << point.transpose());
        INFO("point_search: " << point_searched.transpose());
        REQUIRE(std::abs((point - point_searched).norm()) < 1e-10);
    }

    SECTION("u = 1.0 & v = 0.1234567")
    {
        auto point = surface.Evaluate(1, 0.1234567);
        auto [u, v] = surface.BinarySearchParameterOn(point, 0, 1, 1e-10);
        auto point_searched = surface.Evaluate(u, v);
        INFO("point: " << point.transpose());
        INFO("point_search: " << point_searched.transpose());
        REQUIRE(std::abs((point - point_searched).norm()) < 1e-10);
    }
}

TEST_CASE("Surface/InsertKnotU 1", "[surface][insert_knot]")
{
    auto geom_rect1 = GeomRect::Make({0, 0, 0},
                                     {1, 0, 0},
                                     {0, 2, 0},
                                     {1, 2, 0});
    geom_rect1.DegreeU = 3;
    geom_rect1.DegreeV = 3;
    geom_rect1.ControlPointCountU = 5;
    geom_rect1.ControlPointCountV = 5;

    auto surface = geom_rect1.GetSurface();

    for (Numeric knot : vector{0.1, 0.2, 0.345, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9})
    {
        auto new_surface = surface.InsertKnotU(knot);
        for (Numeric u : vector{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 1.0})
        {
            for (Numeric v : vector{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 1.0})
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

    for (Numeric knot : vector{0.1, 0.2, 0.345, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9})
    {
        auto new_surface = surface.InsertKnotU(knot);
        for (Numeric u : vector{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 1.0})
        {
            for (Numeric v : vector{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 1.0})
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

    for (Numeric knot1 : vector{0.1, 0.2, 0.345, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9})
    {
        for (Numeric knot2 : vector{0.1, 0.2, 0.345, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9})
        {
            auto new_surface = surface.InsertKnotU(knot1).InsertKnotU(knot2);
            for (Numeric u : vector{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 1.0})
            {
                for (Numeric v : vector{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 1.0})
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
    auto geom_rect1 = GeomRect::Make({0, 0, 0},
                                     {1, 0, 0},
                                     {0, 2, 0},
                                     {1, 2, 0});
    geom_rect1.DegreeU = 3;
    geom_rect1.DegreeV = 3;
    geom_rect1.ControlPointCountU = 5;
    geom_rect1.ControlPointCountV = 5;

    auto surface = geom_rect1.GetSurface();

    for (Numeric knot : vector{0.1, 0.2, 0.345, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9})
    {
        auto new_surface = surface.InsertKnotV(knot);
        for (Numeric u : vector{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 1.0})
        {
            for (Numeric v : vector{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 1.0})
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

    for (Numeric knot : vector{0.1, 0.2, 0.345, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9})
    {
        auto new_surface = surface.InsertKnotV(knot);
        for (Numeric u : vector{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 1.0})
        {
            for (Numeric v : vector{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 1.0})
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

    for (Numeric knot1 : vector{0.1, 0.2, 0.345, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9})
    {
        for (Numeric knot2 : vector{0.1, 0.2, 0.345, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9})
        {
            auto new_surface = surface.InsertKnotV(knot1).InsertKnotV(knot2);
            for (Numeric u : vector{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 1.0})
            {
                for (Numeric v : vector{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 1.0})
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
        for (Numeric knot1 : vector{0.1, 0.2, 0.345, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9})
        {
            for (Numeric knot2 : vector{0.1, 0.2, 0.345, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9})
            {
                auto new_surface = surface.InsertKnotU(knot1).InsertKnotV(knot2);
                for (Numeric u : vector{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 1.0})
                {
                    for (Numeric v : vector{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 1.0})
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
        for (Numeric knot1 : vector{0.1, 0.2, 0.345, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9})
        {
            for (Numeric knot2 : vector{0.1, 0.2, 0.345, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9})
            {
                auto new_surface = surface.InsertKnotV(knot1).InsertKnotU(knot2);
                for (Numeric u : vector{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 1.0})
                {
                    for (Numeric v : vector{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 1.0})
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


TEST_CASE("Surface/RemoveKnot", "[surface][remove_knot]")
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

    SECTION("times=1,1")
    {
        for (Numeric knot : vector{0.1, 0.2, 0.345, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9})
        {
            auto [new_surface, t] = surface.InsertKnotU(knot, 1).RemoveKnotU(knot, 1);
            REQUIRE(t == 1);
            for (Numeric u : vector{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 1.0})
            {
                for (Numeric v : vector{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 1.0})
                {
                    Vec3 value = surface.Evaluate(u, v);
                    Vec3 new_value = new_surface.Evaluate(u, v);
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

    SECTION("times=2,1")
    {
        for (Numeric knot : vector{0.1, 0.2, 0.345, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9})
        {
            auto [new_surface, t] = surface.InsertKnotU(knot, 2).RemoveKnotU(knot, 1);
            REQUIRE(t == 1);
            for (Numeric u : vector{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 1.0})
            {
                for (Numeric v : vector{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 1.0})
                {
                    Vec3 value = surface.Evaluate(u, v);
                    Vec3 new_value = new_surface.Evaluate(u, v);
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

    SECTION("times=2,2")
    {
        for (Numeric knot : vector{0.1, 0.2, 0.345, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9})
        {
            auto [new_surface, t] = surface.InsertKnotU(knot, 2).RemoveKnotU(knot, 2);
            REQUIRE(t == 2);
            for (Numeric u : vector{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 1.0})
            {
                for (Numeric v : vector{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 1.0})
                {
                    Vec3 value = surface.Evaluate(u, v);
                    Vec3 new_value = new_surface.Evaluate(u, v);
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


TEST_CASE("Surface/ElevateDegree", "[surface][elevate_degree]")
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

    SECTION("U")
    {
        auto new_surface = surface.ElevateDegreeU();
        for (Numeric u : vector{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 1.0})
        {
            for (Numeric v : vector{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 1.0})
            {
                Vec3 value = surface.Evaluate(u, v);
                Vec3 new_value = new_surface.Evaluate(u, v);
                INFO("u: " << u << ", v: " << v);
                INFO("value: " << value.transpose());
                INFO("new_value: " << new_value.transpose());
                REQUIRE(value.x() == Approx(new_value.x()));
                REQUIRE(value.y() == Approx(new_value.y()));
                REQUIRE(value.z() == Approx(new_value.z()));
            }
        }
    }

    SECTION("V")
    {
        auto new_surface = surface.ElevateDegreeV();
        for (Numeric u : vector{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 1.0})
        {
            for (Numeric v : vector{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 1.0})
            {
                Vec3 value = surface.Evaluate(u, v);
                Vec3 new_value = new_surface.Evaluate(u, v);
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


TEST_CASE("Surface::LoadFromFile - Valid TXT Input", "[LoadFromFile]")
{
    std::string valid_input = R"(LIBNURBS TXT SURFACE
*DegreeU: 3
*DegreeV: 2
*KnotsU:
0, 0, 0, 1, 2, 3, 3, 3
*KnotsV:
0, 0, 1, 1
*ControlPoints:
4 3
1.0 2.0 3.0 1.0
4.0 5.0 6.0 1.0
7.0 8.0 9.0 1.0
1.1 2.1 3.1 1.0
4.1 5.1 6.1 1.0
7.1 8.1 9.1 1.0
1.2 2.2 3.2 1.0
4.2 5.2 6.2 1.0
7.2 8.2 9.2 1.0
1.3 2.3 3.3 1.0
4.3 5.3 6.3 1.0
7.3 8.3 9.3 1.0
)";
    std::istringstream iss(valid_input);
    Surface surface;
    REQUIRE_NOTHROW(surface.LoadFromFile(iss));

    // Check Degrees
    REQUIRE(surface.DegreeU == 3);
    REQUIRE(surface.DegreeV == 2);

    // Check Knot Vectors
    REQUIRE(surface.KnotsU.Values().size() == 8);
    REQUIRE(surface.KnotsV.Values().size() == 4);

    // Check Control Points Grid Dimensions
    REQUIRE(surface.ControlPoints.UCount == 4);
    REQUIRE(surface.ControlPoints.VCount == 3);

    // Check Control Points Values
    REQUIRE(surface.ControlPoints.Values.size() == 12);
    REQUIRE(surface.ControlPoints.Get(0, 0).x() == Approx(1.0));
    REQUIRE(surface.ControlPoints.Get(3, 2).z() == Approx(9.3));
}

TEST_CASE("Surface::LoadFromFile - Missing Header", "[LoadFromFile]")
{
    std::string missing_header = R"(*DegreeU: 3
*DegreeV: 2
)";
    std::istringstream iss(missing_header);
    Surface surface;
    REQUIRE_THROWS_WITH(surface.LoadFromFile(iss), Catch::Matchers::ContainsSubstring("Invalid file magic"));
}

TEST_CASE("Surface::LoadFromFile - Invalid Magic Number", "[LoadFromFile]")
{
    std::string invalid_magic = R"(INVALIDMAGIC TXT SURFACE
*DegreeU: 3
*DegreeV: 2
)";
    std::istringstream iss(invalid_magic);
    Surface surface;
    REQUIRE_THROWS_WITH(surface.LoadFromFile(iss), Catch::Matchers::ContainsSubstring("Invalid file magic"));
}

TEST_CASE("Surface::LoadFromFile - Unknown Serialization Mode", "[LoadFromFile]")
{
    std::string unknown_mode = R"(LIBNURBS UNKNOWN SURFACE
*DegreeU: 3
*DegreeV: 2
)";
    std::istringstream iss(unknown_mode);
    Surface surface;
    REQUIRE_THROWS_WITH(surface.LoadFromFile(iss), Catch::Matchers::ContainsSubstring("Unknown serialization mode"));
}

TEST_CASE("Surface::LoadFromFile - Type Mismatch", "[LoadFromFile]")
{
    std::string type_mismatch = R"(LIBNURBS TXT CURVE
*DegreeU: 3
*DegreeV: 2
)";
    std::istringstream iss(type_mismatch);
    Surface surface;
    REQUIRE_THROWS_WITH(surface.LoadFromFile(iss), Catch::Matchers::ContainsSubstring("Object type mismatch"));
}

TEST_CASE("Surface::LoadFromFile - Invalid Degrees", "[LoadFromFile]")
{
    std::string invalid_degrees = R"(LIBNURBS TXT SURFACE
*DegreeU: invalid
*DegreeV: 2
)";
    std::istringstream iss(invalid_degrees);
    Surface surface;
    REQUIRE_THROWS_WITH(surface.LoadFromFile(iss), Catch::Matchers::ContainsSubstring("Failed to parse DegreeU value"));
}

TEST_CASE("Surface::LoadFromFile - Invalid Knot Values", "[LoadFromFile]")
{
    std::string invalid_knots = R"(LIBNURBS TXT SURFACE
*DegreeU: 3
*DegreeV: 2
*KnotsU:
0, 0, invalid, 1, 2, 3, 3, 3
)";
    std::istringstream iss(invalid_knots);
    Surface surface;
    REQUIRE_THROWS_WITH(surface.LoadFromFile(iss), Catch::Matchers::ContainsSubstring("Failed to parse KnotsU value"));
}

TEST_CASE("Surface::LoadFromFile - Invalid Control Point Data", "[LoadFromFile]")
{
    std::string invalid_cp = R"(LIBNURBS TXT SURFACE
*DegreeU: 3
*DegreeV: 2
*ControlPoints:
2 2
1.0 2.0 invalid 1.0
4.0 5.0 6.0 1.0
)";
    std::istringstream iss(invalid_cp);
    Surface surface;
    REQUIRE_THROWS_WITH(surface.LoadFromFile(iss),
                        Catch::Matchers::ContainsSubstring("ControlPoints count does not match grid dimensions."));
}

TEST_CASE("Surface::LoadFromFile - Extra Spaces and Formatting", "[LoadFromFile]")
{
    std::string extra_spaces = R"(LIBNURBS     TXT    SURFACE

    *DegreeU:   3

    *DegreeV:   2

    *KnotsU:
    0,  0,0,   1,   2,3, 3,3

    *KnotsV:
    0, 0, 1, 1

    *ControlPoints:
    4    3
    1.0 2.0 3.0 1.0
    4.0    5.0 6.0 1.0
    7.0 8.0 9.0 1.0
    1.1 2.1 3.1 1.0
    4.1 5.1 6.1 1.0
    7.1 8.1 9.1 1.0
    1.2 2.2 3.2 1.0
    4.2 5.2 6.2 1.0
    7.2 8.2 9.2 1.0
    1.3 2.3 3.3 1.0
    4.3 5.3 6.3 1.0
    7.3 8.3 9.3 1.0
)";
    std::istringstream iss(extra_spaces);
    Surface surface;
    REQUIRE_NOTHROW(surface.LoadFromFile(iss));
    REQUIRE(surface.DegreeU == 3);
    REQUIRE(surface.DegreeV == 2);
    REQUIRE(surface.KnotsU.Values().size() == 8);
    REQUIRE(surface.KnotsV.Values().size() == 4);
    REQUIRE(surface.ControlPoints.UCount == 4);
    REQUIRE(surface.ControlPoints.VCount == 3);
    REQUIRE(surface.ControlPoints.Values.size() == 12);
}

TEST_CASE("Surface::LoadFromFile - Missing Key", "[LoadFromFile]")
{
    std::string missing_key = R"(LIBNURBS TXT SURFACE
3
)";
    std::istringstream iss(missing_key);
    Surface surface;
    REQUIRE_THROWS_WITH(surface.LoadFromFile(iss), Catch::Matchers::ContainsSubstring("Content found before any key"));
}

TEST_CASE("Surface::LoadFromFile - Empty File", "[LoadFromFile]")
{
    std::string empty_file = "";
    std::istringstream iss(empty_file);
    Surface surface;
    REQUIRE_THROWS_WITH(surface.LoadFromFile(iss), Catch::Matchers::ContainsSubstring("Failed to read header line"));
}

TEST_CASE("Surface::LoadFromFile and SaveToFile - Round-trip Test", "[SaveLoad]")
{
    Surface original_surface;
    // Initialize original_surface with test data
    original_surface.DegreeU = 3;
    original_surface.DegreeV = 2;
    original_surface.KnotsU.Values() = {0, 0, 0, 1, 2, 3, 3, 3};
    original_surface.KnotsV.Values() = {0, 0, 1, 1};
    original_surface.ControlPoints.UCount = 4;
    original_surface.ControlPoints.VCount = 3;
    original_surface.ControlPoints.Values = {
        Vec4{1.0, 2.0, 3.0, 1.0},
        Vec4{4.0, 5.0, 6.0, 1.0},
        Vec4{7.0, 8.0, 9.0, 1.0},
        Vec4{1.1, 2.1, 3.1, 1.0},
        Vec4{4.1, 5.1, 6.1, 1.0},
        Vec4{7.1, 8.1, 9.1, 1.0},
        Vec4{1.2, 2.2, 3.2, 1.0},
        Vec4{4.2, 5.2, 6.2, 1.0},
        Vec4{7.2, 8.2, 9.2, 1.0},
        Vec4{1.3, 2.3, 3.3, 1.0},
        Vec4{4.3, 5.3, 6.3, 1.0},
        Vec4{7.3, 8.3, 9.3, 1.0}
    };

    // Save to string stream
    std::ostringstream oss;
    REQUIRE_NOTHROW(original_surface.SaveToFile(oss));

    // Load from string stream
    std::istringstream iss(oss.str());
    Surface loaded_surface;
    REQUIRE_NOTHROW(loaded_surface.LoadFromFile(iss));

    // Compare original and loaded surfaces
    REQUIRE(loaded_surface.DegreeU == original_surface.DegreeU);
    REQUIRE(loaded_surface.DegreeV == original_surface.DegreeV);
    REQUIRE(loaded_surface.KnotsU.Values() == original_surface.KnotsU.Values());
    REQUIRE(loaded_surface.KnotsV.Values() == original_surface.KnotsV.Values());
    REQUIRE(loaded_surface.ControlPoints.UCount == original_surface.ControlPoints.UCount);
    REQUIRE(loaded_surface.ControlPoints.VCount == original_surface.ControlPoints.VCount);
    REQUIRE(loaded_surface.ControlPoints.Values.size() == original_surface.ControlPoints.Values.size());

    for (size_t i = 0; i < original_surface.ControlPoints.Values.size(); ++i)
    {
        const auto& original_cp = original_surface.ControlPoints.Values[i];
        const auto& loaded_cp = loaded_surface.ControlPoints.Values[i];
        REQUIRE(loaded_cp.x() == Approx(original_cp.x()));
        REQUIRE(loaded_cp.y() == Approx(original_cp.y()));
        REQUIRE(loaded_cp.z() == Approx(original_cp.z()));
        REQUIRE(loaded_cp.w() == Approx(original_cp.w()));
    }
}

TEST_CASE("Surface::SaveToFile and LoadFromFile - Binary Mode", "[SaveLoad]")
{
    Surface original_surface;
    // Initialize original_surface with test data (same as previous test)

    original_surface.DegreeU = 3;
    original_surface.DegreeV = 2;
    original_surface.KnotsU.Values() = {0, 0, 0, 1, 2, 3, 3, 3};
    original_surface.KnotsV.Values() = {0, 0, 1, 1};
    original_surface.ControlPoints.UCount = 4;
    original_surface.ControlPoints.VCount = 3;
    original_surface.ControlPoints.Values = {
        Vec4{1.0, 2.0, 3.0, 1.0},
        Vec4{4.0, 5.0, 6.0, 1.0},
        Vec4{7.0, 8.0, 9.0, 1.0},
        Vec4{1.1, 2.1, 3.1, 1.0},
        Vec4{4.1, 5.1, 6.1, 1.0},
        Vec4{7.1, 8.1, 9.1, 1.0},
        Vec4{1.2, 2.2, 3.2, 1.0},
        Vec4{4.2, 5.2, 6.2, 1.0},
        Vec4{7.2, 8.2, 9.2, 1.0},
        Vec4{1.3, 2.3, 3.3, 1.0},
        Vec4{4.3, 5.3, 6.3, 1.0},
        Vec4{7.3, 8.3, 9.3, 1.0}
    };

    // Save to binary string stream
    std::ostringstream oss(std::ios::binary);
    REQUIRE_NOTHROW(original_surface.SaveToFile(oss, true));

    // Load from binary string stream
    std::istringstream iss(oss.str(), std::ios::binary);
    Surface loaded_surface;
    REQUIRE_NOTHROW(loaded_surface.LoadFromFile(iss));

    // Compare original and loaded surfaces (same as previous test)
    REQUIRE(loaded_surface.DegreeU == original_surface.DegreeU);
    REQUIRE(loaded_surface.DegreeV == original_surface.DegreeV);
    REQUIRE(loaded_surface.KnotsU.Values() == original_surface.KnotsU.Values());
    REQUIRE(loaded_surface.KnotsV.Values() == original_surface.KnotsV.Values());
    REQUIRE(loaded_surface.ControlPoints.UCount == original_surface.ControlPoints.UCount);
    REQUIRE(loaded_surface.ControlPoints.VCount == original_surface.ControlPoints.VCount);
    REQUIRE(loaded_surface.ControlPoints.Values.size() == original_surface.ControlPoints.Values.size());

    for (size_t i = 0; i < original_surface.ControlPoints.Values.size(); ++i)
    {
        const auto& original_cp = original_surface.ControlPoints.Values[i];
        const auto& loaded_cp = loaded_surface.ControlPoints.Values[i];
        REQUIRE(loaded_cp.x() == Approx(original_cp.x()));
        REQUIRE(loaded_cp.y() == Approx(original_cp.y()));
        REQUIRE(loaded_cp.z() == Approx(original_cp.z()));
        REQUIRE(loaded_cp.w() == Approx(original_cp.w()));
    }
}

TEST_CASE("Surface::LoadFromFile - Large Surface Data", "[LoadFromFile]")
{
    // Generate a large surface
    Surface surface;
    surface.DegreeU = 4;
    surface.DegreeV = 4;
    int u_count = 50;
    int v_count = 50;
    surface.KnotsU.Values().resize(u_count + surface.DegreeU + 1, 0.0);
    surface.KnotsV.Values().resize(v_count + surface.DegreeV + 1, 0.0);

    for (int i = 0; i < u_count + surface.DegreeU + 1; ++i)
    {
        surface.KnotsU.Values()[i] = static_cast<Numeric>(i);
    }

    for (int i = 0; i < v_count + surface.DegreeV + 1; ++i)
    {
        surface.KnotsV.Values()[i] = static_cast<Numeric>(i);
    }

    surface.ControlPoints.UCount = u_count;
    surface.ControlPoints.VCount = v_count;
    surface.ControlPoints.Values.resize(u_count * v_count);

    for (int v = 0; v < v_count; ++v)
    {
        for (int u = 0; u < u_count; ++u)
        {
            surface.ControlPoints.Get(u, v) = Vec4{
                static_cast<Numeric>(u),
                static_cast<Numeric>(v),
                static_cast<Numeric>(u + v),
                1.0
            };
        }
    }

    // Save to string stream
    std::ostringstream oss;
    REQUIRE_NOTHROW(surface.SaveToFile(oss));

    // Load from string stream
    std::istringstream iss(oss.str());
    Surface loaded_surface;
    REQUIRE_NOTHROW(loaded_surface.LoadFromFile(iss));

    // Check dimensions
    REQUIRE(loaded_surface.ControlPoints.UCount == surface.ControlPoints.UCount);
    REQUIRE(loaded_surface.ControlPoints.VCount == surface.ControlPoints.VCount);

    // Spot-check some control points
    REQUIRE(loaded_surface.ControlPoints.Get(0, 0).x() == Approx(0.0));
    REQUIRE(loaded_surface.ControlPoints.Get(u_count - 1, v_count - 1).z() == Approx(2 * (u_count - 1)));
}

TEST_CASE("Surface::LoadFromFile - Unrecognized Key", "[LoadFromFile]")
{
    std::string unrecognized_key = R"(LIBNURBS TXT SURFACE
*DegreeU: 3
*DegreeV: 2
*UnknownKey:
SomeValue
)";
    std::istringstream iss(unrecognized_key);
    Surface surface;
    REQUIRE_THROWS_WITH(surface.LoadFromFile(iss), Catch::Matchers::ContainsSubstring("Unknown key encountered"));
}

TEST_CASE("Surface::LoadFromFile - Values Spanning Multiple Lines", "[LoadFromFile]")
{
    std::string multi_line_values = R"(LIBNURBS TXT SURFACE
*DegreeU:
3
*DegreeV:
2
*KnotsU:
0,0,0,
1,2,3,
3,3
*KnotsV:
0,0,
1,1
*ControlPoints:
4 3
1.0 2.0 3.0 1.0
4.0 5.0 6.0 1.0
7.0 8.0 9.0 1.0
1.1 2.1 3.1 1.0
4.1 5.1 6.1 1.0
7.1 8.1 9.1 1.0
1.2 2.2 3.2 1.0
4.2 5.2 6.2 1.0
7.2 8.2 9.2 1.0
1.3 2.3 3.3 1.0
4.3 5.3 6.3 1.0
7.3 8.3 9.3 1.0
)";
    std::istringstream iss(multi_line_values);
    Surface surface;
    REQUIRE_NOTHROW(surface.LoadFromFile(iss));
    REQUIRE(surface.KnotsU.Values().size() == 8);
    REQUIRE(surface.KnotsV.Values().size() == 4);
    REQUIRE(surface.ControlPoints.Values.size() == 12);
}



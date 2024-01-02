#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>

#include <libnurbs/Geometry/GeomRect.hpp>
#include <libnurbs/Surface/Surface.hpp>

using namespace libnurbs;


TEST_CASE("Geometry/GeomRect", "[libnurbs_Geometry]")
{
    SECTION("GetSurface 1")
    {
        GeomRect rect;
        rect.DegreeU = 1;
        rect.DegreeV = 1;
        rect.ControlPointCountU = 2;
        rect.ControlPointCountV = 2;
        rect.LeftBottomPoint = {0.0, 0.0, 0.0};
        rect.RightBottomPoint = {1.0, 0.0, 0.0};
        rect.LeftTopPoint = {0.0, 1.0, 0.0};
        rect.RightTopPoint = {1.0, 1.0, 0.0};
        Surface surface = rect.GetSurface();
        REQUIRE(surface.DegreeU == 1);
        REQUIRE(surface.DegreeV == 1);
        REQUIRE(surface.KnotsU.Count() == (rect.ControlPointCountU + rect.DegreeU + 1));
        REQUIRE(surface.KnotsV.Count() == (rect.ControlPointCountV + rect.DegreeV + 1));
        REQUIRE(surface.ControlPoints.Size() == (rect.ControlPointCountU * rect.ControlPointCountV));
        REQUIRE(surface.Evaluate(0.0, 0.0) == rect.LeftBottomPoint);
        REQUIRE(surface.Evaluate(1.0, 0.0) == rect.RightBottomPoint);
        REQUIRE(surface.Evaluate(0.0, 1.0) == rect.LeftTopPoint);
        REQUIRE(surface.Evaluate(1.0, 1.0) == rect.RightTopPoint);
    }
}

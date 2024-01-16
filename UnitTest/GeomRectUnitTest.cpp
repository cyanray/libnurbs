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
        REQUIRE(surface.DegreeU == rect.DegreeU);
        REQUIRE(surface.DegreeV == rect.DegreeV);
        REQUIRE(surface.KnotsU.Count() == (rect.ControlPointCountU + rect.DegreeU + 1));
        REQUIRE(surface.KnotsV.Count() == (rect.ControlPointCountV + rect.DegreeV + 1));
        REQUIRE(surface.ControlPoints.Count() == (rect.ControlPointCountU * rect.ControlPointCountV));
        REQUIRE(surface.Evaluate(0.0, 0.0) == rect.LeftBottomPoint);
        REQUIRE(surface.Evaluate(1.0, 0.0) == rect.RightBottomPoint);
        REQUIRE(surface.Evaluate(0.0, 1.0) == rect.LeftTopPoint);
        REQUIRE(surface.Evaluate(1.0, 1.0) == rect.RightTopPoint);
    }

    SECTION("GetSurface 2")
    {
        Vec3 N1{0.00, 0.00, 0.00};
        Vec3 N2{0.00, 1.00, 0.00};
        Vec3 N3{1.00, 0.00, 0.00};
        Vec3 N4{1.00, 1.00, 0.00};

        auto rect = GeomRect::Make(N1, N3, N2, N4);
        rect.DegreeU = 1;
        rect.DegreeV = 2;
        rect.ControlPointCountU = 3;
        rect.ControlPointCountV = 4;

        auto surface = rect.GetSurface();
        REQUIRE(surface.DegreeU == rect.DegreeU);
        REQUIRE(surface.DegreeV == rect.DegreeV);
        REQUIRE(surface.KnotsU.Count() == (rect.ControlPointCountU + rect.DegreeU + 1));
        REQUIRE(surface.KnotsV.Count() == (rect.ControlPointCountV + rect.DegreeV + 1));
        REQUIRE(surface.ControlPoints.Count() == (rect.ControlPointCountU * rect.ControlPointCountV));
        REQUIRE(surface.Evaluate(0.0, 0.0) == rect.LeftBottomPoint);
        REQUIRE(surface.Evaluate(1.0, 0.0) == rect.RightBottomPoint);
        REQUIRE(surface.Evaluate(0.0, 1.0) == rect.LeftTopPoint);
        REQUIRE(surface.Evaluate(1.0, 1.0) == rect.RightTopPoint);

        REQUIRE(surface.ControlPoints.Get(0, 0) == Vec4{0, 0, 0, 1});
        REQUIRE(surface.ControlPoints.Get(0, 1) == Vec4{0, 1.0/3, 0, 1});
        REQUIRE(surface.ControlPoints.Get(0, 2) == Vec4{0, 2.0/3, 0, 1});
        REQUIRE(surface.ControlPoints.Get(0, 3) == Vec4{0, 1.0, 0, 1});

        REQUIRE(surface.ControlPoints.Get(1, 0) == Vec4{0.5, 0, 0, 1});
        REQUIRE(surface.ControlPoints.Get(1, 1) == Vec4{0.5, 1.0/3, 0, 1});
        REQUIRE(surface.ControlPoints.Get(1, 2) == Vec4{0.5, 2.0/3, 0, 1});
        REQUIRE(surface.ControlPoints.Get(1, 3) == Vec4{0.5, 1.0, 0, 1});

        REQUIRE(surface.ControlPoints.Get(2, 0) == Vec4{1.0, 0, 0, 1});
        REQUIRE(surface.ControlPoints.Get(2, 1) == Vec4{1.0, 1.0/3, 0, 1});
        REQUIRE(surface.ControlPoints.Get(2, 2) == Vec4{1.0, 2.0/3, 0, 1});
        REQUIRE(surface.ControlPoints.Get(2, 3) == Vec4{1.0, 1.0, 0, 1});
    }
}

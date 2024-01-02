#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>

#include <libnurbs/Geometry/GeomSegment.hpp>
#include <libnurbs/Curve/Curve.hpp>

using namespace libnurbs;


Vec3 MiddlePoint(const Vec3& left, const Vec3& right)
{
    return (right - left) / 2;
}

TEST_CASE("Geometry/GeomSegment", "[libnurbs_Geometry]")
{
    SECTION("ToCurve 1")
    {
        GeomSegment seg;
        seg.Degree = 1;
        seg.ControlPointCount = 2;
        seg.LeftPoint = {0, 0, 0};
        seg.RightPoint = {10, 10, 10};
        Curve curve = seg.GetCurve();
        REQUIRE(curve.Degree == seg.Degree);
        REQUIRE(curve.ControlPoints.size() == seg.ControlPointCount);
        REQUIRE(curve.Knots.Count() == (seg.ControlPointCount + seg.Degree + 1));
        REQUIRE(curve.Evaluate(0.0) == seg.LeftPoint);
        REQUIRE(curve.Evaluate(1.0) == seg.RightPoint);
        // Vec3 mid_res = (curve.Evaluate(0.5) - MiddlePoint(seg.LeftPoint, seg.RightPoint));
        // REQUIRE(mid_res.x() == Catch::Approx(0.0));
        // REQUIRE(mid_res.y() == Catch::Approx(0.0));
        // REQUIRE(mid_res.z() == Catch::Approx(0.0));
    }
}

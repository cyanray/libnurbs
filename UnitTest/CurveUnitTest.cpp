#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>

#include <libnurbs/Curve/Curve.hpp>

#include <stdexcept>

using namespace Catch;
using namespace libnurbs;
using namespace std;


template<typename T>
std::vector<Numeric> ToStdVec(const T& mat)
{
    std::vector<Numeric> result(mat.size());
    std::copy(mat.data(), mat.data() + mat.size(), result.begin());
    return result;
}


void Check(const Curve& curve, Numeric x, const vector<Numeric>& result)
{
    Vec3 r = curve.Evaluate(x);
    REQUIRE_THAT(ToStdVec(r), Matchers::Approx(result));
}

TEST_CASE("Curve/Curve", "[libnurbs_BSplineCurve]")
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
        Check(curve, 0.0, {0.0, 0.0, 0.0});
        Check(curve, 1.0, {1.0, 1.0, 0.0});
        Check(curve, 0.5, {0.75, 0.25, 0.0});
        Check(curve, 0.3, {0.51, 0.09, 0.0});
        Check(curve, 0.678, {0.896316, 0.459684, 0.0});
    }

    SECTION("Evaluate (p=3)")
    {
        int degree = 3;
        KnotVector U{{0,0,0,0, 0.23, 0.67 ,1,1,1,1}};
        vector<Vec4> controlPoints{{0.0, 0.0, 0.0, 1.0},
                                   {1.0, 0.0, 0.0, 1.0},
                                   {1.0, 1.0, 0.0, 1.0},
                                   {2.0, 1.0, 0.0, 1.0},
                                   {2.0, -2.0, 0.0, 1.0},
                                   {0.0, -3.0, 0.0, 1.0}};
        Curve curve;
        curve.Degree = degree;
        curve.Knots = U;
        curve.ControlPoints = controlPoints;
        Check(curve, 0.0, {0.0, 0.0, 0.0});
        Check(curve, 1.0, {0.0, -3.0, 0.0});
        Check(curve, 0.123, {0.911390003336650, 0.224002988190835, 0.0});
        Check(curve, 0.345, {1.24695646619324, 0.808711157707757, 0.0});
        Check(curve, 0.5, {1.55856931399672, 0.748777149009726, 0.0});
        Check(curve, 0.789, {1.86924650801601, -0.950119518113616, 0.0});
        Check(curve, 0.987, {0.227165964922742, -2.88053915042935, 0.0});
    }

}
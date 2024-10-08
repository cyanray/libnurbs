#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>

#include <libnurbs/Curve/Curve.hpp>

#include <stdexcept>

using namespace Catch;
using namespace libnurbs;
using namespace std;


TEST_CASE("Curve/Evaluate (p=2)", "[curve][non_rational]")
{
    int degree = 2;
    KnotVector U{{0.0, 0.0, 0.0, 1.0, 1.0, 1.0}};
    vector<Vec4> controlPoints
    {
        {0.0, 0.0, 0.0, 1.0},
        {1.0, 0.0, 0.0, 1.0},
        {1.0, 1.0, 0.0, 1.0}
    };
    Curve curve;
    curve.Degree = degree;
    curve.Knots = U;
    curve.ControlPoints = controlPoints;

    SECTION("x = 0.0")
    {
        Vec3 value = curve.Evaluate(0.0);
        INFO("value: " << value.transpose());
        REQUIRE(value.x() == Approx(0.0));
        REQUIRE(value.y() == Approx(0.0));
        REQUIRE(value.z() == Approx(0.0));
    }

    SECTION("x = 1.0")
    {
        Vec3 value = curve.Evaluate(1.0);
        INFO("value: " << value.transpose());
        REQUIRE(value.x() == Approx(1.0));
        REQUIRE(value.y() == Approx(1.0));
        REQUIRE(value.z() == Approx(0.0));
    }

    SECTION("x = 0.5")
    {
        Vec3 value = curve.Evaluate(0.5);
        INFO("value: " << value.transpose());
        REQUIRE(value.x() == Approx(0.75));
        REQUIRE(value.y() == Approx(0.25));
        REQUIRE(value.z() == Approx(0.0));
    }

    SECTION("x = 0.3")
    {
        Vec3 value = curve.Evaluate(0.3);
        INFO("value: " << value.transpose());
        REQUIRE(value.x() == Approx(0.51));
        REQUIRE(value.y() == Approx(0.09));
        REQUIRE(value.z() == Approx(0.0));
    }

    SECTION("x = 0.678")
    {
        Vec3 value = curve.Evaluate(0.678);
        INFO("value: " << value.transpose());
        REQUIRE(value.x() == Approx(0.896316));
        REQUIRE(value.y() == Approx(0.459684));
        REQUIRE(value.z() == Approx(0.0));
    }
}


TEST_CASE("Curve/Evaluate (p=3)", "[curve][non_rational]")
{
    int degree = 3;
    KnotVector U{{0, 0, 0, 0, 0.23, 0.67, 1, 1, 1, 1}};
    vector<Vec4> controlPoints
    {
        {0.0, 0.0, 0.0, 1.0},
        {1.0, 0.0, 0.0, 1.0},
        {1.0, 1.0, 0.0, 1.0},
        {2.0, 1.0, 0.0, 1.0},
        {2.0, -2.0, 0.0, 1.0},
        {0.0, -3.0, 0.0, 1.0}
    };
    Curve curve;
    curve.Degree = degree;
    curve.Knots = U;
    curve.ControlPoints = controlPoints;
    SECTION("x = 0.0")
    {
        Vec3 value = curve.Evaluate(0.0);
        INFO("value: " << value.transpose());
        REQUIRE(value.x() == Approx(0.0));
        REQUIRE(value.y() == Approx(0.0));
        REQUIRE(value.z() == Approx(0.0));
    }

    SECTION("x = 1.0")
    {
        Vec3 value = curve.Evaluate(1.0);
        INFO("value: " << value.transpose());
        REQUIRE(value.x() == Approx(0.0));
        REQUIRE(value.y() == Approx(-3.0));
        REQUIRE(value.z() == Approx(0.0));
    }

    SECTION("x = 0.123")
    {
        Vec3 value = curve.Evaluate(0.123);
        INFO("value: " << value.transpose());
        REQUIRE(value.x() == Approx(0.911390003336650));
        REQUIRE(value.y() == Approx(0.224002988190835));
        REQUIRE(value.z() == Approx(0.0));
    }

    SECTION("x = 0.345")
    {
        Vec3 value = curve.Evaluate(0.345);
        INFO("value: " << value.transpose());
        REQUIRE(value.x() == Approx(1.24695646619324));
        REQUIRE(value.y() == Approx(0.808711157707757));
        REQUIRE(value.z() == Approx(0.0));
    }

    SECTION("x = 0.5")
    {
        Vec3 value = curve.Evaluate(0.5);
        INFO("value: " << value.transpose());
        REQUIRE(value.x() == Approx(1.55856931399672));
        REQUIRE(value.y() == Approx(0.748777149009726));
        REQUIRE(value.z() == Approx(0.0));
    }

    SECTION("x = 0.789")
    {
        Vec3 value = curve.Evaluate(0.789);
        INFO("value: " << value.transpose());
        REQUIRE(value.x() == Approx(1.86924650801601));
        REQUIRE(value.y() == Approx(-0.950119518113616));
        REQUIRE(value.z() == Approx(0.0));
    }

    SECTION("x = 0.987")
    {
        Vec3 value = curve.Evaluate(0.987);
        INFO("value: " << value.transpose());
        REQUIRE(value.x() == Approx(0.227165964922742));
        REQUIRE(value.y() == Approx(-2.88053915042935));
        REQUIRE(value.z() == Approx(0.0));
    }
}


TEST_CASE("Curve/Evaluate (circle arc)", "[curve][rational]")
{
    int degree = 2;
    KnotVector U{{0, 0, 0, 1, 1, 1}};
    vector<Vec4> controlPoints
    {
        {0.0, 1.0, 0.0, 2.0},
        {1.0, 1.0, 0.0, 1.0},
        {1.0, 0.0, 0.0, 1.0}
    };
    Curve curve;
    curve.Degree = degree;
    curve.Knots = U;
    curve.ControlPoints = controlPoints;
    SECTION("x = 0.0")
    {
        Vec3 value = curve.Evaluate(0.0);
        INFO("value: " << value.transpose());
        REQUIRE(value.x() == Approx(0.0));
        REQUIRE(value.y() == Approx(1.0));
        REQUIRE(value.z() == Approx(0.0));
    }

    SECTION("x = 1.0")
    {
        Vec3 value = curve.Evaluate(1.0);
        INFO("value: " << value.transpose());
        REQUIRE(value.x() == Approx(1.0));
        REQUIRE(value.y() == Approx(0.0));
        REQUIRE(value.z() == Approx(0.0));
    }

    SECTION("x = 0.123")
    {
        Vec3 value = curve.Evaluate(0.123);
        INFO("value: " << value.transpose());
        REQUIRE(value.x() == Approx(0.130499810923907));
        REQUIRE(value.y() == Approx(0.991448334180266));
        REQUIRE(value.z() == Approx(0.0));
    }

    SECTION("x = 0.345")
    {
        Vec3 value = curve.Evaluate(0.345);
        INFO("value: " << value.transpose());
        REQUIRE(value.x() == Approx(0.399555641083956));
        REQUIRE(value.y() == Approx(0.916708944909991));
        REQUIRE(value.z() == Approx(0.0));
    }

    SECTION("x = 0.5")
    {
        Vec3 value = curve.Evaluate(0.5);
        INFO("value: " << value.transpose());
        REQUIRE(value.x() == Approx(0.600000000000000));
        REQUIRE(value.y() == Approx(0.800000000000000));
        REQUIRE(value.z() == Approx(0.0));
    }

    SECTION("x = 0.789")
    {
        Vec3 value = curve.Evaluate(0.789);
        INFO("value: " << value.transpose());
        REQUIRE(value.x() == Approx(0.914753269680552));
        REQUIRE(value.y() == Approx(0.404012939902596));
        REQUIRE(value.z() == Approx(0.0));
    }

    SECTION("x = 0.987")
    {
        Vec3 value = curve.Evaluate(0.987);
        INFO("value: " << value.transpose());
        REQUIRE(value.x() == Approx(0.999662057112348));
        REQUIRE(value.y() == Approx(0.0259956067424605));
        REQUIRE(value.z() == Approx(0.0));
    }
}


TEST_CASE("Curve/EvaluateDerivative (non-rational)", "[curve][non_rational]")
{
    Curve curve;
    curve.Degree = 2;
    curve.Knots = KnotVector{{0, 0, 0, 1, 1, 1}};
    curve.ControlPoints = vector<Vec4>
    {
        {-1.0, 0.0, 0.0, 1.0},
        {1.0, 1.0, 0.0, 1.0},
        {1.0, 0.0, 0.0, 1.0}
    };

    SECTION("x = 0.5")
    {
        Vec3 tangent = curve.EvaluateDerivative(0.5, 1);
        REQUIRE(tangent.x() == Approx(2.0));
        REQUIRE(tangent.y() == Approx(0.0));
        REQUIRE(tangent.z() == Approx(0.0));
    }
}

TEST_CASE("Curve/EvaluateDerivative order=2 (non-rational)", "[curve][non_rational]")
{
    Curve curve;
    curve.Degree = 4;
    curve.Knots = KnotVector{{0.0, 0.0, 0.0, 0.0, 0.0, 0.333, 0.667, 1.0, 1.0, 1.0, 1.0, 1.0}};
    curve.ControlPoints =
    {
        {0.00, 0.00, 0.00, 1},
        {0.25, 0.00, 0.00, 1},
        {0.50, 0.00, 0.00, 1},
        {1.00, 0.00, 0.00, 1},
        {1.50, 0.00, 0.00, 1},
        {1.75, 0.00, 0.00, 1},
        {2.00, 0.00, 0.00, 1}
    };

    SECTION("x = 0.0")
    {
        Vec3 der = curve.EvaluateDerivative(0.0, 2);
        REQUIRE(der.x() == Approx(-13.547321));
        REQUIRE(der.y() == Approx(0.0));
        REQUIRE(der.z() == Approx(0.0));
    }

    SECTION("x = 0.3050847")
    {
        Vec3 der = curve.EvaluateDerivative(0.3050847, 2);
        REQUIRE(der.x() == Approx(1.197402021));
        REQUIRE(der.y() == Approx(0.0));
        REQUIRE(der.z() == Approx(0.0));
    }

    SECTION("x = 0.5084745763")
    {
        Vec3 der = curve.EvaluateDerivative(0.5084745763, 2);
        REQUIRE(der.x() == Approx(-0.057231905282));
        REQUIRE(der.y() == Approx(0.0));
        REQUIRE(der.z() == Approx(0.0));
    }
}

TEST_CASE("Curve/EvaluateDerivative order=2 (rational)", "[curve][rational]")
{
    const Numeric k = std::sqrt(2.0) / 2.0;
    Curve curve;
    curve.Degree = 2;
    curve.Knots = KnotVector{{0, 0, 0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1, 1, 1}};
    curve.ControlPoints = vector<Vec4>
    {
        {1, 0, 0, 1},
        {1, 1, 0, k},
        {0, 1, 0, 1},
        {-1, 1, 0, k},
        {-1, 0, 0, 1},
        {-1, -1, 0, k},
        {0, -1, 0, 1},
        {1, -1, 0, k},
        {1, 0, 0, 1},
    };

    SECTION("x = 0.0")
    {
        Vec3 der = curve.EvaluateDerivative(0.0, 1);
        REQUIRE(der.x() == Approx(0.0));
        REQUIRE(der.y() == Approx(5.6568542495));
        REQUIRE(der.z() == Approx(0.0));
    }

    SECTION("x = 0.5")
    {
        Vec3 der = curve.EvaluateDerivative(0.5, 1);
        REQUIRE(der.x() == Approx(0.0));
        REQUIRE(der.y() == Approx(-5.6568542495));
        REQUIRE(der.z() == Approx(0.0));
    }
}

TEST_CASE("Curve/SearchParameter", "[curve][non_rational]")
{
    int degree = 2;
    KnotVector U{{0.0, 0.0, 0.0, 1.0, 1.0, 1.0}};
    vector<Vec4> controlPoints
    {
        {0.0, 0.0, 0.0, 1.0},
        {1.0, 0.0, 0.0, 1.0},
        {1.0, 1.0, 0.0, 1.0}
    };
    Curve curve;
    curve.Degree = degree;
    curve.Knots = U;
    curve.ControlPoints = controlPoints;

    SECTION("x = 0.0")
    {
        Numeric u = curve.SearchParameter({0, 0, 0});;
        INFO("u: " << u);
        REQUIRE(u == Approx(0.0));
    }

    SECTION("x = 1.0")
    {
        Numeric u = curve.SearchParameter({1, 1, 0});;
        INFO("u: " << u);
        REQUIRE(u == Approx(1.0));
    }

    SECTION("x = 0.5")
    {
        Numeric u = curve.SearchParameter({0.75, 0.25, 0});;
        INFO("u: " << u);
        REQUIRE(u == Approx(0.5));
    }

    SECTION("x = 0.3")
    {
        Numeric u = curve.SearchParameter({0.51, 0.09, 0});;
        INFO("u: " << u);
        REQUIRE(u == Approx(0.3));
    }

    SECTION("x = 0.678")
    {
        Numeric u = curve.SearchParameter({0.896316, 0.459684, 0});;
        INFO("u: " << u);
        REQUIRE(u == Approx(0.678));
    }

    SECTION("x = 0.00123")
    {
        Numeric u = curve.SearchParameter(curve.Evaluate(0.00123));;
        INFO("u: " << u);
        REQUIRE(u == Approx(0.00123));
    }

    SECTION("x = 1e-15")
    {
        Numeric u = curve.SearchParameter(curve.Evaluate(1e-15), 0, 1e-15);
        INFO("u: " << u);
        REQUIRE(std::abs(u - 1e-15) < 1e-2);
    }
}

TEST_CASE("Curve/BinarySearchParameter", "[curve][non_rational]")
{
    int degree = 2;
    KnotVector U{{0.0, 0.0, 0.0, 1.0, 1.0, 1.0}};
    vector<Vec4> controlPoints
    {
        {0.0, 0.0, 0.0, 1.0},
        {1.0, 0.0, 0.0, 1.0},
        {1.0, 1.0, 0.0, 1.0}
    };
    Curve curve;
    curve.Degree = degree;
    curve.Knots = U;
    curve.ControlPoints = controlPoints;

    SECTION("x = 0.0")
    {
        Numeric u = curve.BinarySearchParameter({0, 0, 0});;
        INFO("u: " << u);
        REQUIRE(std::abs(u) < 1e-20);
    }

    SECTION("x = 1.0")
    {
        Numeric u = curve.BinarySearchParameter({1, 1, 0});;
        INFO("u: " << u);
        REQUIRE(u == Approx(1.0));
    }

    SECTION("x = 0.5")
    {
        Numeric u = curve.BinarySearchParameter({0.75, 0.25, 0});;
        INFO("u: " << u);
        REQUIRE(u == Approx(0.5));
    }

    SECTION("x = 0.3")
    {
        Numeric u = curve.BinarySearchParameter({0.51, 0.09, 0});;
        INFO("u: " << u);
        REQUIRE(u == Approx(0.3));
    }

    SECTION("x = 0.678")
    {
        Numeric u = curve.BinarySearchParameter({0.896316, 0.459684, 0});;
        INFO("u: " << u);
        REQUIRE(u == Approx(0.678));
    }

    SECTION("x = 0.00123")
    {
        Numeric u = curve.BinarySearchParameter(curve.Evaluate(0.00123));;
        INFO("u: " << u);
        REQUIRE(u == Approx(0.00123));
    }

    SECTION("x = 1e-15")
    {
        Numeric u = curve.BinarySearchParameter(curve.Evaluate(1e-15));
        INFO("u: " << u);
        REQUIRE(std::abs(u - 1e-15) < 1e-2);
    }
}

TEST_CASE("Curve/BinarySearchParameter Arc", "[curve][non_rational]")
{
    int degree = 2;
    KnotVector U{{0, 0, 0, 1, 1, 1}};
    vector<Vec4> controlPoints
    {
            {0.0, 1.0, 0.0, 2.0},
            {1.0, 1.0, 0.0, 1.0},
            {1.0, 0.0, 0.0, 1.0}
    };
    Curve curve;
    curve.Degree = degree;
    curve.Knots = U;
    curve.ControlPoints = controlPoints;

    SECTION("x = 0.0")
    {
        Numeric u = curve.BinarySearchParameter(curve.Evaluate(0.0));
        INFO("u: " << u);
        REQUIRE(std::abs(u) < 1e-20);
    }

    SECTION("x = 1.0")
    {
        Numeric u = curve.BinarySearchParameter(curve.Evaluate(1.0));;
        INFO("u: " << u);
        REQUIRE(u == Approx(1.0));
    }

    SECTION("x = 0.5")
    {
        Numeric u = curve.BinarySearchParameter(curve.Evaluate(0.5));;
        INFO("u: " << u);
        REQUIRE(u == Approx(0.5));
    }

    SECTION("x = 0.3")
    {
        Numeric u = curve.BinarySearchParameter(curve.Evaluate(0.3));;
        INFO("u: " << u);
        REQUIRE(u == Approx(0.3));
    }

    SECTION("x = 0.678")
    {
        Numeric u = curve.BinarySearchParameter(curve.Evaluate(0.678));;
        INFO("u: " << u);
        REQUIRE(u == Approx(0.678));
    }

    SECTION("x = 0.00123")
    {
        Numeric u = curve.BinarySearchParameter(curve.Evaluate(0.00123));;
        INFO("u: " << u);
        REQUIRE(u == Approx(0.00123));
    }

    SECTION("x = 1e-15")
    {
        Numeric u = curve.BinarySearchParameter(curve.Evaluate(1e-15));
        INFO("u: " << u);
        REQUIRE(std::abs(u - 1e-15) < 1e-2);
    }
}

TEST_CASE("Curve/InsertKnot", "[curve][non_rational]")
{
    int degree = 2;
    KnotVector U{{0.0, 0.0, 0.0, 1.0, 1.0, 1.0}};
    vector<Vec4> controlPoints
    {
        {0.0, 0.0, 0.0, 1.0},
        {1.0, 0.0, 0.0, 1.0},
        {1.0, 1.0, 0.0, 1.0}
    };
    Curve curve;
    curve.Degree = degree;
    curve.Knots = U;
    curve.ControlPoints = controlPoints;

    SECTION("1")
    {
        Curve new_curve = curve.InsertKnot(0.5);
        for (Numeric x: vector{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 1.0})
        {
            Vec3 value = curve.Evaluate(x);
            Vec3 new_value = new_curve.Evaluate(x);
            INFO("x: " << x);
            INFO("value: " << value.transpose());
            INFO("new_value: " << new_value.transpose());
            REQUIRE(value.x() == Approx(new_value.x()));
            REQUIRE(value.y() == Approx(new_value.y()));
            REQUIRE(value.z() == Approx(new_value.z()));
        }
    }

    SECTION("2")
    {
        for (Numeric knot: vector{0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9})
        {
            Curve new_curve = curve.InsertKnot(knot);
            for (Numeric x: vector{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 1.0})
            {
                Vec3 value = curve.Evaluate(x);
                Vec3 new_value = new_curve.Evaluate(x);
                INFO("knot: "<< knot);
                INFO("x: " << x);
                INFO("value: " << value.transpose());
                INFO("new_value: " << new_value.transpose());
                REQUIRE(value.x() == Approx(new_value.x()));
                REQUIRE(value.y() == Approx(new_value.y()));
                REQUIRE(value.z() == Approx(new_value.z()));
            }
        }
    }
}


TEST_CASE("Curve/InsertKnot 2", "[curve][rational]")
{
    const Numeric k = std::sqrt(2.0) / 2.0;
    Curve curve;
    curve.Degree = 2;
    curve.Knots = KnotVector{{0, 0, 0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1, 1, 1}};
    curve.ControlPoints = vector<Vec4>
    {
        {1, 0, 0, 1},
        {1, 1, 0, k},
        {0, 1, 0, 1},
        {-1, 1, 0, k},
        {-1, 0, 0, 1},
        {-1, -1, 0, k},
        {0, -1, 0, 1},
        {1, -1, 0, k},
        {1, 0, 0, 1},
    };

    SECTION("1")
    {
        for (Numeric knot: vector{0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9})
        {
            Curve new_curve = curve.InsertKnot(knot);
            for (Numeric x: vector{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 1.0})
            {
                Vec3 value = curve.Evaluate(x);
                Vec3 new_value = new_curve.Evaluate(x);
                INFO("knot: "<< knot);
                INFO("x: " << x);
                INFO("value: " << value.transpose());
                INFO("new_value: " << new_value.transpose());
                REQUIRE(value.x() == Approx(new_value.x()));
                REQUIRE(value.y() == Approx(new_value.y()));
                REQUIRE(value.z() == Approx(new_value.z()));
            }
        }
    }
}


TEST_CASE("Curve/InsertKnot 3", "[curve][rational]")
{
    const Numeric k = std::sqrt(2.0) / 2.0;
    Curve curve;
    curve.Degree = 2;
    curve.Knots = KnotVector{{0, 0, 0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1, 1, 1}};
    curve.ControlPoints = vector<Vec4>
    {
            {1, 0, 0, 1},
            {1, 1, 0, k},
            {0, 1, 0, 1},
            {-1, 1, 0, k},
            {-1, 0, 0, 1},
            {-1, -1, 0, k},
            {0, -1, 0, 1},
            {1, -1, 0, k},
            {1, 0, 0, 1},
        };

    SECTION("1")
    {
        for (Numeric knot2: vector{0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9})
        {
            for (Numeric knot1: vector{0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9})
            {
                Curve new_curve = curve.InsertKnot(knot1).InsertKnot(knot2);
                for (Numeric x: vector{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 1.0})
                {
                    Vec3 value = curve.Evaluate(x);
                    Vec3 new_value = new_curve.Evaluate(x);
                    INFO("knot 1: "<< knot1);
                    INFO("knot 2: "<< knot2);
                    INFO("x: " << x);
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


TEST_CASE("Curve/InsertKnot(list)", "[curve][rational]")
{
    const Numeric k = std::sqrt(2.0) / 2.0;
    Curve curve;
    curve.Degree = 2;
    curve.Knots = KnotVector{{0, 0, 0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1, 1, 1}};
    curve.ControlPoints = vector<Vec4>
    {
                {1, 0, 0, 1},
                {1, 1, 0, k},
                {0, 1, 0, 1},
                {-1, 1, 0, k},
                {-1, 0, 0, 1},
                {-1, -1, 0, k},
                {0, -1, 0, 1},
                {1, -1, 0, k},
                {1, 0, 0, 1},
            };

    SECTION("[0.2, 0.8]")
    {
        vector knots{0.2 , 0.8};
        auto new_curve = curve.InsertKnot(knots);
        for (Numeric x: vector{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 1.0})
        {
            Vec3 value = curve.Evaluate(x);
            Vec3 new_value = new_curve.Evaluate(x);
            INFO("x: " << x);
            INFO("value: " << value.transpose());
            INFO("new_value: " << new_value.transpose());
            REQUIRE(value.x() == Approx(new_value.x()));
            REQUIRE(value.y() == Approx(new_value.y()));
            REQUIRE(value.z() == Approx(new_value.z()));
        }
    }

    SECTION("[0.2, 0.2]")
    {
        vector knots{0.2 , 0.2};
        auto new_curve = curve.InsertKnot(knots);
        for (Numeric x: vector{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 1.0})
        {
            Vec3 value = curve.Evaluate(x);
            Vec3 new_value = new_curve.Evaluate(x);
            INFO("x: " << x);
            INFO("value: " << value.transpose());
            INFO("new_value: " << new_value.transpose());
            REQUIRE(value.x() == Approx(new_value.x()));
            REQUIRE(value.y() == Approx(new_value.y()));
            REQUIRE(value.z() == Approx(new_value.z()));
        }
    }

    SECTION("[0.2, 0.2, 0.2]")
    {
        vector knots{0.2 , 0.2, 0.2};
        auto new_curve = curve.InsertKnot(knots);
        for (Numeric x: vector{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 1.0})
        {
            Vec3 value = curve.Evaluate(x);
            Vec3 new_value = new_curve.Evaluate(x);
            INFO("x: " << x);
            INFO("value: " << value.transpose());
            INFO("new_value: " << new_value.transpose());
            REQUIRE(value.x() == Approx(new_value.x()));
            REQUIRE(value.y() == Approx(new_value.y()));
            REQUIRE(value.z() == Approx(new_value.z()));
        }
    }
}

TEST_CASE("Curve/InsertKnot times", "[curve][rational]")
{
    const Numeric k = std::sqrt(2.0) / 2.0;
    Curve curve;
    curve.Degree = 2;
    curve.Knots = KnotVector{{0, 0, 0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1, 1, 1}};
    curve.ControlPoints = vector<Vec4>
    {
                {1, 0, 0, 1},
                {1, 1, 0, k},
                {0, 1, 0, 1},
                {-1, 1, 0, k},
                {-1, 0, 0, 1},
                {-1, -1, 0, k},
                {0, -1, 0, 1},
                {1, -1, 0, k},
                {1, 0, 0, 1},
            };


    SECTION("[0.2] x3")
    {
        auto new_curve = curve.InsertKnot(0.2, 3);
        for (Numeric x: vector{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 1.0})
        {
            Vec3 value = curve.Evaluate(x);
            Vec3 new_value = new_curve.Evaluate(x);
            INFO("x: " << x);
            INFO("value: " << value.transpose());
            INFO("new_value: " << new_value.transpose());
            REQUIRE(value.x() == Approx(new_value.x()));
            REQUIRE(value.y() == Approx(new_value.y()));
            REQUIRE(value.z() == Approx(new_value.z()));
        }
    }
}


TEST_CASE("Curve/RemoveKnot", "[curve][rational]")
{
    const Numeric k = std::sqrt(2.0) / 2.0;
    Curve curve;
    curve.Degree = 2;
    curve.Knots = KnotVector{{0, 0, 0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1, 1, 1}};
    curve.ControlPoints = vector<Vec4>
    {
                    {1, 0, 0, 1},
                    {1, 1, 0, k},
                    {0, 1, 0, 1},
                    {-1, 1, 0, k},
                    {-1, 0, 0, 1},
                    {-1, -1, 0, k},
                    {0, -1, 0, 1},
                    {1, -1, 0, k},
                    {1, 0, 0, 1},
                };


    SECTION("1")
    {
        auto new_curve = curve.InsertKnot(0.2, 2);
        for (Numeric x: vector{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 1.0})
        {
            Vec3 value = curve.Evaluate(x);
            Vec3 new_value = new_curve.Evaluate(x);
            INFO("x: " << x);
            INFO("value: " << value.transpose());
            INFO("new_value: " << new_value.transpose());
            REQUIRE(value.x() == Approx(new_value.x()));
            REQUIRE(value.y() == Approx(new_value.y()));
            REQUIRE(value.z() == Approx(new_value.z()));
        }

        auto [removed_curve, t] = new_curve.RemoveKnot(0.2);
        REQUIRE(t == 1);
        for (Numeric x: vector{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 1.0})
        {
            Vec3 value = removed_curve.Evaluate(x);
            Vec3 new_value = new_curve.Evaluate(x);
            INFO("x: " << x);
            INFO("value: " << value.transpose());
            INFO("new_value: " << new_value.transpose());
            REQUIRE(value.x() == Approx(new_value.x()));
            REQUIRE(value.y() == Approx(new_value.y()));
            REQUIRE(value.z() == Approx(new_value.z()));
        }
    }
}


TEST_CASE("Curve/ElevateDegree", "[curve][rational]")
{
    const Numeric k = std::sqrt(2.0) / 2.0;
    Curve curve;
    curve.Degree = 2;
    curve.Knots = KnotVector{{0, 0, 0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1, 1, 1}};
    curve.ControlPoints = vector<Vec4>
    {
                        {1, 0, 0, 1},
                        {1, 1, 0, k},
                        {0, 1, 0, 1},
                        {-1, 1, 0, k},
                        {-1, 0, 0, 1},
                        {-1, -1, 0, k},
                        {0, -1, 0, 1},
                        {1, -1, 0, k},
                        {1, 0, 0, 1},
                    };


    SECTION("1")
    {
        auto new_curve = curve.ElevateDegree();
        for (Numeric x: vector{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 1.0})
        {
            Vec3 value = curve.Evaluate(x);
            Vec3 new_value = new_curve.Evaluate(x);
            INFO("x: " << x);
            INFO("value: " << value.transpose());
            INFO("new_value: " << new_value.transpose());
            REQUIRE(value.x() == Approx(new_value.x()));
            REQUIRE(value.y() == Approx(new_value.y()));
            REQUIRE(value.z() == Approx(new_value.z()));
        }


    }
}
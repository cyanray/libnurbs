#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <libnurbs/Core/KnotVector.hpp>

#include <stdexcept>

using namespace Catch;
using namespace libnurbs;
using namespace std;


TEST_CASE("KnotVector/Constructor", "[knot_vector, constructor]")
{
    SECTION("simple")
    {
        REQUIRE_NOTHROW(KnotVector{{0.0, 0.0, 0.0, 1.0, 1.0, 1.0}});
        REQUIRE_NOTHROW(KnotVector{{0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0}});
        REQUIRE_NOTHROW(KnotVector{{0.0, 1.0}});
        REQUIRE_NOTHROW(KnotVector{{0.0, 0.0, 1.0, 1.0}});

        REQUIRE_THROWS_AS((KnotVector{{0.0, 0.0, 1.0, 2.0}}), invalid_argument);
        REQUIRE_THROWS_AS((KnotVector{{0.0, 0.0, 2.0, 1.0}}), invalid_argument);
    }

    SECTION("static KnotVector::Uniform()")
    {
        KnotVector knots = KnotVector::Uniform(3, 10);
        REQUIRE(knots.Count() == 10);

        auto result = knots.GetKnotPairs();
        REQUIRE(result.size() == 4);
        REQUIRE(result[0].Value == 0.0);
        REQUIRE(result[0].Multiplicity == 4);

        REQUIRE(result[1].Value == Approx(0.33333333333));
        REQUIRE(result[1].Multiplicity == 1);

        REQUIRE(result[2].Value == Approx(0.66666666666));
        REQUIRE(result[2].Multiplicity == 1);

        REQUIRE(result[3].Value == 1.0);
        REQUIRE(result[3].Multiplicity == 4);
    }

}


TEST_CASE("KnotVector/GetKnotPairs", "[knot_vector, get_knot_pairs]")
{
    SECTION("simple")
    {
        KnotVector U{{0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0}};
        auto result = U.GetKnotPairs();
        REQUIRE(result.size() == 3);
        REQUIRE(result[0].Value == 0.0);
        REQUIRE(result[0].Multiplicity == 3);
        REQUIRE(result[1].Value == 0.5);
        REQUIRE(result[1].Multiplicity == 1);
        REQUIRE(result[2].Value == 1.0);
        REQUIRE(result[2].Multiplicity == 3);
    }
}


TEST_CASE("KnotVector/FindSpanIndex", "[knot_vector, find_span_index]")
{
    SECTION("simple")
    {
        KnotVector U{{0.0, 0.0, 0.0, 0.3, 0.5, 0.75, 1.0, 1.0, 1.0}};
        REQUIRE(U.FindSpanIndex(2, 0.0) == 2);
        REQUIRE(U.FindSpanIndex(2, 0.1) == 2);
        REQUIRE(U.FindSpanIndex(2, 0.2) == 2);
        REQUIRE(U.FindSpanIndex(2, 0.3) == 3);
        REQUIRE(U.FindSpanIndex(2, 0.4) == 3);
        REQUIRE(U.FindSpanIndex(2, 0.5) == 4);
        REQUIRE(U.FindSpanIndex(2, 0.6) == 4);
        REQUIRE(U.FindSpanIndex(2, 0.7) == 4);
        REQUIRE(U.FindSpanIndex(2, 0.8) == 5);
        REQUIRE(U.FindSpanIndex(2, 0.9) == 5);
        REQUIRE(U.FindSpanIndex(2, 1.0) == 5);
        REQUIRE(U.FindSpanIndex(2, 1.1) == 5); // u=1.1 is an illegal input, but acceptable
    }

    SECTION("simple2")
    {
        KnotVector U{{0.0, 0.0, 0.0, 0.3, 0.5, 0.5, 0.75, 1.0, 1.0, 1.0}};
        REQUIRE(U.FindSpanIndex(2, 0.0) == 2);
        REQUIRE(U.FindSpanIndex(2, 0.1) == 2);
        REQUIRE(U.FindSpanIndex(2, 0.2) == 2);
        REQUIRE(U.FindSpanIndex(2, 0.3) == 3);
        REQUIRE(U.FindSpanIndex(2, 0.4) == 3);
        REQUIRE(U.FindSpanIndex(2, 0.5) == 5);
        REQUIRE(U.FindSpanIndex(2, 0.6) == 5);
        REQUIRE(U.FindSpanIndex(2, 0.7) == 5);
        REQUIRE(U.FindSpanIndex(2, 0.8) == 6);
        REQUIRE(U.FindSpanIndex(2, 0.9) == 6);
        REQUIRE(U.FindSpanIndex(2, 1.0) == 6);
        REQUIRE(U.FindSpanIndex(2, 1.1) == 6); // u=1.1 is an illegal input, but acceptable
    }
}

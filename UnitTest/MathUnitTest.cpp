#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

using namespace Catch;

TEST_CASE("Math funtions", "[libnurbs_math]")
{
    SECTION("GetTriangleArea")
    {
        REQUIRE(0.5 == Approx(0.5));
    }
}
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <libnurbs/Algorithm/MathUtils.hpp>

using namespace std;
using namespace libnurbs;
using namespace Catch;

TEST_CASE("Algorithm/MathUtils", "[libnurbs_MathUtils]")
{
    SECTION("Factorial")
    {
        REQUIRE(Factorial(0) == 1);
        REQUIRE(Factorial(1) == 1);
        REQUIRE(Factorial(2) == 2);
        REQUIRE(Factorial(3) == 6);
        REQUIRE(Factorial(4) == 24);
        REQUIRE(Factorial(5) == 120);
        REQUIRE(Factorial(6) == 720);
        REQUIRE(Factorial(7) == 5040);
        REQUIRE(Factorial(8) == 40320);
        REQUIRE(Factorial(9) == 362880);
    }
}
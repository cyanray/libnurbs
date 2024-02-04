#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <libnurbs/Core/Grid.hpp>

using namespace std;
using namespace libnurbs;
using namespace Catch;

TEST_CASE("Core/Grid", "[grid]")
{
    Grid<double> grid(3, 3);
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            grid.Get(i, j) = i * 3.0 + j;
        }
    }

    SECTION("InsertV")
    {
        grid.InsertV(2, -1.0);
        // check new values after V insertion
        REQUIRE(grid.Get(0, 2) == Approx(-1.0));
        REQUIRE(grid.Get(1, 2) == Approx(-1.0));
        REQUIRE(grid.Get(2, 2) == Approx(-1.0));
        // check original values
        REQUIRE(grid.Get(0, 0) == Approx(0.0));
        REQUIRE(grid.Get(1, 1) == Approx(4.0));
    }

    SECTION("InsertU")
    {
        grid.InsertU(2, -2.0);
        // check new values after U insertion
        REQUIRE(grid.Get(2, 0) == Approx(-2.0));
        REQUIRE(grid.Get(2, 1) == Approx(-2.0));
        REQUIRE(grid.Get(2, 2) == Approx(-2.0));
        // check original values
        REQUIRE(grid.Get(0, 0) == Approx(0.0));
        REQUIRE(grid.Get(1, 1) == Approx(4.0));
    }

    SECTION("InsertRow and InsertColumn")
    {
        grid.InsertV(2, -1.0);
        // check new values after V insertion
        REQUIRE(grid.Get(0, 2) == Approx(-1.0));
        REQUIRE(grid.Get(1, 2) == Approx(-1.0));
        REQUIRE(grid.Get(2, 2) == Approx(-1.0));
        // check original values
        REQUIRE(grid.Get(0, 0) == Approx(0.0));
        REQUIRE(grid.Get(1, 1) == Approx(4.0));


        grid.InsertU(2, -2.0);
        // check new values after U insertion
        REQUIRE(grid.Get(2, 0) == Approx(-2.0));
        REQUIRE(grid.Get(2, 1) == Approx(-2.0));
        REQUIRE(grid.Get(2, 2) == Approx(-2.0));
        // check original values
        REQUIRE(grid.Get(0, 0) == Approx(0.0));
        REQUIRE(grid.Get(1, 1) == Approx(4.0));
    }
}

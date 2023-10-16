#pragma once

#include <cassert>

namespace libnurbs
{

    inline constexpr int Factorial(int n)
    {
        constexpr const int result_array[10]{1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880};
        assert(n >= 0 && n <= 9);
        return result_array[n];
    }

}
#pragma once

#include <libnurbs/Core/Typedefs.hpp>
#include <libnurbs/Core/KnotVector.hpp>

namespace libnurbs
{
    struct BSplineBasis
    {
        static VecX Evaluate(int degree, const KnotVector& knot_vec, Numeric x);
    };
}
#pragma once

#include <libnurbs/Core/Typedefs.hpp>
#include <libnurbs/Core/KnotVector.hpp>
#include <vector>


namespace libnurbs
{
    struct BSplineBasis
    {
        static VecX Evaluate(int degree, const vector<Numeric >& knots, int index_span, Numeric x);
    };
}
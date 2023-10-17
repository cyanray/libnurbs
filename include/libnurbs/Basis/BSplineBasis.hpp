#pragma once

#include <libnurbs/Core/Typedefs.hpp>
#include <vector>

using std::vector;

namespace libnurbs
{
    class KnotVector;

    struct BSplineBasis
    {
        static VecX Evaluate(int degree, const KnotVector& knot_vec, Numeric x);

        static VecX Evaluate(int degree, const vector<Numeric>& knots, int index_span, Numeric x);

        static MatX EvaluateDerivative(int degree, const KnotVector& knot_vec, Numeric x, int order = 1);

        static MatX EvaluateDerivative(int degree, const std::vector<Numeric>& knots, int index_span, Numeric x, int order);

    };
}
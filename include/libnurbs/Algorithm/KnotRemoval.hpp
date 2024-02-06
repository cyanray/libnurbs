#pragma once
#include "libnurbs/Core/Typedefs.hpp"
#include <vector>

namespace libnurbs
{
    class KnotVector;

    int KnotRemoval(KnotVector& knot_vector,
                    std::vector<Vec4>& points,
                    int degree,
                    Numeric knot_remove,
                    int times,
                    Numeric tolerance);
}

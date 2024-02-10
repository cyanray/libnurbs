#pragma once
#include <vector>

#include "libnurbs/Core/Typedefs.hpp"

namespace libnurbs
{
    class KnotVector;

    void DegreeElevation(KnotVector& knot_vector,
                      std::vector<Vec4>& points,
                      int& degree,
                      int times);

    int KnotRemoval(KnotVector& knot_vector,
                std::vector<Vec4>& points,
                int degree,
                Numeric knot_remove,
                int times,
                Numeric tolerance);
}

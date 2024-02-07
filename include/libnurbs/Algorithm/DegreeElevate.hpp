#pragma once
#include <vector>

#include "libnurbs/Core/Typedefs.hpp"

namespace libnurbs
{
    class KnotVector;

    void DegreeElevate(KnotVector& knot_vector,
                      std::vector<Vec4>& points,
                      int degree,
                      int times,
                      Numeric tolerance);

}

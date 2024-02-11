#pragma once
#include <vector>

#include "libnurbs/Core/Typedefs.hpp"

namespace libnurbs
{
    class KnotVector;

    void NurbsDegreeElevation(KnotVector& knot_vector,
                      std::vector<Vec4>& points,
                      int& degree,
                      int times);

    bool BezierDegreeReduction(std::vector<Vec4>& bpts, int degree);

    bool NurbsDegreeReduction(KnotVector& knot_vector,
                std::vector<Vec4>& points,
                Numeric tolerance);
}

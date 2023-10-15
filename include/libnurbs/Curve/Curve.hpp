#pragma once

#include <libnurbs/Core/Typedefs.hpp>
#include <libnurbs/Core/KnotVector.hpp>
#include <vector>

using std::vector;

namespace libnurbs
{
    class Curve
    {
    public:
        int Degree{INVALID_DEGREE};
        KnotVector Knots{};
        vector<Vec4> ControlPoints{};
    public:
        Curve() = default;

        bool Validate() const;

        [[nodiscard]] Vec3 Evaluate(Numeric x) const;
    };
}
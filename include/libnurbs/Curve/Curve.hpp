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

        [[nodiscard]] Vec3 EvaluateDerivative(Numeric x, int order) const;

        [[nodiscard]] vector<Vec3> EvaluateAll(Numeric x, int order) const;

        [[nodiscard]] Numeric SearchParameter(const Vec3& point) const;

    private:
        [[nodiscard]] vector<Vec4> HomogeneousDerivative(Numeric x, int order) const;
    };
}
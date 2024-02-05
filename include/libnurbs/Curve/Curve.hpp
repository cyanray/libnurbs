#pragma once

#include <span>
#include <vector>
#include <libnurbs/Core/KnotVector.hpp>
#include <libnurbs/Core/Typedefs.hpp>

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

        [[nodiscard]] Numeric SearchParameter(const Vec3& point,
                                              Numeric init = 0.5,
                                              Numeric epsion = 1e-8,
                                              Numeric max_iteration_count = 512) const;

        [[nodiscard]] Curve InsertKnot(Numeric knot_value) const;

        [[nodiscard]] Curve InsertKnot(std::span<Numeric> knots_to_insert) const;

    private:
        [[nodiscard]] vector<Vec4> HomogeneousDerivative(Numeric x, int order) const;
    };
}

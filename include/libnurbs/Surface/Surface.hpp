#pragma once

#include <libnurbs/Core/Typedefs.hpp>
#include <libnurbs/Core/KnotVector.hpp>
#include <libnurbs/Core/Grid.hpp>

namespace libnurbs
{
    class Surface
    {
    public:
        int DegreeU{INVALID_DEGREE};
        int DegreeV{INVALID_DEGREE};
        KnotVector KnotsU{};
        KnotVector KnotsV{};
        ControlPointGrid ControlPoints{};
    public:
        Surface() = default;

        bool Validate() const;

        [[nodiscard]] Vec3 Evaluate(Numeric u, Numeric v) const;

        [[nodiscard]] Vec3 EvaluateDerivative(Numeric u, Numeric v, int order_u, int order_v) const;

        [[nodiscard]] Grid<Vec3> EvaluateAll(Numeric u, Numeric v, int order_u, int order_v) const;

    private:
        [[nodiscard]] Grid<Vec4> HomogeneousDerivative(Numeric u, Numeric v, int order_u, int order_v) const;

    };


}
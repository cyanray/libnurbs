#pragma once

#include <span>
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

        [[nodiscard]] Vec3 Evaluate(Numeric u, Numeric v) const;

        [[nodiscard]] Vec3 EvaluateDerivative(Numeric u, Numeric v, int order_u, int order_v) const;

        [[nodiscard]] Grid<Vec3> EvaluateAll(Numeric u, Numeric v, int order_u, int order_v) const;

        [[nodiscard]] std::pair<Numeric, Numeric> SearchParameter(const Vec3& point,
                                                                  Numeric init_u = 0.5,
                                                                  Numeric init_v = 0.5,
                                                                  Numeric epsion = 1e-8,
                                                                  Numeric max_iteration_count = 512) const;

        [[nodiscard]] Surface InsertKnotU(Numeric knot_value) const;

        [[nodiscard]] Surface InsertKnotU(Numeric knot_value, int times) const;

        [[nodiscard]] Surface InsertKnotU(std::span<Numeric> knots_to_insert) const;

        [[nodiscard]] Surface InsertKnotV(Numeric knot_value) const;

        [[nodiscard]] Surface InsertKnotV(Numeric knot_value, int times) const;

        [[nodiscard]] Surface InsertKnotV(std::span<Numeric> knots_to_insert) const;

        [[nodiscard]] std::tuple<Surface, int> RemoveKnotU(Numeric knot_remove,
                                                           int times = 1,
                                                           Numeric tolerance = 1e-12) const;

        [[nodiscard]] std::tuple<Surface, int> RemoveKnotV(Numeric knot_remove,
                                                   int times = 1,
                                                   Numeric tolerance = 1e-12) const;

        [[nodiscard]] Surface ElevateDegreeU(int times = 1) const;

        [[nodiscard]] Surface ElevateDegreeV(int times = 1) const;

    private:
        [[nodiscard]] Grid<Vec4> HomogeneousDerivative(Numeric u, Numeric v, int order_u, int order_v) const;
    };
}

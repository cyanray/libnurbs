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

        Surface& LoadFromFile(const string& filename);

        Surface& LoadFromFile(std::istream& is);

        void SaveToFile(const string& filename, bool binary_mode = false) const;

        void SaveToFile(std::ostream& os, bool binary_mode = false) const;

        [[nodiscard]] Vec3 Evaluate(Numeric u, Numeric v) const;

        [[nodiscard]] Vec3 EvaluateDerivative(Numeric u, Numeric v, int order_u, int order_v) const;

        [[nodiscard]] Grid<Vec3> EvaluateAll(Numeric u, Numeric v, int order_u, int order_v) const;

        [[nodiscard]] bool IsRational() const;

        /**
         * @brief Searches for the parameter value
         *        corresponding to a point on the surface using the BFGS method.
         * @param point Point on surface
         * @param init_u Initial value for parameter u, default is 0.5
         * @param init_v Initial value for parameter v, default is 0.5
         * @param epsilon Precision, default is 1e-8
         * @param max_iteration_count Maximum number of iterations allowed, default is 512
         * @return pair(u, v)
         */
        [[nodiscard]] auto SearchParameter(const Vec3& point,
                                           Numeric init_u = 0.5,
                                           Numeric init_v = 0.5,
                                           Numeric epsilon = 1e-8,
                                           Numeric max_iteration_count = 512) const -> std::pair<Numeric, Numeric>;

        [[nodiscard]] auto SearchParameterOn(const Vec3& point, int direction, Numeric constant,
                                             Numeric init_value = 0.5,
                                             Numeric epsilon = 1e-8,
                                             Numeric max_iteration_count = 512) const -> std::pair<Numeric, Numeric>;


        [[nodiscard]] auto BinarySearchParameterOn(const Vec3& point,
                                                   int direction,
                                                   Numeric constant,
                                                   Numeric epsilon = 1e-8,
                                                   Numeric max_iteration_count = 512) const -> std::pair<Numeric, Numeric>;


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

        [[nodiscard]] Surface Transform(const Mat3x3& R) const;


        enum class AlignAxis
        {
            XAxis,
            YAxis,
            ZAxis
        };

        [[nodiscard]] Surface AlignParameterDomain(AlignAxis u_axis, AlignAxis v_axis);

    private:
        [[nodiscard]] Grid<Vec4> HomogeneousDerivative(Numeric u, Numeric v, int order_u, int order_v) const;
    };
}

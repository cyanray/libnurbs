#pragma once

#include <span>
#include <vector>
#include <tuple>
#include <string>
#include <libnurbs/Core/KnotVector.hpp>
#include <libnurbs/Core/Typedefs.hpp>

using std::vector;
using std::string;

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

        Curve& LoadFromFile(const string& filename);

        Curve& LoadFromFile(std::istream& is);

        void SaveToFile(const string& filename, bool binary_mode = false) const;

        void SaveToFile(std::ostream& os, bool binary_mode = false) const;


        [[nodiscard]] Vec3 Evaluate(Numeric x) const;

        [[nodiscard]] Vec3 EvaluateDerivative(Numeric x, int order) const;

        [[nodiscard]] vector<Vec3> EvaluateAll(Numeric x, int order) const;

        [[nodiscard]] bool IsRational() const;

        [[nodiscard]] Numeric SearchParameter(const Vec3& point,
                                              Numeric init = 0.5,
                                              Numeric epsilon = 1e-8,
                                              Numeric max_iteration_count = 512) const;

        [[nodiscard]] Numeric BinarySearchParameter(const Vec3& point,
                                                    Numeric epsilon = 1e-8,
                                                    int max_iteration_count = 512) const;


        [[nodiscard]] Curve InsertKnot(Numeric knot_value) const;

        [[nodiscard]] Curve InsertKnot(Numeric knot_value, int times) const;

        [[nodiscard]] Curve InsertKnot(std::span<Numeric> knots_to_insert) const;

        [[nodiscard]] std::tuple<Curve, int> RemoveKnot(Numeric knot_remove,
                                                        int times = 1,
                                                        Numeric tolerance = 1e-12) const;

        [[nodiscard]] Curve ElevateDegree(int times = 1) const;

        [[nodiscard]] Curve Transform(const Mat3x3& R) const;

    private:
        [[nodiscard]] vector<Vec4> HomogeneousDerivative(Numeric x, int order) const;
    };
}

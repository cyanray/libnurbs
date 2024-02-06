#include <libnurbs/Curve/Curve.hpp>
#include <libnurbs/Basis/BSplineBasis.hpp>
#include <cassert>

#include "libnurbs/Algorithm/MathUtils.hpp"

namespace libnurbs
{
    Vec3 Curve::Evaluate(Numeric x) const
    {
        assert(x >= 0 && x <= 1);
        int index_span = Knots.FindSpanIndex(Degree, x);
        VecX basis = BSplineBasis::Evaluate(Degree, Knots.Values(), index_span, x);
        assert(basis.size() == Degree + 1);
        Vec4 result = Vec4::Zero();
        for (int i = 0; i <= Degree; i++)
        {
            auto point = ToHomo(ControlPoints[index_span - Degree + i]);
            result += basis(i) * point;
        }
        return result.head<3>() / result.w();
    }


    std::vector<Vec4> Curve::HomogeneousDerivative(Numeric x, int order) const
    {
        assert(x >= 0 && x <= 1);
        int index_span = Knots.FindSpanIndex(Degree, x);
        MatX basis = BSplineBasis::EvaluateAll(Degree, Knots.Values(), index_span, x, order);
        assert(basis.rows() == order + 1);
        assert(basis.cols() == Degree + 1);
        std::vector<Vec4> result(order + 1);
        for (int k = 0; k <= order; ++k)
        {
            Vec4 tmp = Vec4::Zero();
            for (int i = 0; i <= Degree; i++)
            {
                auto point = ToHomo(ControlPoints[index_span - Degree + i]);
                tmp.noalias() += point * basis(k, i);
            }
            result[k] = tmp;
        }
        return result;
    }

    Vec3 Curve::EvaluateDerivative(Numeric x, int order) const
    {
        return EvaluateAll(x, order)[order];
    }

    vector<Vec3> Curve::EvaluateAll(Numeric x, int order) const
    {
        auto homo_ders = HomogeneousDerivative(x, order);
        vector<Vec3> result(order + 1, Vec3::Zero());

        // Compute rational derivatives
        Numeric Wders0 = homo_ders[0].w();
        for (int k = 0; k <= order; k++)
        {
            Vec3 Aders = homo_ders[k].head<3>();
            for (int i = 1; i <= k; i++)
            {
                Numeric Wders = homo_ders[i].w();
                Aders.noalias() -= Binomial(k, i) * Wders * result[k - i];
            }
            result[k] = (Aders / Wders0);
        }
        return result;
    }

    Numeric Curve::SearchParameter(const Vec3& point, Numeric init, Numeric epsion, Numeric max_iteration_count) const
    {
        auto Ri = [&point, this](Numeric u) -> Vec3 { return Evaluate(u) - point; };
        auto fi = [this, &Ri](Numeric u) -> Numeric
        {
            auto ri = Ri(u);
            auto Cu = EvaluateDerivative(u, 1);
            return ri.dot(Cu);
        };
        auto Ji = [this, &Ri](Numeric u) -> Numeric
        {
            auto Cu = EvaluateDerivative(u, 1);
            auto Cuu = EvaluateDerivative(u, 2);
            auto ri = Ri(u);
            return Cu.dot(Cu) + ri.dot(Cuu);
        };

        Numeric u_last = init;
        Numeric res = 1;
        int count = 0;
        // coefficient *s* is used to prevent uv_last from going out of range
        Numeric s = 1.0;
        while (res >= epsion && (count++ < max_iteration_count))
        {
            auto f = fi(u_last);
            auto j = Ji(u_last);
            Numeric u_new = u_last - f / j * s;

            res = std::abs(u_new - u_last);
            u_last = u_new;

            if (u_last < 0 || u_last > 1)
            {
                s /= 2;
                u_last = u_last > 1.0 ? 1.0 : u_last;
                u_last = u_last < 0.0 ? 0.0 : u_last;
            }
        }
        return u_last;
    }

    Curve Curve::InsertKnot(Numeric knot_value) const
    {
        Curve result{*this};
        int k = result.Knots.InsertKnot(knot_value);
        result.ControlPoints.insert(result.ControlPoints.begin() + k, Vec4::Zero());
        const auto& knots = this->Knots.Values();
        for (int i = k - Degree + 1; i <= k; i++)
        {
            Numeric alpha = (knot_value - knots[i]) / (knots[i + Degree] - knots[i]);
            auto point_i = ToHomo(this->ControlPoints[i]);
            auto point_im1 = ToHomo(this->ControlPoints[i - 1]);
            result.ControlPoints[i] = FromHomo(alpha * point_i + (1 - alpha) * point_im1);
        }
        return result;
    }

    Curve Curve::InsertKnot(Numeric knot_value, int times) const
    {
        // TODO: better implementation
        std::vector list(times, knot_value);
        return InsertKnot(list);
    }

    Curve Curve::InsertKnot(std::span<Numeric> knots_to_insert) const
    {
        Curve result = *this;
        auto& knots = result.Knots.Values();
        knots.reserve(knots.size() + knots_to_insert.size());
        auto& points = result.ControlPoints;
        points.reserve(points.size() + knots_to_insert.size());
        for (auto knot_value: knots_to_insert)
        {
            int k = result.Knots.FindSpanIndex(Degree, knot_value);
            points.insert(points.begin() + k, points[k]);
            Vec4 point_last = ToHomo(points[k - Degree]);
            for (int i = k - Degree + 1; i <= k; i++)
            {
                Numeric alpha = (knot_value - knots[i]) / (knots[i + Degree] - knots[i]);
                auto point_i = ToHomo(points[i]);
                Vec4 point_new = FromHomo(alpha * point_i + (1 - alpha) * point_last);
                point_last = point_i;
                points[i] = point_new;
            }
            knots.insert(knots.begin() + k + 1, knot_value);
        }
        return result;
    }

    std::tuple<Curve, int> Curve::RemoveKnot(Numeric knot_remove, int times, Numeric tolerance) const
    {
        constexpr static auto CalcTOL = [](const vector<Vec4>& points, Numeric epsilon)-> Numeric
        {
            Numeric w_min = 1e12;
            Numeric dist_max = 0.0;
            for (const auto& point: points)
            {
                w_min = std::min(w_min, point.w());
                dist_max = std::max(dist_max, point.head<3>().norm());
            }
            return (epsilon * w_min) / (1 + dist_max);
        };

        assert(knot_remove > 0 && knot_remove < 1);
        assert(times > 0);

        Curve result{*this};
        int degree = result.Degree;
        auto& points = result.ControlPoints;
        const auto& points_original = this->ControlPoints;
        auto& knots = result.Knots.Values();
        const auto& knots_original = this->Knots.Values();

        int s = result.Knots.GetMultiplicity(knot_remove);
        int r = result.Knots.FindSpanIndex(degree, knot_remove);

        double tol = CalcTOL(result.ControlPoints, tolerance);

        int n = (int)result.ControlPoints.size() - 1;
        int m = n + degree + 1;
        int order = degree + 1;
        int last = r - s;
        int first = r - degree;

        std::vector<Vec4> temp(2 * degree + 1);
        int t;
        for (t = 0; t < times; t++)
        {
            int off = first - 1;
            temp[0] = points[off];
            temp[last + 1 - off] = points[last + 1];
            int i = first;
            int j = last;
            int ii = 1;
            int jj = last - off;

            while (j - i >= t)
            {
                double alphai = (knot_remove - knots[i]) / (knots[i + order + t] - knots[i]);
                double alphaj = (knot_remove - knots[j - t]) / (knots[j + order] - knots[j - t]);

                temp[ii] = (points[i] - (1.0 - alphai) * temp[ii - 1]) / alphai;
                temp[jj] = (points[j] - alphaj * temp[jj + 1]) / (1.0 - alphaj);

                i = i + 1;
                ii = ii + 1;

                j = j - 1;
                jj = jj - 1;
            }

            Numeric dist;
            if (j - i < t)
            {
                dist = (ToHomo(temp[ii - 1]) - ToHomo(temp[jj + 1])).norm();
            }
            else
            {
                double alphai = (knot_remove - knots[i]) / (knots[i + order + t] - knots[i]);
                dist = (alphai * ToHomo(temp[ii + t + 1]) + (1.0 - alphai) * ToHomo(temp[ii - 1])).norm();
            }
            if (dist > tol) break;

            i = first;
            j = last;
            while (j - i > t)
            {
                points[i] = temp[i - off];
                points[j] = temp[j - off];
                i = i + 1;
                j = j - 1;
            }

            first = first - 1;
            last = last + 1;
        }

        // no move operation success, return the original curve
        // TODO: return operation status
        if (t == 0) return {result, t};

        int j = (2 * r - s - degree) / 2, i = j;
        for (int k = 1; k < t; k++)
        {
            i += (k % 2 == 1) ? 1 : 0;
            j -= (k % 2 == 1) ? 0 : 1;
        }

        /* Update control points */
        for (int k = i + 1; k <= n; k++, j++)
        {
            points[j] = points_original[k];
        }
        points.resize(points.size() - times);

        /* Update knot vector */
        for (int k = r + 1; k <= m; k++)
        {
            knots[k - times] = knots_original[k];
        }
        knots.resize(knots.size() - times);

        return {result, t};
    }
}

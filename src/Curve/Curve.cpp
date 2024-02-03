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
            auto point = ControlPoints[index_span - Degree + i];
            point.head<3>() *= point(3);
            result += basis(i) * point;
        }
        return result.head<3>() / result(3);
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
                auto point = ControlPoints[index_span - Degree + i];
                point.head<3>() *= point.w();
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
            Numeric u_new = u_last -f / j * s;

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
        result.ControlPoints.insert(result.ControlPoints.begin() + k - 1 ,Vec4::Zero());
        const auto& knots = this->Knots.Values();
        for(int i = k - Degree + 1; i <= k; i++)
        {
            Numeric alpha = (knot_value - knots[i]) / (knots[i + Degree] - knots[i]);
            auto point_i = this->ControlPoints[i];
            point_i.head<3>() *= point_i.w();
            auto point_im1 = this->ControlPoints[i - 1];
            point_im1.head<3>() *= point_im1.w();
            result.ControlPoints[i] = alpha * point_i + (1 - alpha) * point_im1;
            result.ControlPoints[i].head<3>() /= result.ControlPoints[i].w();
        }
        return result;
    }
}

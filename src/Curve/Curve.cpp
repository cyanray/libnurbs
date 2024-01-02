#include <libnurbs/Curve/Curve.hpp>
#include <libnurbs/Basis/BSplineBasis.hpp>
#include <cassert>

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
                tmp += point * basis(k, i);
            }
            result[k] = tmp;
        }
        return result;
    }

    Vec3 Curve::EvaluateDerivative(Numeric x, int order) const
    {
        return EvaluateAll(x, order).row(order);
    }

    MatX Curve::EvaluateAll(Numeric x, int order) const
    {
        static auto binomial = [](int n, int k) -> int
        {
            int result = 1;
            if (k > n) return 0;
            for (int i = 1; i <= k; i++)
            {
                result *= (n - i + 1);
                result /= i;
            }
            return result;
        };

        auto homo_ders = HomogeneousDerivative(x, order);
        MatX result = MatX::Zero(order + 1, 3);
        // Compute rational derivatives
        Numeric Wders0 = homo_ders[0].w();
        for (int k = 0; k <= order; k++)
        {
            Vec3 Aders = homo_ders[k].head<3>();
            for (int i = 1; i <= k; i++)
            {
                Numeric Wders = homo_ders[i].w();
                Aders -= binomial(k, i) * Wders * result.row(k - i);
            }
            result.row(k) = (Aders / Wders0);
        }
        return result;
    }
}

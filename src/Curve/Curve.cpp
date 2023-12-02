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

    Vec3 Curve::EvaluateDerivative(Numeric x, int order) const
    {
        // TODO: support NURBS derivative
        assert(x >= 0 && x <= 1);
        int index_span = Knots.FindSpanIndex(Degree, x);
        VecX basis = BSplineBasis::EvaluateDerivative(Degree, Knots.Values(), index_span, x, order);
        assert(basis.size() == Degree + 1);
        Vec4 result = Vec4::Zero();
        for (int i = 0; i <= Degree; i++)
        {
            auto point = ControlPoints[index_span - Degree + i];
            point.head<3>() *= basis(i);
            result += point;
        }
        return result.head<3>();
    }
}

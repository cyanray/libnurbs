#include <libnurbs/Curve/Curve.hpp>
#include <libnurbs/Basis/BSplineBasis.hpp>

namespace libnurbs
{
    Vec3 Curve::Evaluate(Numeric x) const
    {
        int index_span = Knots.FindSpanIndex(Degree, x);
        VecX basis = BSplineBasis::Evaluate(Degree, Knots.Values(), index_span, x);
        Vec4 result = Vec4::Zero();
        for (int i = 0; i <= Degree; i++)
        {
            result += basis(i) * ControlPoints[index_span - Degree + i];
        }
        return result.head<3>() / result(3);
    }
}

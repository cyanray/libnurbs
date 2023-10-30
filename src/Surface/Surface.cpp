#include <libnurbs/Surface/Surface.hpp>
#include <libnurbs/Basis/BSplineBasis.hpp>

using namespace std;


namespace libnurbs
{
    Vec3 Surface::Evaluate(Numeric u, Numeric v) const
    {
        assert(u >= 0 && u <= 1);
        assert(v >= 0 && v <= 1);
        int index_span_u = KnotsU.FindSpanIndex(DegreeU, u);
        int index_span_v = KnotsV.FindSpanIndex(DegreeV, v);
        VecX basis_u = BSplineBasis::Evaluate(DegreeU, KnotsU.Values(), index_span_u, u);
        VecX basis_v = BSplineBasis::Evaluate(DegreeV, KnotsV.Values(), index_span_v, v);
        assert(basis_u.size() == DegreeU + 1);
        assert(basis_v.size() == DegreeV + 1);
        Vec4 result = Vec4::Zero();
        for (int i = 0; i <= DegreeU; i++)
        {
            for (int j = 0; j <= DegreeV; j++)
            {
                auto point = ControlPoints.get(index_span_u - DegreeU + i, index_span_v - DegreeV + j);
                point.head<3>() *= point(3);
                result += basis_u(i) * basis_v(j) * point;
            }
        }
        return result.head<3>() / result(3);
    }
}
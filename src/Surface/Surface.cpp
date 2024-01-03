#include <libnurbs/Surface/Surface.hpp>
#include <libnurbs/Basis/BSplineBasis.hpp>

#include "libnurbs/Algorithm/MathUtils.hpp"

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
                auto point = ControlPoints.Get(index_span_u - DegreeU + i, index_span_v - DegreeV + j);
                point.head<3>() *= point(3);
                result += basis_u(i) * basis_v(j) * point;
            }
        }
        return result.head<3>() / result(3);
    }

    Vec3 Surface::EvaluateDerivative(Numeric u, Numeric v, int order_u, int order_v) const
    {
        assert(u >= 0 && u <= 1);
        assert(v >= 0 && v <= 1);
        int index_span_u = KnotsU.FindSpanIndex(DegreeU, u);
        int index_span_v = KnotsV.FindSpanIndex(DegreeV, v);
        VecX basis_u = BSplineBasis::EvaluateDerivative(DegreeU, KnotsU.Values(), index_span_u, u, order_u);
        VecX basis_v = BSplineBasis::EvaluateDerivative(DegreeV, KnotsV.Values(), index_span_v, v, order_v);
        assert(basis_u.size() == DegreeU + 1);
        assert(basis_v.size() == DegreeV + 1);
        Vec4 result = Vec4::Zero();
        for (int i = 0; i <= DegreeU; i++)
        {
            for (int j = 0; j <= DegreeV; j++)
            {
                auto point = ControlPoints.Get(index_span_u - DegreeU + i, index_span_v - DegreeV + j);
                point.head<3>() *= basis_u(i) * basis_v(j);
                result += point;
            }
        }
        return result.head<3>();
    }

    vvector<Vec3> Surface::EvaluateAll(Numeric u, Numeric v, int order_u, int order_v) const
    {
        auto homo_ders = HomogeneousDerivative(u, v, order_u, order_v);
        std::vector result(order_u + 1, std::vector<Vec3>(order_v + 1, Vec3::Zero()));
        Numeric Wders00 = homo_ders[0][0].w();

        for (int ou = 0; ou <= order_u; ++ou)
        {
            for (int ov = 0; ov <= order_v; ++ov)
            {
                Vec3 Aders = homo_ders[ou][ov].head<3>();
                for (int i = 1; i <= ou; ++i)
                {
                    Numeric Wders = homo_ders[i][0].w();
                    Aders -= Binomial(ou, i) * Wders * result[ou - i][ov];
                }

                for (int j = 1; j <= ov; ++j)
                {
                    Numeric Wders = homo_ders[0][j].w();
                    Aders -= Binomial(ov, j) * Wders * result[ou][ov - j];
                }

                for (int i = 1; i <= ou; ++ou)
                {
                    const int bi = Binomial(ou, i);
                    for (int j = 1; j <= ov; ++ov)
                    {
                        Numeric Wders = homo_ders[i][j].w();
                        Aders -= bi * Binomial(ov, j) * Wders * result[ou - i][ov - j];
                    }
                }
                result[ou][ov] = Aders / Wders00;
            }
        }

        return result;
    }

    vvector<Vec4> Surface::HomogeneousDerivative(Numeric u, Numeric v, int order_u, int order_v) const
    {
        assert(u >= 0 && u <= 1);
        assert(v >= 0 && v <= 1);
        int index_span_u = KnotsU.FindSpanIndex(DegreeU, u);
        int index_span_v = KnotsV.FindSpanIndex(DegreeV, v);
        MatX basis_u = BSplineBasis::EvaluateAll(DegreeU, KnotsU.Values(), index_span_u, u, order_u);
        MatX basis_v = BSplineBasis::EvaluateAll(DegreeV, KnotsV.Values(), index_span_v, v, order_v);
        assert(basis_u.size() == DegreeU + 1);
        assert(basis_v.size() == DegreeV + 1);

        std::vector result(order_u + 1, std::vector<Vec4>(order_v + 1, Vec4::Zero()));
        for (int k = 0; k <= order_u; ++k)
        {
            for (int l = 0; l <= order_v; ++l)
            {
                for (int i = 0; i <= DegreeU; i++)
                {
                    for (int j = 0; j <= DegreeV; j++)
                    {
                        int index_u = index_span_u - DegreeU + i;
                        int index_v = index_span_v - DegreeV + j;
                        auto point = ControlPoints.Get(index_u, index_v);
                        point.head<3>() *= point(3);
                        result[k][l] += point * basis_u(k, i) * basis_v(l, j);
                    }
                }
            }
        }
        return result;
    }
}

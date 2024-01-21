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

    Grid<Vec3> Surface::EvaluateAll(Numeric u, Numeric v, int order_u, int order_v) const
    {
        auto homo_ders = HomogeneousDerivative(u, v, order_u, order_v);
        Grid<Vec3> result(order_u + 1, order_v + 1, Vec3::Zero());

        Numeric Wders00 = homo_ders.Get(0, 0).w();

        for (int ou = 0; ou <= order_u; ++ou)
        {
            for (int ov = 0; ov <= order_v; ++ov)
            {
                Vec3 Aders = homo_ders.Get(ou, ov).head<3>();
                for (int i = 1; i <= ou; ++i)
                {
                    Numeric Wders = homo_ders.Get(i, 0).w();
                    Aders -= Binomial(ou, i) * Wders * result.Get(ou - i, ov);
                }

                for (int j = 1; j <= ov; ++j)
                {
                    Numeric Wders = homo_ders.Get(0, j).w();
                    Aders -= Binomial(ov, j) * Wders * result.Get(ou, ov - j);
                }

                for (int i = 1; i <= ou; ++ou)
                {
                    const int bi = Binomial(ou, i);
                    for (int j = 1; j <= ov; ++ov)
                    {
                        Numeric Wders = homo_ders.Get(i, j).w();
                        Aders -= bi * Binomial(ov, j) * Wders * result.Get(ou - i, ov - j);
                    }
                }
                result.Get(ou, ov) = Aders / Wders00;
            }
        }

        return result;
    }

    std::pair<Numeric, Numeric> Surface::SearchParameter(const Vec3& point) const
    {
        using Vec2 = Eigen::Vector2<Numeric>;
        using Mat2x2 = Eigen::Matrix<Numeric, 2, 2>;

        constexpr static int MAX_ITERATION_COUNT = 512;
        constexpr static Numeric EPSILON = 1e-18;

        auto Ri = [&point, this](Numeric u, Numeric v) -> Vec3
        {
            return Evaluate(u, v) - point;
        };

        auto Ki = [this, &Ri](Numeric u, Numeric v) -> Vec2
        {
            auto ri = Ri(u, v);
            auto Su = EvaluateDerivative(u, v, 1, 0);
            auto Sv = EvaluateDerivative(u, v, 0, 1);
            return Vec2{ri.dot(Su), ri.dot(Sv)};
        };

        auto Ji = [this, &Ri](Numeric u, Numeric v) -> Mat2x2
        {
            auto Su = EvaluateDerivative(u, v, 1, 0);
            auto Sv = EvaluateDerivative(u, v, 0, 1);
            auto Suu = EvaluateDerivative(u, v, 2, 0);
            auto Svv = EvaluateDerivative(u, v, 0, 2);
            auto Suv = EvaluateDerivative(u, v, 1, 1);
            auto ri = Ri(u, v);
            auto a = Su.dot(Su) + ri.dot(Suu);
            auto c = Su.dot(Sv) + ri.dot(Suv);
            auto d = Sv.dot(Sv) + ri.dot(Svv);
            Mat2x2 result;
            result << a, c, c, d;
            return result;
        };

        Numeric u_last = 0.5, v_last = 0.5;
        Vec3 res = Ri(u_last, v_last);

        // coefficient *s* is used to prevent uv_last from going out of range
        Numeric s = 1.0;

        int count = 0;
        while (res.norm() >= EPSILON && (count++ < MAX_ITERATION_COUNT))
        {
            auto k = Ki(u_last, v_last);
            auto j = Ji(u_last, v_last);

            Vec2 uv_last{u_last, v_last};
            uv_last += -j.inverse() * k * s;

            u_last = uv_last.x();
            v_last = uv_last.y();

            if (u_last < 0 || u_last > 1 || v_last < 0 || v_last > 1)
            {
                s /= 2;
                u_last = u_last > 1.0 ? 1.0 : u_last;
                u_last = u_last < 0.0 ? 0.0 : u_last;
                v_last = v_last > 1.0 ? 1.0 : v_last;
                v_last = v_last < 0.0 ? 0.0 : v_last;
            }
            res = Ri(u_last, v_last);
        }
        return {u_last, v_last};
    }

    Grid<Vec4> Surface::HomogeneousDerivative(Numeric u, Numeric v, int order_u, int order_v) const
    {
        assert(u >= 0 && u <= 1);
        assert(v >= 0 && v <= 1);
        int index_span_u = KnotsU.FindSpanIndex(DegreeU, u);
        int index_span_v = KnotsV.FindSpanIndex(DegreeV, v);
        MatX basis_u = BSplineBasis::EvaluateAll(DegreeU, KnotsU.Values(), index_span_u, u, order_u);
        MatX basis_v = BSplineBasis::EvaluateAll(DegreeV, KnotsV.Values(), index_span_v, v, order_v);
        assert(basis_u.size() == DegreeU + 1);
        assert(basis_v.size() == DegreeV + 1);

        Grid<Vec4> result(order_u + 1, order_v + 1, Vec4::Zero());

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
                        result.Get(k, l) += point * basis_u(k, i) * basis_v(l, j);
                    }
                }
            }
        }
        return result;
    }
}

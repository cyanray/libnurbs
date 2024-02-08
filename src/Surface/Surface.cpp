#include <libnurbs/Surface/Surface.hpp>
#include <libnurbs/Basis/BSplineBasis.hpp>

#include "libnurbs/Algorithm/DegreeElevate.hpp"
#include "libnurbs/Algorithm/KnotRemoval.hpp"
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
        int index_pre_u = index_span_u - DegreeU;
        int index_pre_v = index_span_v - DegreeV;

        for (int j = 0; j <= DegreeV; j++)
        {
            int index_v = index_pre_v + j;
            Vec4 tmp = Vec4::Zero();
            for (int i = 0; i <= DegreeU; i++)
            {
                int index_u = index_pre_u + i;
                auto point = ToHomo(ControlPoints.Get(index_u, index_v));
                tmp.noalias() += basis_u(i) * point;
            }
            result.noalias() += basis_v(j) * tmp;
        }
        return result.head<3>() / result.w();
    }

    Vec3 Surface::EvaluateDerivative(Numeric u, Numeric v, int order_u, int order_v) const
    {
        return EvaluateAll(u, v, order_u, order_v).Get(order_u, order_v);
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
                    Aders.noalias() -= Binomial(ou, i) * Wders * result.Get(ou - i, ov);
                }

                for (int j = 1; j <= ov; ++j)
                {
                    Numeric Wders = homo_ders.Get(0, j).w();
                    Aders.noalias() -= Binomial(ov, j) * Wders * result.Get(ou, ov - j);
                }

                for (int i = 1; i <= ou; ++i)
                {
                    const int bi = Binomial(ou, i);
                    for (int j = 1; j <= ov; ++j)
                    {
                        Numeric Wders = homo_ders.Get(i, j).w();
                        Aders.noalias() -= bi * Binomial(ov, j) * Wders * result.Get(ou - i, ov - j);
                    }
                }
                result.Get(ou, ov) = Aders / Wders00;
            }
        }

        return result;
    }

    std::pair<Numeric, Numeric> Surface::SearchParameter(const Vec3& point, Numeric init_u, Numeric init_v,
                                                         Numeric epsion, Numeric max_iteration_count) const
    {
        using Vec2 = Eigen::Vector2<Numeric>;
        using Mat2x2 = Eigen::Matrix<Numeric, 2, 2>;

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

        Numeric u_last = init_u, v_last = init_v;
        Numeric u_res = 1, v_res = 1;

        // coefficient *s* is used to prevent uv_last from going out of range
        Numeric s = 1.0;

        int count = 0;
        while ((u_res >= epsion || v_res >= epsion) && (count++ < max_iteration_count))
        {
            auto k = Ki(u_last, v_last);
            auto j = Ji(u_last, v_last);

            Vec2 uv_last{u_last, v_last};
            uv_last.noalias() += -j.inverse() * k * s;

            u_res = abs(uv_last.x() - u_last);
            v_res = abs(uv_last.y() - v_last);

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
        }
        return {u_last, v_last};
    }

    Surface Surface::InsertKnotU(Numeric knot_value) const
    {
        Surface result{*this};
        int k = result.KnotsU.InsertKnot(knot_value);
        result.ControlPoints.InsertU(k, Vec4::Zero());
        const auto& knots = this->KnotsU.Values();
        vector<Numeric> alpha_list(result.DegreeU);
        for (int i = k - DegreeU + 1; i <= k; i++)
        {
            int idx = i - k + DegreeU - 1;
            alpha_list[idx] = (knot_value - knots[i]) / (knots[i + DegreeU] - knots[i]);
        }

        for (int j = 0; j < result.ControlPoints.VCount; j++)
        {
            for (int i = k - DegreeU + 1; i <= k; i++)
            {
                Vec4 point_i = ToHomo(this->ControlPoints.Get(i, j));
                Vec4 point_im1 = ToHomo(this->ControlPoints.Get(i - 1, j));
                int idx = i - k + DegreeU - 1;
                Numeric alpha = alpha_list[idx];
                Vec4& point_new = result.ControlPoints.Get(i, j);
                point_new = FromHomo(alpha * point_i + (1 - alpha) * point_im1);
            }
        }
        return result;
    }

    Surface Surface::InsertKnotU(Numeric knot_value, int times) const
    {
        std::vector list(times, knot_value);
        return InsertKnotU(list);
    }

    Surface Surface::InsertKnotU(std::span<Numeric> knots_to_insert) const
    {
        Surface result{*this};
        for (auto knot: knots_to_insert)
        {
            result = result.InsertKnotU(knot);
        }
        return result;
    }

    Surface Surface::InsertKnotV(Numeric knot_value) const
    {
        Surface result{*this};
        int k = result.KnotsV.InsertKnot(knot_value);
        result.ControlPoints.InsertV(k, Vec4::Zero());
        const auto& knots = this->KnotsV.Values();
        vector<Numeric> alpha_list(result.DegreeV);
        for (int i = k - DegreeV + 1; i <= k; i++)
        {
            int idx = i - k + DegreeV - 1;
            alpha_list[idx] = (knot_value - knots[i]) / (knots[i + DegreeV] - knots[i]);
        }

        for (int j = 0; j < result.ControlPoints.UCount; j++)
        {
            for (int i = k - DegreeV + 1; i <= k; i++)
            {
                Vec4 point_i = ToHomo(this->ControlPoints.Get(j, i));
                Vec4 point_im1 = ToHomo(this->ControlPoints.Get(j, i - 1));

                int idx = i - k + DegreeV - 1;
                Numeric alpha = alpha_list[idx];
                Vec4& point_new = result.ControlPoints.Get(j, i);
                point_new = FromHomo(alpha * point_i + (1 - alpha) * point_im1);
            }
        }
        return result;
    }

    Surface Surface::InsertKnotV(Numeric knot_value, int times) const
    {
        std::vector list(times, knot_value);
        return InsertKnotV(list);
    }

    Surface Surface::InsertKnotV(std::span<Numeric> knots_to_insert) const
    {
        Surface result{*this};
        for (auto knot: knots_to_insert)
        {
            result = result.InsertKnotV(knot);
        }
        return result;
    }

    std::tuple<Surface, int> Surface::RemoveKnotU(Numeric knot_remove, int times, Numeric tolerance) const
    {
        Surface result{*this};
        auto& points_ref = result.ControlPoints;
        points_ref = {points_ref.UCount - times, points_ref.VCount};
        int t = 0;
        for (int it_v = 0; it_v < ControlPoints.VCount; ++it_v)
        {
            auto knots = result.KnotsU;
            auto points = this->ControlPoints.GetV(it_v);
            int t_tmp = KnotRemoval(knots, points, result.DegreeU, knot_remove, times, tolerance);
            if (t_tmp == 0 || it_v != 0 && t != t_tmp) return {*this, 0};
            if (it_v == ControlPoints.VCount - 1) result.KnotsU = knots;
            points_ref.SetV(it_v, points);
            t = t_tmp;
        }
        return {result, t};
    }

    std::tuple<Surface, int> Surface::RemoveKnotV(Numeric knot_remove, int times, Numeric tolerance) const
    {
        Surface result{*this};
        auto& points_ref = result.ControlPoints;
        points_ref = {points_ref.UCount, points_ref.VCount - times};
        int t = 0;
        for (int it_u = 0; it_u < ControlPoints.UCount; ++it_u)
        {
            auto knots = result.KnotsV;
            auto points = this->ControlPoints.GetU(it_u);
            int t_tmp = KnotRemoval(knots, points, result.DegreeV, knot_remove, times, tolerance);
            if (t_tmp == 0 || it_u != 0 && t != t_tmp) return {*this, 0};
            if (it_u == ControlPoints.UCount - 1) result.KnotsV = knots;
            points_ref.SetU(it_u, points);
            t = t_tmp;
        }
        return {result, t};
    }

    Surface Surface::ElevateDegreeU(int times) const
    {
        Surface result{*this};
        auto& points_ref = result.ControlPoints;
        points_ref = {points_ref.UCount + times, points_ref.VCount};
        for (int it_v = 0; it_v < ControlPoints.VCount; ++it_v)
        {
            auto degree = result.DegreeU;
            auto knots = result.KnotsU;
            auto points = this->ControlPoints.GetV(it_v);
            DegreeElevate(knots, points, degree, times);
            if (it_v == ControlPoints.VCount - 1)
            {
                result.KnotsU = knots;
                result.DegreeU = degree;
            }
            points_ref.SetV(it_v, points);
        }
        return result;
    }

    Surface Surface::ElevateDegreeV(int times) const
    {
        Surface result{*this};
        auto& points_ref = result.ControlPoints;
        points_ref = {points_ref.UCount, points_ref.VCount + times};
        for (int it_u = 0; it_u < ControlPoints.UCount; ++it_u)
        {
            auto degree = result.DegreeV;
            auto knots = result.KnotsV;
            auto points = this->ControlPoints.GetU(it_u);
            DegreeElevate(knots, points, degree, times);
            if (it_u == ControlPoints.UCount - 1)
            {
                result.KnotsV = knots;
                result.DegreeV = degree;
            }
            points_ref.SetU(it_u, points);
        }
        return result;
    }

    Grid<Vec4> Surface::HomogeneousDerivative(Numeric u, Numeric v, int order_u, int order_v) const
    {
        assert(u >= 0 && u <= 1);
        assert(v >= 0 && v <= 1);
        int index_span_u = KnotsU.FindSpanIndex(DegreeU, u);
        int index_span_v = KnotsV.FindSpanIndex(DegreeV, v);
        MatX basis_u = BSplineBasis::EvaluateAll(DegreeU, KnotsU.Values(), index_span_u, u, order_u);
        MatX basis_v = BSplineBasis::EvaluateAll(DegreeV, KnotsV.Values(), index_span_v, v, order_v);
        assert(basis_u.cols() == DegreeU + 1);
        assert(basis_v.cols() == DegreeV + 1);

        int index_pre_u = index_span_u - DegreeU;
        int index_pre_v = index_span_v - DegreeV;

        Grid<Vec4> result(order_u + 1, order_v + 1, Vec4::Zero());

        for (int l = 0; l <= order_v; ++l)
        {
            for (int k = 0; k <= order_u; ++k)
            {
                for (int j = 0; j <= DegreeV; j++)
                {
                    int index_v = index_pre_v + j;
                    Vec4 tmp = Vec4::Zero();
                    for (int i = 0; i <= DegreeU; i++)
                    {
                        int index_u = index_pre_u + i;
                        auto point = ToHomo(ControlPoints.Get(index_u, index_v));
                        tmp.noalias() += basis_u(k, i) * point;
                    }
                    result.Get(k, l).noalias() += basis_v(l, j) * tmp;
                }
            }
        }
        return result;
    }
}

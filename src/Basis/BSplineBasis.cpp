#include <libnurbs/Basis/BSplineBasis.hpp>
#include <libnurbs/Core/KnotVector.hpp>
#include <libnurbs/Algorithm/MathUtils.hpp>

namespace libnurbs
{
    // from https://github.com/pradeep-pyro/tinynurbs
    VecX BSplineBasis::Evaluate(int degree, const vector<Numeric>& knots, int index_span, Numeric x)
    {
        VecX result = VecX::Zero(degree + 1);
        VecX left = VecX::Zero(degree + 1);
        VecX right = VecX::Zero(degree + 1);
        Numeric saved, temp;
        result[0] = 1.0;
        for (int j = 1; j <= degree; j++)
        {
            left[j] = (x - knots[index_span + 1 - j]);
            right[j] = knots[index_span + j] - x;
            saved = 0.0;
            for (int r = 0; r < j; r++)
            {
                temp = result[r] / (right[r + 1] + left[j - r]);
                result[r] = saved + right[r + 1] * temp;
                saved = left[j - r] * temp;
            }
            result[j] = saved;
        }
        return result;
    }


    VecX BSplineBasis::Evaluate(int degree, const KnotVector& knot_vec, Numeric x)
    {
        auto& knots = knot_vec.Values();
        auto index_span = knot_vec.FindSpanIndex(degree, x);
        return BSplineBasis::Evaluate(degree, knots, index_span, x);
    }


    VecX BSplineBasis::EvaluateDerivative(int degree, const KnotVector& knot_vec, Numeric x, int order)
    {
        auto& knots = knot_vec.Values();
        auto index_span = knot_vec.FindSpanIndex(degree, x);
        return BSplineBasis::EvaluateDerivative(degree, knots, index_span, x, order);
    }


    VecX BSplineBasis::EvaluateDerivative(int degree, const std::vector<Numeric>& knots,
                                          int index_span, Numeric x, int order)
    {
        return EvaluateAll(degree, knots, index_span, x, order).row(order);
    }

    // from https://github.com/pradeep-pyro/tinynurbs
    MatX BSplineBasis::EvaluateAll(int degree, const vector<Numeric>& knots, int index_span, Numeric x, int order)
    {
        std::vector<Numeric> left, right;
        left.resize(degree + 1, 0.0);
        right.resize(degree + 1, 0.0);
        Numeric saved, temp;

        MatX ndu = MatX::Zero(degree + 1, degree + 1);
        ndu(0, 0) = 1.0;

        for (int j = 1; j <= degree; j++)
        {
            left[j] = x - knots[index_span + 1 - j];
            right[j] = knots[index_span + j] - x;
            saved = 0.0;
            for (int r = 0; r < j; r++)
            {
                // Lower triangle
                ndu(j, r) = right[r + 1] + left[j - r];
                temp = ndu(r, j - 1) / ndu(j, r);
                // Upper triangle
                ndu(r, j) = saved + right[r + 1] * temp;
                saved = left[j - r] * temp;
            }
            ndu(j, j) = saved;
        }

        MatX result = MatX::Zero(order + 1, degree + 1);
        result.row(0) = ndu.col(degree);

        MatX a = MatX::Zero(2, degree + 1);
        for (int r = 0; r <= degree; r++)
        {
            int s1 = 0, s2 = 1;
            a(0, 0) = 1.0;
            for (int k = 1; k <= order; k++)
            {
                Numeric d = 0.0;
                int rk = r - k;
                int pk = degree - k;
                int j1, j2;

                if (r >= k)
                {
                    a(s2, 0) = a(s1, 0) / ndu(pk + 1, rk);
                    d = a(s2, 0) * ndu(rk, pk);
                }

                j1 = (rk >= -1) ? 1 : -rk;
                j2 = (r - 1 <= pk) ? k - 1 : degree - r;

                for (int j = j1; j <= j2; j++)
                {
                    a(s2, j) = (a(s1, j) - a(s1, j - 1)) / ndu(pk + 1, rk + j);
                    d += a(s2, j) * ndu(rk + j, pk);
                }

                if (r <= pk)
                {
                    a(s2, k) = -a(s1, k - 1) / ndu(pk + 1, r);
                    d += a(s2, k) * ndu(r, pk);
                }

                result(k, r) = d;
                std::swap(s1, s2);
            }
        }

        Numeric fac = degree;
        for (int k = 1; k <= order; k++)
        {
            result.row(k) *= fac;
            fac *= (degree - k);
        }

        return result;
    }

    MatX BSplineBasis::EvaluateAll(int degree, const KnotVector& knot_vec, Numeric x, int order)
    {
        auto& knots = knot_vec.Values();
        auto index_span = knot_vec.FindSpanIndex(degree, x);
        return BSplineBasis::EvaluateAll(degree, knots, index_span, x, order);
    }

}
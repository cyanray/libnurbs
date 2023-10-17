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

//    VecX BSplineBasis::EvaluateDerivative(int degree, const KnotVector& knot_vec, Numeric x, int order)
//    {
//        auto& knots = knot_vec.Values();
//        auto index_span = knot_vec.FindSpanIndex(degree, x);
//        auto calc_a = [p = degree, i = index_span, &u = knots](auto self, int k, int j) -> Numeric
//        {
//            if (k == 0 && j == 0) return 1;
//            if (j == 0)
//            {
//                Numeric m = (u[i + p - k + 1] - u[i]);
//                return (m == 0) ? 0.0 : self(self, k - 1, 0) / m;
//            }
//            if (j == k)
//            {
//                Numeric m = (u[i + p + 1] - u[i + k]);
//                return (m == 0) ? 0.0 : (-self(self, k - 1, k - 1)) / m;
//            }
//            Numeric m = (u[i + p + j - k + 1] - u[i + j]);
//            return (m == 0) ? 0.0 : (self(self, k - 1, j) - self(self, k - 1, j - 1)) / m;
//        };
//        //VecX result = VecX::Zero(degree + 1);
//        VecX N = BSplineBasis::Evaluate(degree - order, knot_vec, x);
//        Numeric R = 0.0;
//        for (int j = 0; j <= order; ++j)
//        {
//            Numeric a = calc_a(calc_a, order, j);
//            R += a * N(j);
//        }
//        Numeric t = (Numeric) Factorial(degree) / (Numeric) Factorial(degree - order);
//        R *= t;
//        return {};
//    }

    MatX BSplineBasis::EvaluateDerivative(int degree, const KnotVector& knot_vec, Numeric x, int order)
    {
        auto& knots = knot_vec.Values();
        auto index_span = knot_vec.FindSpanIndex(degree, x);
        return BSplineBasis::EvaluateDerivative(degree, knots, index_span, x, order);
    }

    // from https://github.com/pradeep-pyro/tinynurbs
    MatX BSplineBasis::EvaluateDerivative(int degree, const std::vector<Numeric>& knots, int index_span, Numeric x,
                                          int order)
    {
        std::vector<Numeric> left, right;
        left.resize(degree + 1, 0.0);
        right.resize(degree + 1, 0.0);
        Numeric saved = 0.0, temp = 0.0;

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

        MatX ders = MatX::Zero(order + 1, degree + 1);

        for (int j = 0; j <= degree; j++)
        {
            ders(0, j) = ndu(j, degree);
        }

        MatX a = MatX::Zero(2, degree + 1);

        for (int r = 0; r <= degree; r++)
        {
            int s1 = 0;
            int s2 = 1;
            a(0, 0) = 1.0;

            for (int k = 1; k <= order; k++)
            {
                Numeric d = 0.0;
                int rk = r - k;
                int pk = degree - k;
                int j1 = 0;
                int j2 = 0;

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

                ders(k, r) = d;

                std::swap(s1, s2);
            }
        }

        Numeric fac = degree;
        for (int k = 1; k <= order; k++)
        {
            for (int j = 0; j <= degree; j++)
            {
                ders(k, j) *= fac;
            }
            fac *= (degree - k);
        }

        return ders;
    }

}
#include <libnurbs/Basis/BSplineBasis.hpp>
#include <libnurbs/Algorithm/MathUtils.hpp>

namespace libnurbs
{
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

    VecX BSplineBasis::EvaluateDerivative(int degree, const vector<Numeric>& knots,
                                          int index_span, Numeric x, int order)
    {
        Numeric t = (Numeric) Factorial(degree) / (Numeric) Factorial(degree - order);
        Numeric a_cache;
        auto calc_a = [p = degree, i = index_span, &u = knots](auto self, int k, int j) -> Numeric
        {
            if (k == 0 && j == 0) return 1;
            if (j == 0)
            {
                return self(self, k - 1, 0) / (u[i + p - k + 1] - u[i]);
            }
            if (j == k)
            {
                return (-self(self, k - 1, k - 1)) / (u[i + p + 1] - u[i + k]);
            }
            return (self(self, k - 1, j) - self(self, k - 1, j - 1)) / (u[i + p + j - k + 1] - u[i + j]);
        };
        VecX result = VecX::Zero(degree + 1);
        for (int j = 0; j <= order; ++j)
        {
            VecX N = BSplineBasis::Evaluate(degree - order, knots, index_span + j, x);
            Numeric a = calc_a(calc_a, order, j);
            result += a * N;
        }
        return t * result;
    }
}
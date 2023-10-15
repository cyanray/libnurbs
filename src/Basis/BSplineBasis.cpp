#include <libnurbs/Basis/BSplineBasis.hpp>

namespace libnurbs
{
    VecX BSplineBasis::Evaluate(int degree, const KnotVector& knot_vec, Numeric x)
    {
        VecX result; result.resize(degree + 1); result.fill(0.0);
        VecX left; left.resize(degree + 1); left.fill(0.0);
        VecX right; right.resize(degree + 1); right.fill(0.0);
        int span = knot_vec.FindSpanIndex(degree, x);
        auto& knots = knot_vec.Values();
        Numeric saved, temp;
        result[0] = 1.0;
        for (int j = 1; j <= degree; j++)
        {
            left[j] = (x - knots[span + 1 - j]);
            right[j] = knots[span + j] - x;
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
}
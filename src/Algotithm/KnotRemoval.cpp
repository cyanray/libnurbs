#include "libnurbs/Algorithm/KnotRemoval.hpp"
#include "libnurbs/Core/KnotVector.hpp"

namespace
{
    using namespace libnurbs;
    constexpr Numeric CalcTOL(const vector<Vec4>& points, Numeric epsilon)
    {
        Numeric w_min = 1e12;
        Numeric dist_max = 0.0;
        for (const auto& point: points)
        {
            w_min = std::min(w_min, point.w());
            dist_max = std::max(dist_max, point.head<3>().norm());
        }
        return (epsilon * w_min) / (1 + dist_max);
    };
}

namespace libnurbs
{
    int KnotRemoval(KnotVector& knot_vector,
                    std::vector<Vec4>& points,
                    int degree, Numeric knot_remove, int times, Numeric tolerance)
    {
        auto& knots = knot_vector.Values();
        int s = knot_vector.GetMultiplicity(knot_remove);
        int r = knot_vector.FindSpanIndex(degree, knot_remove);

        double tol = CalcTOL(points, tolerance);

        int n = (int)points.size() - 1;
        int m = n + degree + 1;
        int order = degree + 1;
        int last = r - s;
        int first = r - degree;

        std::vector<Vec4> temp(2 * degree + 1);
        int t;
        for (t = 0; t < times; t++)
        {
            int off = first - 1;
            temp[0] = points[off];
            temp[last + 1 - off] = points[last + 1];
            int i = first;
            int j = last;
            int ii = 1;
            int jj = last - off;

            while (j - i >= t)
            {
                double alphai = (knot_remove - knots[i]) / (knots[i + order + t] - knots[i]);
                double alphaj = (knot_remove - knots[j - t]) / (knots[j + order] - knots[j - t]);

                temp[ii] = (points[i] - (1.0 - alphai) * temp[ii - 1]) / alphai;
                temp[jj] = (points[j] - alphaj * temp[jj + 1]) / (1.0 - alphaj);

                i = i + 1;
                ii = ii + 1;

                j = j - 1;
                jj = jj - 1;
            }

            Numeric dist;
            if (j - i < t)
            {
                dist = (ToHomo(temp[ii - 1]) - ToHomo(temp[jj + 1])).norm();
            }
            else
            {
                double alphai = (knot_remove - knots[i]) / (knots[i + order + t] - knots[i]);
                dist = (alphai * ToHomo(temp[ii + t + 1]) + (1.0 - alphai) * ToHomo(temp[ii - 1])).norm();
            }
            if (dist > tol) break;

            i = first;
            j = last;
            while (j - i > t)
            {
                points[i] = temp[i - off];
                points[j] = temp[j - off];
                i = i + 1;
                j = j - 1;
            }

            first = first - 1;
            last = last + 1;
        }

        // no remove operation succeeded
        if (t == 0) return t;

        int j = (2 * r - s - degree) / 2, i = j;
        for (int k = 1; k < t; k++)
        {
            i += (k % 2 == 1) ? 1 : 0;
            j -= (k % 2 == 1) ? 0 : 1;
        }

        /* Update control points */
        for (int k = i + 1; k <= n; k++, j++)
        {
            points[j] = points[k];
        }
        points.resize(points.size() - times);

        /* Update knot vector */
        for (int k = r + 1; k <= m; k++)
        {
            knots[k - times] = knots[k];
        }
        knots.resize(knots.size() - times);

        return t;
    }
}

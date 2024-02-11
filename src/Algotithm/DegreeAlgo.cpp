#include "libnurbs/Algorithm/DegreeAlgo.hpp"

#include "libnurbs/Algorithm/MathUtils.hpp"
#include "libnurbs/Core/KnotVector.hpp"

namespace libnurbs
{
    void NurbsDegreeElevation(KnotVector& knot_vector, std::vector<Vec4>& points, int& degree, int times)
    {
        assert(times >= 1);

        // Convert to Homogeneous Coordinate
        for (auto& point: points)
        {
            point.head<3>() *= point.w();
        }

        auto& knots = knot_vector.Values();

        int n = (int)points.size() - 1;
        int m = n + degree + 1;
        int ph = degree + times;
        int ph2 = ph / 2;

        std::vector bezalfs(degree + times + 1, std::vector<Numeric>(degree + 1));
        bezalfs[0][0] = bezalfs[ph][degree] = 1.0;

        for (int i = 1; i <= ph2; i++)
        {
            double inv = 1.0 / Binomial(ph, i);
            int mpi = std::min(degree, i);

            for (int j = std::max(0, i - times); j <= mpi; j++)
            {
                bezalfs[i][j] = inv * Binomial(degree, j) * Binomial(times, i - j);
            }
        }

        for (int i = ph2 + 1; i <= ph - 1; i++)
        {
            int mpi = std::min(degree, i);
            for (int j = std::max(0, i - times); j <= mpi; j++)
            {
                bezalfs[i][j] = bezalfs[ph - i][degree - j];
            }
        }

        int mh = ph;
        int kind = ph + 1;
        int r = -1;
        int a = degree;
        int b = degree + 1;
        int cind = 1;
        double ua = knots[0];

        std::vector<Vec4> Qw((n + times) * 2);
        Qw[0] = points[0];

        std::vector<Numeric> Uh((n + times) * 2 + ph + 1);
        std::fill(Uh.begin(), Uh.begin() + ph + 1, ua);

        std::vector<Vec4> bpts(points.begin(), points.begin() + degree + 1);
        std::vector<Vec4> nextbpts(degree - 1);
        std::vector<Vec4> ebpts(degree + times + 1);

        while (b < m)
        {
            int ii = b;
            while (b < m && Approx(knots[b], knots[b + 1]))
            {
                b = b + 1;
            }
            int mul = b - ii + 1;
            mh += mul + times;
            double ub = knots[b];

            int oldr = r;
            r = degree - mul;

            int lbz = oldr > 0 ? (oldr + 2) / 2 : 1;
            int rbz = r > 0 ? ph - (r + 1) / 2 : ph;

            if (r > 0)
            {
                double numer = ub - ua;
                std::vector<Numeric> alfs(degree - 1); /// ??: move to out of the loop scope?
                for (int k = degree; k > mul; k--)
                {
                    alfs[k - mul - 1] = numer / (knots[a + k] - ua);
                }
                for (int j = 1; j <= r; j++)
                {
                    int save = r - j;
                    int s = mul + j;

                    for (int k = degree; k >= s; k--)
                    {
                        bpts[k] = alfs[k - s] * bpts[k] + (1.0 - alfs[k - s]) * bpts[k - 1];
                    }
                    nextbpts[save] = bpts[degree];
                }
            }

            for (int i = lbz; i <= ph; i++)
            {
                ebpts[i] = Vec4::Zero();
                int mpi = std::min(degree, i);
                for (int j = std::max(0, i - times); j <= mpi; j++)
                {
                    ebpts[i].noalias() += bezalfs[i][j] * bpts[j];
                }
            }

            if (oldr > 1)
            {
                int first = kind - 2;
                int last = kind;
                double den = ub - ua;
                double bet = (ub - Uh[kind - 1]) / den;

                for (int tr = 1; tr < oldr; tr++)
                {
                    int i = first;
                    int j = last;
                    int kj = j - kind + 1;

                    while (j - i > tr)
                    {
                        if (i < cind)
                        {
                            Numeric alf = (ub - Uh[i]) / (ua - Uh[i]);
                            Qw[i].noalias() = alf * Qw[i] + (1.0 - alf) * Qw[i - 1];
                        }

                        if (j >= lbz)
                        {
                            if (j - tr <= kind - ph + oldr)
                            {
                                double gam = (ub - Uh[j - tr]) / den;
                                ebpts[kj].noalias() = gam * ebpts[kj] + (1.0 - gam) * ebpts[kj + 1];
                            }
                            else
                            {
                                ebpts[kj].noalias() = bet * ebpts[kj] + (1.0 - bet) * ebpts[kj + 1];
                            }
                        }

                        i = i + 1;
                        j = j - 1;
                        kj = kj - 1;
                    } // end of: while (j - i > tr)

                    first =  first - 1;
                    last = last + 1;
                } // end of: for (int tr = 1; tr < oldr; tr++)

            } // end of: if (oldr > 1)

            if (a != degree)
            {
                for (int i = 0; i < ph - oldr; i++)
                {
                    Uh[kind] = ua;
                    kind = kind + 1;
                }
            }

            for (int j = lbz; j <= rbz; j++)
            {
                Qw[cind] = ebpts[j];
                cind = cind + 1;
            }

            if (b < m)
            {
                for (int j = 0; j < r; j++)
                {
                    bpts[j] = nextbpts[j];
                }
                for (int j = r; j <= degree; j++)
                {
                    bpts[j] = points[b - degree + j];
                }

                a = b;
                b = b + 1;
                ua = ub;
            }
            else
            {
                for (int i = 0; i <= ph; i++)
                {
                    Uh[kind + i] = ub;
                }
            }

        } // end of: while (b < m)

        int nh = mh - ph - 1;
        Qw.resize(nh + 1);
        Uh.resize(mh + 1);

        // Convert to Cartesian Coordinate
        for (auto& point: Qw)
        {
            point.head<3>() /= point.w();
        }

        degree = ph;
        knot_vector = KnotVector(Uh);
        points = Qw;
    }
}

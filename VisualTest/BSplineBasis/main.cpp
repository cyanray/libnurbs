#include <iostream>
#include <libnurbs/Basis/BSplineBasis.hpp>
#include <libnurbs/Core/KnotVector.hpp>
#include <Eigen/Dense>

#include <matplot/matplot.h>

#include <vector>

using namespace std;
using namespace libnurbs;
using namespace matplot;

int main()
{
    int degree = 3;
    KnotVector U{{0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0, 1.0}};
    VecX x = VecX::LinSpaced(30, 0.0, 1.0);
    MatX points;
    points.resize(degree + 1, x.size());
    for (int i = 0; i < x.size(); ++i)
    {
        VecX y = BSplineBasis::Evaluate(degree, U, x(i));
        points.col(i) = y;
    }
    vector<double> px{x.data(), x.data() + x.size()};
    for(int i = 0; i < points.rows(); ++i)
    {
        vector<double> py(px.size());
        for(int j = 0; j < px.size(); ++j)
        {
            py[j] = points(i, j);
        }

        plot(px, py, "ko");
        hold(on);
    }
    show();
    return 0;
}
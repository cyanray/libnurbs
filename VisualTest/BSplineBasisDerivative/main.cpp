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
    VecX x = VecX::LinSpaced(20, 0.0, 1.0);
    for (int i = 0; i < x.size(); ++i)
    {
        VecX y = BSplineBasis::EvaluateDerivative(degree, U.Values(), U.FindSpanIndex(degree, x(i)), x(i), 1);
        for (int j = 0; j < y.size(); ++j)
        {
            plot({x(i)}, {y(j)}, "ko");
            hold(on);
        }
    }

    cin.get();
    return 0;
}
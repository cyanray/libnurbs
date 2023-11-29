#include <iostream>
#include <libnurbs/Basis/BSplineBasis.hpp>
#include <libnurbs/Core/KnotVector.hpp>
#include <Eigen/Dense>

#include <matplot/matplot.h>

#include <vector>

#include <libnurbs/Curve/Curve.hpp>

using namespace std;
using namespace libnurbs;
using namespace matplot;

void CaseA()
{
    Curve curve;
    curve.Degree = 4;
    curve.Knots = KnotVector{{0.0, 0.0, 0.0, 0.0, 0.0, 0.333, 0.667, 1.0, 1.0, 1.0, 1.0, 1.0}};
    curve.ControlPoints =
    {
        {0.00, 0.00, 0.00, 1},
        {0.25, 0.00, 0.00, 1},
        {0.50, 0.00, 0.00, 1},
        {1.00, 0.00, 0.00, 1},
        {1.50, 0.00, 0.00, 1},
        {1.75, 0.00, 0.00, 1},
        {2.00, 0.00, 0.00, 1}
    };

    VecX s = VecX::LinSpaced(60, 0.0, 1.0);

    vector<double> px;
    px.reserve(s.size());
    vector<double> py;
    py.reserve(s.size());
    vector<double> pz;
    pz.reserve(s.size());

    vector<double> vx;
    vx.reserve(s.size());
    vector<double> vy;
    vy.reserve(s.size());

    for (auto x: s)
    {
        auto value = curve.Evaluate(x);
        auto de = curve.EvaluateDerivative(x, 2);
        px.push_back(value.x());
        py.push_back(value.y());
        pz.push_back(value.z());
        vx.push_back(0.001);
        vy.push_back(de.x());
    }
    plot(s, vy, "r.");
    show();
    quiver(px, py, vx,vy);
    show();
}


void CaseB()
{
    Curve curve;
    curve.Degree = 3;
    curve.Knots = KnotVector{{0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0, 1.0}};
    curve.ControlPoints =
    {
        {0.00, 0.00, 0.00, 1},
        {0.50, 0.00, 0.00, 1},
        {1.00, 0.00, 0.00, 1},
        {1.50, 0.00, 0.00, 1},
        {2.00, 0.00, 0.00, 1},
    };

    VecX s = VecX::LinSpaced(60, 0.0, 1.0);

    vector<double> px;
    px.reserve(s.size());
    vector<double> py;
    py.reserve(s.size());
    vector<double> pz;
    pz.reserve(s.size());

    vector<double> vx;
    vx.reserve(s.size());
    vector<double> vy;
    vy.reserve(s.size());

    for (auto x: s)
    {
        auto value = curve.Evaluate(x);
        auto de = curve.EvaluateDerivative(x, 1);
        px.push_back(value.x());
        py.push_back(value.y());
        pz.push_back(value.z());
        vx.push_back(0.001);
        vy.push_back(de.x());
    }
    plot(px, vy, "r.");
    show();
}

void CaseC()
{
    Curve curve;
    curve.Degree = 2;
    curve.Knots = KnotVector{{0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0}};
    curve.ControlPoints =
    {
        {0.00, 0.00, 0.00, 1},
        {0.50, 0.00, 0.00, 1},
        {1.50, 0.00, 0.00, 1},
        {2.00, 0.00, 0.00, 1},
    };

    VecX s = VecX::LinSpaced(60, 0.0, 1.0);

    vector<double> px;
    px.reserve(s.size());
    vector<double> py;
    py.reserve(s.size());
    vector<double> pz;
    pz.reserve(s.size());

    vector<double> vx;
    vx.reserve(s.size());
    vector<double> vy;
    vy.reserve(s.size());

    for (auto x: s)
    {
        auto value = curve.Evaluate(x);
        auto de = curve.EvaluateDerivative(x, 1);
        px.push_back(value.x());
        py.push_back(value.y());
        pz.push_back(value.z());
        vx.push_back(0.001);
        vy.push_back(de.x());
    }
    plot(px, vy, "r.");
    show();
}

int main()
{
    CaseA();
    // CaseB();
    // CaseC();

    return 0;
}

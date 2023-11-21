#include <iostream>
#include <libnurbs/Basis/BSplineBasis.hpp>
#include <libnurbs/Core/KnotVector.hpp>
#include <Eigen/Dense>

#include <matplot/matplot.h>

#include <vector>

#include <libnurbs/Core/ControlPointGrid.hpp>

#include <libnurbs/Surface/Surface.hpp>

using namespace std;
using namespace libnurbs;
using namespace matplot;

void CaseA()
{
    ControlPointGrid grid;
    grid.UCount = 3;
    grid.VCount = 2;
    grid.ControlPoints.emplace_back(1, 0, 0, 1);
    grid.ControlPoints.emplace_back(1, 1, 0, 1);
    grid.ControlPoints.emplace_back(0, 1, 0, 2);
    grid.ControlPoints.emplace_back(1, 0, 1, 1);
    grid.ControlPoints.emplace_back(1, 1, 1, 1);
    grid.ControlPoints.emplace_back(0, 1, 1, 2);

    Surface surface;
    surface.DegreeU = 2;
    surface.DegreeV = 1;
    surface.KnotsU = KnotVector{{0.0, 0.0, 0.0, 1.0, 1.0, 1.0}};
    surface.KnotsV = KnotVector{{0.0, 0.0, 1.0, 1.0}};
    surface.ControlPoints = grid;

    VecX s = VecX::LinSpaced(20, 0.0, 1.0);
    VecX t = VecX::LinSpaced(20, 0.0, 1.0);

    vector<double> px;
    px.reserve(s.size());
    vector<double> py;
    py.reserve(s.size());
    vector<double> pz;
    pz.reserve(s.size());

    for (auto x: s)
    {
        for (auto y: t)
        {
            auto value = surface.Evaluate(x, y);
            px.push_back(value.x());
            py.push_back(value.y());
            pz.push_back(value.z());
        }
    }
    plot3(px, py, pz, "r.");
    show();
}

void CaseB()
{
    ControlPointGrid grid;
    grid.UCount = 3;
    grid.VCount = 2;
    grid.ControlPoints.emplace_back(1, 0, 0, 1);
    grid.ControlPoints.emplace_back(1, 1, 0, 1);
    grid.ControlPoints.emplace_back(0, 1, 0, 2);
    grid.ControlPoints.emplace_back(1, 0, 1, 1);
    grid.ControlPoints.emplace_back(1, 1, 1, 1);
    grid.ControlPoints.emplace_back(0, 1, 1, 2);

    Surface surface;
    surface.DegreeU = 2;
    surface.DegreeV = 1;
    surface.KnotsU = KnotVector{{0.0, 0.0, 0.0, 1.0, 1.0, 1.0}};
    surface.KnotsV = KnotVector{{0.0, 0.0, 1.0, 1.0}};
    surface.ControlPoints = grid;

    VecX s = VecX::LinSpaced(20, 0.0, 1.0);
    VecX t = VecX::LinSpaced(20, 0.0, 1.0);

    vector<double> px;
    px.reserve(s.size());
    vector<double> py;
    py.reserve(s.size());
    vector<double> pz;
    pz.reserve(s.size());

    for (const auto x: s)
    {
        auto value = surface.Evaluate(x, 0.0);
        px.push_back(value.x());
        py.push_back(value.y());
        pz.push_back(value.z());
        value = surface.Evaluate(x, 1.0);
        px.push_back(value.x());
        py.push_back(value.y());
        pz.push_back(value.z());
    }
    for (const auto y: t)
    {
        auto value = surface.Evaluate(0.0, y);
        px.push_back(value.x());
        py.push_back(value.y());
        pz.push_back(value.z());
        value = surface.Evaluate(1.0, y);
        px.push_back(value.x());
        py.push_back(value.y());
        pz.push_back(value.z());
    }
    plot3(px, py, pz, "r.");
    show();
}

void CaseC()
{
    double w = sqrt(2.0) / 2.0;

    ControlPointGrid grid;
    grid.UCount = 9;
    grid.VCount = 2;
    grid.ControlPoints.emplace_back(0.5, 0.0, 0.0, 1);
    grid.ControlPoints.emplace_back(1.0, 0.0, 0.0, w);
    grid.ControlPoints.emplace_back(1.0, 0.5, 0.0, 1);
    grid.ControlPoints.emplace_back(1.0, 1.0, 0.0, w);
    grid.ControlPoints.emplace_back(0.5, 1.0, 0.0, 1);
    grid.ControlPoints.emplace_back(0.0, 1.0, 0.0, w);
    grid.ControlPoints.emplace_back(0.0, 0.5, 0.0, 1);
    grid.ControlPoints.emplace_back(0.0, 0.0, 0.0, w);
    grid.ControlPoints.emplace_back(0.5, 0.0, 0.0, 1);

    grid.ControlPoints.emplace_back(0.5, 0.0, 3.0, 1);
    grid.ControlPoints.emplace_back(1.0, 0.0, 3.0, w);
    grid.ControlPoints.emplace_back(1.0, 0.5, 3.0, 1);
    grid.ControlPoints.emplace_back(1.0, 1.0, 3.0, w);
    grid.ControlPoints.emplace_back(0.5, 1.0, 3.0, 1);
    grid.ControlPoints.emplace_back(0.0, 1.0, 3.0, w);
    grid.ControlPoints.emplace_back(0.0, 0.5, 3.0, 1);
    grid.ControlPoints.emplace_back(0.0, 0.0, 3.0, w);
    grid.ControlPoints.emplace_back(0.5, 0.0, 3.0, 1);

    Surface surface;
    surface.DegreeU = 2;
    surface.DegreeV = 1;
    surface.KnotsU = KnotVector{{0.0, 0.0, 0.0, 1 / 4.0, 1 / 4.0, 1 / 2.0, 1 / 2.0, 3 / 4.0, 3 / 4.0, 1.0, 1.0, 1.0}};
    surface.KnotsV = KnotVector{{0.0, 0.0, 1.0, 1.0}};
    surface.ControlPoints = grid;

    VecX s = VecX::LinSpaced(100, 0.0, 1.0);
    VecX t = VecX::LinSpaced(8, 0.0, 1.0);

    vector<double> px;
    px.reserve(s.size());
    vector<double> py;
    py.reserve(s.size());
    vector<double> pz;
    pz.reserve(s.size());

    for (auto x: s)
    {
        for (auto y: t)
        {
            auto value = surface.Evaluate(x, y);
            px.push_back(value.x());
            py.push_back(value.y());
            pz.push_back(value.z());
        }
    }
    plot3(px, py, pz, "r.");
    show();
}

void CaseD()
{
    double w = sqrt(2.0) / 2.0;

    ControlPointGrid grid;
    grid.UCount = 9;
    grid.VCount = 2;
    grid.ControlPoints.emplace_back(0.5, 0.0, 0.0, 1);
    grid.ControlPoints.emplace_back(1.0, 0.0, 0.0, w);
    grid.ControlPoints.emplace_back(1.0, 0.5, 0.0, 1);
    grid.ControlPoints.emplace_back(1.0, 1.0, 0.0, w);
    grid.ControlPoints.emplace_back(0.5, 1.0, 0.0, 1);
    grid.ControlPoints.emplace_back(0.0, 1.0, 0.0, w);
    grid.ControlPoints.emplace_back(0.0, 0.5, 0.0, 1);
    grid.ControlPoints.emplace_back(0.0, 0.0, 0.0, w);
    grid.ControlPoints.emplace_back(0.5, 0.0, 0.0, 1);

    grid.ControlPoints.emplace_back(0.5, 0.0, 3.0, 1);
    grid.ControlPoints.emplace_back(1.0, 0.0, 3.0, w);
    grid.ControlPoints.emplace_back(1.0, 0.5, 3.0, 1);
    grid.ControlPoints.emplace_back(1.0, 1.0, 3.0, w);
    grid.ControlPoints.emplace_back(0.5, 1.0, 3.0, 1);
    grid.ControlPoints.emplace_back(0.0, 1.0, 3.0, w);
    grid.ControlPoints.emplace_back(0.0, 0.5, 3.0, 1);
    grid.ControlPoints.emplace_back(0.0, 0.0, 3.0, w);
    grid.ControlPoints.emplace_back(0.5, 0.0, 3.0, 1);

    Surface surface;
    surface.DegreeU = 2;
    surface.DegreeV = 1;
    surface.KnotsU = KnotVector{{0.0, 0.0, 0.0, 1 / 4.0, 1 / 4.0, 1 / 2.0, 1 / 2.0, 3 / 4.0, 3 / 4.0, 1.0, 1.0, 1.0}};
    surface.KnotsV = KnotVector{{0.0, 0.0, 1.0, 1.0}};
    surface.ControlPoints = grid;

    VecX s = VecX::LinSpaced(100, 0.0, 1.0);
    VecX t = VecX::LinSpaced(30, 0.0, 1.0);

    vector<double> px;
    px.reserve(s.size());
    vector<double> py;
    py.reserve(s.size());
    vector<double> pz;
    pz.reserve(s.size());

    for (const auto x: s)
    {
        auto value = surface.Evaluate(x, 0.0);
        px.push_back(value.x());
        py.push_back(value.y());
        pz.push_back(value.z());
        value = surface.Evaluate(x, 1.0);
        px.push_back(value.x());
        py.push_back(value.y());
        pz.push_back(value.z());
    }
    plot3(px, py, pz, "r.");
    hold(on);
    px.clear();
    py.clear();
    pz.clear();
    for (const auto y: t)
    {
        auto value = surface.Evaluate(0.0, y);
        px.push_back(value.x());
        py.push_back(value.y());
        pz.push_back(value.z());
        value = surface.Evaluate(1.0, y);
        px.push_back(value.x() + 0.1);
        py.push_back(value.y() + 0.1);
        pz.push_back(value.z());
    }
    plot3(px, py, pz, "b.");
    show();
}

int main()
{
    //CaseA();
    //CaseB();
    //CaseC();
    CaseD();

    return 0;
}

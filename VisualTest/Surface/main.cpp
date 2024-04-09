#include <vector>
#include <libnurbs/Core/Grid.hpp>
#include <libnurbs/Core/KnotVector.hpp>
#include <libnurbs/Surface/Surface.hpp>
#include <matplot/matplot.h>

using namespace std;
using namespace libnurbs;
using namespace matplot;

void ShowSurface(const Surface& surface, double step = 0.1)
{
    using namespace matplot;
    auto [xi, eta] = meshgrid(iota(0, step, 1));
    size_t len = xi.size();
    std::vector X(len, std::vector<double>(len));
    std::vector Y(len, std::vector<double>(len));
    std::vector Z(len, std::vector<double>(len));

    for (size_t i = 0; i < len; ++i)
    {
        for (size_t j = 0; j < len; ++j)
        {
            auto eval_result = surface.Evaluate(xi[i][j], eta[i][j]);
            X[i][j] = eval_result.x();
            Y[i][j] = eval_result.y();
            Z[i][j] = eval_result.z();
        }
    }
    surf(X, Y, Z)->face_alpha(0.5).edge_color("none");
    show();
}


Surface CaseSemiCylinder()
{
    ControlPointGrid grid;
    grid.UCount = 3;
    grid.VCount = 2;
    grid.Values.emplace_back(1, 0, 0, 1);
    grid.Values.emplace_back(1, 1, 0, 1);
    grid.Values.emplace_back(0, 1, 0, 2);
    grid.Values.emplace_back(1, 0, 1, 1);
    grid.Values.emplace_back(1, 1, 1, 1);
    grid.Values.emplace_back(0, 1, 1, 2);

    Surface surface;
    surface.DegreeU = 2;
    surface.DegreeV = 1;
    surface.KnotsU = KnotVector{{0.0, 0.0, 0.0, 1.0, 1.0, 1.0}};
    surface.KnotsV = KnotVector{{0.0, 0.0, 1.0, 1.0}};
    surface.ControlPoints = grid;

    return surface;
}


Surface CaseCylinder()
{
    double w = sqrt(2.0) / 2.0;

    ControlPointGrid grid;
    grid.UCount = 9;
    grid.VCount = 2;
    grid.Values.emplace_back(0.5, 0.0, 0.0, 1);
    grid.Values.emplace_back(1.0, 0.0, 0.0, w);
    grid.Values.emplace_back(1.0, 0.5, 0.0, 1);
    grid.Values.emplace_back(1.0, 1.0, 0.0, w);
    grid.Values.emplace_back(0.5, 1.0, 0.0, 1);
    grid.Values.emplace_back(0.0, 1.0, 0.0, w);
    grid.Values.emplace_back(0.0, 0.5, 0.0, 1);
    grid.Values.emplace_back(0.0, 0.0, 0.0, w);
    grid.Values.emplace_back(0.5, 0.0, 0.0, 1);

    grid.Values.emplace_back(0.5, 0.0, 3.0, 1);
    grid.Values.emplace_back(1.0, 0.0, 3.0, w);
    grid.Values.emplace_back(1.0, 0.5, 3.0, 1);
    grid.Values.emplace_back(1.0, 1.0, 3.0, w);
    grid.Values.emplace_back(0.5, 1.0, 3.0, 1);
    grid.Values.emplace_back(0.0, 1.0, 3.0, w);
    grid.Values.emplace_back(0.0, 0.5, 3.0, 1);
    grid.Values.emplace_back(0.0, 0.0, 3.0, w);
    grid.Values.emplace_back(0.5, 0.0, 3.0, 1);

    Surface surface;
    surface.DegreeU = 2;
    surface.DegreeV = 1;
    surface.KnotsU = KnotVector{{0.0, 0.0, 0.0, 1 / 4.0, 1 / 4.0, 1 / 2.0, 1 / 2.0, 3 / 4.0, 3 / 4.0, 1.0, 1.0, 1.0}};
    surface.KnotsV = KnotVector{{0.0, 0.0, 1.0, 1.0}};
    surface.ControlPoints = grid;

    return surface;
}



Surface CaseSaddle()
{
    ControlPointGrid grid;
    grid.UCount = 3;
    grid.VCount = 3;
    grid.Values.emplace_back(0, 1000, 0, 1);
    grid.Values.emplace_back(500, 1000, -300, 1);
    grid.Values.emplace_back(1000, 1000, 0, 1);
    grid.Values.emplace_back(0, 500, 0, 1);
    grid.Values.emplace_back(500, 500, 0, 1);
    grid.Values.emplace_back(1000, 500, 0, 1);
    grid.Values.emplace_back(0, 0, 0, 1);
    grid.Values.emplace_back(500, 0, 300, 1);
    grid.Values.emplace_back(1000, 0, 0, 1);

    Surface surface;
    surface.DegreeU = 2;
    surface.DegreeV = 2;
    surface.KnotsU = KnotVector{{0, 0, 0, 1, 1, 1}};
    surface.KnotsV = KnotVector{{0, 0, 0, 1, 1, 1}};
    surface.ControlPoints = grid;

    return surface;
}

int main()
{
    //ShowSurface(CaseSemiCylinder());
    ShowSurface(CaseCylinder(), 0.1);
    //ShowSurface(CaseSaddle());

    return 0;
}

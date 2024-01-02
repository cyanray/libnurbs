#include "libnurbs/Geometry/GeomRect.hpp"

#include "libnurbs/Surface/Surface.hpp"

namespace libnurbs
{
    Surface GeomRect::GetSurface() const
    {
        assert(DegreeU >= 1);
        assert(DegreeV >= 1);
        assert(ControlPointCountU >= DegreeU + 1);
        assert(ControlPointCountV >= DegreeV + 1);
        const int knot_count_u = ControlPointCountU + DegreeU + 1;
        const int knot_count_v = ControlPointCountV + DegreeV + 1;
        Surface surface;
        surface.DegreeU = DegreeU;
        surface.DegreeV = DegreeV;
        surface.KnotsU = KnotVector::Uniform(DegreeU, knot_count_u);
        surface.KnotsV = KnotVector::Uniform(DegreeV, knot_count_v);
        surface.ControlPoints = {ControlPointCountU, ControlPointCountV};
        Vec3 k1 = (RightBottomPoint - LeftBottomPoint) / (ControlPointCountU - 1);
        Vec3 k2 = (LeftTopPoint - LeftBottomPoint) / (ControlPointCountV - 1);
        for (int i = 0; i < ControlPointCountU; ++i)
        {
            for (int j = 0; j < ControlPointCountV; ++j)
            {
                surface.ControlPoints.Get(i, j).head<3>() = LeftBottomPoint + i * k1 + j * k2;
            }
        }
        return surface;
    }
}

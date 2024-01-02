#include "libnurbs/Geometry/GeomSegment.hpp"

#include "libnurbs/Curve/Curve.hpp"

namespace libnurbs
{
    Curve GeomSegment::GetCurve() const
    {
        assert(Degree >= 1);
        assert(ControlPointCount >= Degree + 1);
        const int knot_count = ControlPointCount + Degree + 1;
        Curve curve;
        curve.Degree = Degree;
        curve.Knots = KnotVector::Uniform(Degree, knot_count);
        curve.ControlPoints.resize(ControlPointCount);
        const Vec3 k = (RightPoint - LeftPoint) / (ControlPointCount - 1);
        for (int i = 0; i < ControlPointCount; ++i)
        {
            curve.ControlPoints[i].head<3>() = LeftPoint + i * k;
        }
        return curve;
    }
}

#include "libnurbs/Geometry/GeomSegment.hpp"

#include "libnurbs/Curve/Curve.hpp"

namespace libnurbs
{
    Curve GeomSegment::GetCurve() const
    {
        assert(ControlPointCount >= 2);
        assert(Degree >= 1);
        int knot_count = ControlPointCount + Degree + 1;
        Curve curve;
        curve.Degree = Degree;
        auto& knots = curve.Knots.Values();
        knots.resize(knot_count);
        for (int i = 0; i < Degree + 1; ++i)
        {
            knots[i] = 0.0;
        }
        Numeric t = 1.0 / (knot_count - (Degree + 1) * 2);
        for (int i = Degree + 1; i < knot_count - (Degree + 1); ++i)
        {
            knots[i] = i * t;
        }
        for (int i = knot_count - (Degree + 1); i < knot_count; ++i)
        {
            knots[i] = 1.0;
        }
        curve.ControlPoints.resize(ControlPointCount);
        const Vec3 k = (RightPoint - LeftPoint) / (ControlPointCount - 1);
        for (int i = 0; i < ControlPointCount; ++i)
        {
            curve.ControlPoints[i] = LeftPoint + i * k;
        }
        return curve;
    }
}

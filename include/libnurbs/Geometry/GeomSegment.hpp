#pragma once

#include "libnurbs/Core/Typedefs.hpp"

namespace libnurbs
{
    class Curve;
}

namespace libnurbs
{
    class GeomSegment
    {
    public:
        int Degree{1};
        int ControlPointCount{2};
        Vec3 LeftPoint = Vec3::Zero();
        Vec3 RightPoint = Vec3::Zero();

    public:
        static GeomSegment Make(const Vec3& left, const Vec3& right)
        {
            GeomSegment result;
            result.LeftPoint = left;
            result.RightPoint = right;
            return result;
        }

        [[nodiscard]] Curve GetCurve() const;
    };
}

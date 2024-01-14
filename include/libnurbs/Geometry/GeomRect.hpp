#pragma once

#include "libnurbs/Core/Typedefs.hpp"

namespace libnurbs
{
    class Surface;
}

namespace libnurbs
{
    class GeomRect
    {
    public:
        int DegreeU{1};
        int DegreeV{1};
        int ControlPointCountU{2};
        int ControlPointCountV{2};
        Vec3 LeftBottomPoint = Vec3::Zero();
        Vec3 RightBottomPoint = Vec3::Zero();
        Vec3 LeftTopPoint = Vec3::Zero();
        Vec3 RightTopPoint = Vec3::Zero();

    public:
        static GeomRect Make(const Vec3& left_bottom, const Vec3& right_bottom, const Vec3& left_top, const Vec3& right_top)
        {
            GeomRect result;
            result.LeftBottomPoint = left_bottom;
            result.RightBottomPoint = right_bottom;
            result.LeftTopPoint = left_top;
            result.RightTopPoint = right_top;
            return result;
        }

        [[nodiscard]] Surface GetSurface() const;
    };
}

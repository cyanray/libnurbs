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
        [[nodiscard]] Surface GetSurface() const;
    };
}

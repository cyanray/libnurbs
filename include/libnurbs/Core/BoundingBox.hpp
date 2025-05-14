#pragma once
#include "Typedefs.hpp"

namespace libnurbs
{
    class BoundingBox
    {
    public:
        Vec3 Min{0.0, 0.0, 0.0};
        Vec3 Max{0.0, 0.0, 0.0};

    public:
        BoundingBox() = default;

        BoundingBox(const Vec3& min, const Vec3& max)
        {
            Min = min.cwiseMin(max);
            Max = min.cwiseMax(max);
        }

        ~BoundingBox() = default;

        void ExpandToInclude(const Vec3& point)
        {
            Min = Min.cwiseMin(point);
            Max = Min.cwiseMax(point);
        }

        void ExpandToInclude(const BoundingBox& box)
        {
            Min = Min.cwiseMin(box.Min);
            Max = Max.cwiseMax(box.Max);
        }

        bool Contains(const Vec3& point) const
        {
            return (point.array() >= Min.array()).all() &&
                   (point.array() <= Max.array()).all();
        }

        bool Contains(const BoundingBox& box) const
        {
            return Contains(box.Min) && Contains(box.Max);
        }

        Vec3 Center() const
        {
            return (Min + Max) * 0.5;
        }

        Vec3 Length() const
        {
            return Max - Min;
        }

        double Volume() const
        {
            Vec3 length = Length();
            return length.prod();
        }
    };
}

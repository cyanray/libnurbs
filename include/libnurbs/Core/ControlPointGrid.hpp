#pragma once

#include <libnurbs/Core/Typedefs.hpp>
#include <vector>

namespace libnurbs
{
    class ControlPointGrid
    {
    public:
        int UCount{INVALID_INDEX};
        int VCount{INVALID_INDEX};
        std::vector<Vec4> ControlPoints{};

    public:
        ControlPointGrid() = default;

        ControlPointGrid(int u_count, int v_count)
            : UCount(u_count),
              VCount(v_count),
              ControlPoints(u_count * v_count)
        {
        }

        /*
        // need c++23 features:
        template<typename Self>
        decltype(auto) operator[](this Self&& self, index_u, int index_v)
        {
            return self.values[index_v * self.UCount + self.index_u];
        }
         */

        Vec4& get(int index_u, int index_v)
        {
            assert(index_u < UCount && index_v < VCount);
            return ControlPoints[index_v * UCount + index_u];
        }

        [[nodiscard]] const Vec4& get(int index_u, int index_v) const
        {
            assert(index_u < UCount && index_v < VCount);
            return ControlPoints[index_v * UCount + index_u];
        }

        [[nodiscard]] int size() const noexcept
        {
            return UCount * VCount;
        }
    };
}

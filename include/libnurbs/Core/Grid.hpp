#pragma once

#include <libnurbs/Core/Typedefs.hpp>
#include <vector>

namespace libnurbs
{
    template<typename T>
    class Grid
    {
    public:
        int UCount{INVALID_INDEX};
        int VCount{INVALID_INDEX};
        std::vector<T> Values{};

    public:
        Grid() = default;

        Grid(int u_count, int v_count)
            : UCount(u_count),
              VCount(v_count),
              Values(u_count * v_count)
        {
        }

        Grid(int u_count, int v_count, T&& init_value)
            : UCount(u_count),
              VCount(v_count),
              Values(u_count * v_count, init_value)
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

        T& Get(int index_u, int index_v)
        {
            assert(index_u < UCount && index_v < VCount);
            return Values[index_v * UCount + index_u];
        }

        [[nodiscard]] const T& Get(int index_u, int index_v) const
        {
            assert(index_u < UCount && index_v < VCount);
            return Values[index_v * UCount + index_u];
        }

        [[nodiscard]] int Count() const noexcept
        {
            return UCount * VCount;
        }
    };
}

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

        [[nodiscard]] std::vector<T> GetU(int index_u) const
        {
            assert(index_u < UCount);
            std::vector<T> result(VCount);
            for (int i = 0; i < VCount; ++i)
            {
                result[i] = Values[i * UCount + index_u];
            }
            return result;
        }

        [[nodiscard]] std::vector<T> GetV(int index_v) const
        {
            assert(index_v < VCount);
            std::vector<T> result(UCount);
            for (int i = 0; i < UCount; ++i)
            {
                result[i] = Values[index_v * UCount + i];
            }
            return result;
        }

        void SetU(int index_u, const std::vector<T>& values)
        {
            assert(index_u < UCount);
            assert(values.size() == VCount);
            for (int i = 0; i < VCount; ++i)
            {
                Values[i * UCount + index_u] = values[i];
            }
        }

        void SetV(int index_v, const std::vector<T>& values)
        {
            assert(index_v < VCount);
            assert(values.size() == UCount);
            for (int i = 0; i < UCount; ++i)
            {
                Values[index_v * UCount + i] = values[i];
            }
        }

        [[nodiscard]] int Count() const noexcept
        {
            return UCount * VCount;
        }

        void InsertV(int v_index, const T& default_value = T())
        {
            assert(v_index >= 0 && v_index < VCount);
            auto iter = Values.begin() + v_index * UCount;
            Values.insert(iter, UCount, default_value);
            VCount += 1;
        }

        void InsertU(int u_index, const T& default_value = T())
        {
            assert(u_index >= 0 && u_index < UCount);
            int float_idx = UCount * VCount - 1;
            UCount += 1;
            Values.resize(UCount * VCount);

            for (int it_v = VCount - 1; it_v >= 0; --it_v)
            {
                int start_idx = it_v * UCount;
                int insertion_idx = start_idx + u_index;
                int end_idx = (it_v + 1) * UCount;

                for (int i = end_idx - 1; i > insertion_idx; --i)
                {
                    Values[i] = std::move(Values[float_idx--]);
                }
                for (int i = insertion_idx - 1; i >= start_idx; --i)
                {
                    Values[i] = std::move(Values[float_idx--]);
                }
                Values[insertion_idx] = default_value;
            }
        }
    };
}

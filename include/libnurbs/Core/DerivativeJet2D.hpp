#pragma once

#include <cassert>
#include <vector>

namespace libnurbs
{
    template <typename T>
    class DerivativeJet2D
    {
    public:
        DerivativeJet2D() = default;

        explicit DerivativeJet2D(const int max_order)
        {
            Reset(max_order);
        }

        void Reset(const int max_order)
        {
            assert(max_order >= 0);
            MaxOrderValue = max_order;
            Values.resize(TermCount(max_order));
        }

        [[nodiscard]] int MaxOrder() const noexcept
        {
            return MaxOrderValue;
        }

        [[nodiscard]] int Count() const noexcept
        {
            return static_cast<int>(Values.size());
        }

        [[nodiscard]] bool Contains(const int order_u, const int order_v) const noexcept
        {
            return order_u >= 0 && order_v >= 0 && order_u + order_v <= MaxOrderValue;
        }

        [[nodiscard]] T& Get(const int order_u, const int order_v)
        {
            assert(Contains(order_u, order_v));
            return Values[Index(order_u, order_v)];
        }

        [[nodiscard]] const T& Get(const int order_u, const int order_v) const
        {
            assert(Contains(order_u, order_v));
            return Values[Index(order_u, order_v)];
        }

    private:
        [[nodiscard]] static constexpr int TermCount(const int max_order) noexcept
        {
            return (max_order + 1) * (max_order + 2) / 2;
        }

        [[nodiscard]] static constexpr int Index(const int order_u, const int order_v) noexcept
        {
            const int total_order = order_u + order_v;
            return total_order * (total_order + 1) / 2 + order_v;
        }

    private:
        int MaxOrderValue{-1};
        std::vector<T> Values{};
    };
}

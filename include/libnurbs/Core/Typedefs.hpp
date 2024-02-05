#pragma once

#include <Eigen/Dense>

namespace libnurbs
{
    constexpr int INVALID_DEGREE{-1};

    constexpr int INVALID_VALUE{-1};

    constexpr int INVALID_INDEX{-1};

    using Numeric = double;

    using Vec3 = Eigen::Vector<Numeric, 3>;

    using Vec4 = Eigen::Vector<Numeric, 4>;

    using VecX = Eigen::Vector<Numeric, Eigen::Dynamic>;

    using MatX = Eigen::Matrix<Numeric, Eigen::Dynamic, Eigen::Dynamic>;

    template<typename T>
    class Grid;

    using ControlPointGrid = Grid<Vec4>;

    inline Vec4 ToHomo(const Vec4& vec)
    {
        Vec4 result(vec);
        result.head<3>() *= result.w();
        return result;
    }

    inline Vec4 FromHomo(const Vec4& vec)
    {
        Vec4 result(vec);
        result.head<3>() /= result.w();
        return result;
    }
}
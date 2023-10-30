#pragma once

#include <libnurbs/Core/Typedefs.hpp>
#include <libnurbs/Core/KnotVector.hpp>
#include <libnurbs/Core/ControlPointGrid.hpp>

namespace libnurbs
{
    class Surface
    {
    public:
        int DegreeU{INVALID_DEGREE};
        int DegreeV{INVALID_DEGREE};
        KnotVector KnotsU{};
        KnotVector KnotsV{};
        ControlPointGrid ControlPoints{};
    public:
        Surface() = default;

        bool Validate() const;

        [[nodiscard]] Vec3 Evaluate(Numeric u, Numeric v) const;

    };


}
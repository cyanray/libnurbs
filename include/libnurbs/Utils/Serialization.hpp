#pragma once
#include "libnurbs/Core/KnotVector.hpp"

namespace libnurbs::Utils
{
    void ReadKnotVectorFromStream(KnotVector& knot_vector, std::istream& is);

    void ReadVec4FromStream(Vec4& vec, std::istream& is);

    void ReadVec4FromText(Vec4& vec, const string& text);

    void WriteKnotVectorToStream(const KnotVector& knot_vector, std::ostream& os);

    void WriteVec4ToStream(const Vec4& vec, std::ostream& os);


    void TrimString(string& str);
}

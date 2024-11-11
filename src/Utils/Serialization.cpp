#include "libnurbs/Utils/Serialization.hpp"

using namespace libnurbs;

namespace libnurbs::Utils
{
    void ReadKnotVectorFromStream(KnotVector& knot_vector, std::istream& is)
    {
        int count;
        is.read(reinterpret_cast<char*>(&count), sizeof(int));
        if (is.fail() || count < 0)
        {
            throw std::runtime_error("Failed to read KnotVector count from binary stream.");
        }
        knot_vector.Values().resize(count);
        is.read(reinterpret_cast<char*>(knot_vector.Values().data()), sizeof(Numeric) * count);
        if (is.fail())
        {
            throw std::runtime_error("Failed to read KnotVector data from binary stream.");
        }
    }

    void ReadVec4FromStream(Vec4& vec, std::istream& is)
    {
        is.read(reinterpret_cast<char*>(vec.data()), sizeof(Numeric) * 4);
        if (is.fail())
        {
            throw std::runtime_error("Failed to read Vec4 from binary stream.");
        }
    }

    void ReadVec4FromText(Vec4& vec, const std::string& text)
    {
        std::stringstream ss(text);
        ss >> vec.x() >> vec.y() >> vec.z() >> vec.w();
        if (ss.fail())
        {
            throw std::runtime_error("Failed to parse Vec4 from text: " + text);
        }
    }

    void WriteKnotVectorToStream(const KnotVector& knot_vector, std::ostream& os)
    {
        int count = static_cast<int>(knot_vector.Values().size());
        os.write(reinterpret_cast<const char*>(&count), sizeof(int));
        if (os.fail())
        {
            throw std::runtime_error("Failed to write KnotVector count to binary stream.");
        }

        os.write(reinterpret_cast<const char*>(knot_vector.Values().data()), sizeof(Numeric) * count);
        if (os.fail())
        {
            throw std::runtime_error("Failed to write KnotVector data to binary stream.");
        }
    }

    void WriteVec4ToStream(const Vec4& vec, std::ostream& os)
    {
        os.write(reinterpret_cast<const char*>(vec.data()), sizeof(Numeric) * 4);
        if (os.fail())
        {
            throw std::runtime_error("Failed to write Vec4 to binary stream.");
        }
    }

    void TrimString(std::string& str)
    {
        // Trim from start (left)
        str.erase(str.begin(),
                  std::find_if(str.begin(),
                               str.end(),
                               [](unsigned char ch) { return !std::isspace(ch); }));

        // Trim from end (right)
        str.erase(std::find_if(str.rbegin(),
                               str.rend(),
                               [](unsigned char ch) { return !std::isspace(ch); }).base(),
                  str.end());
    }
}

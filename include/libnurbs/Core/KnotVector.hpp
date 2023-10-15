#pragma once

#include <libnurbs/Core/Typedefs.hpp>
#include <vector>

using std::vector;

namespace libnurbs
{
    class KnotVector
    {
    private:
        vector <Numeric> m_Values{};
    public:
        struct KnotPair;
    public:
        KnotVector() = default;

        explicit KnotVector(const vector <Numeric>& values);

        explicit KnotVector(const vector <KnotPair>& knots);

        [[nodiscard]] vector<KnotPair> GetKnotPairs() const;

        [[nodiscard]] auto& Values() const
        {
            return m_Values;
        }

        [[nodiscard]] bool IsNonperiodic() const;

        [[nodiscard]] bool IsUniform() const;

        /**
         * @brief Find the index of the span contains u into the knot vector.
         * @param degree
         * @param u
         * @return
         */
        [[nodiscard]] int FindSpanIndex(int degree, Numeric u) const;

    };

    struct KnotVector::KnotPair
    {
        Numeric Value{};
        int Multiplicity{};

        KnotPair(Numeric value, int multiplicity) : Value(value), Multiplicity(multiplicity) {}
    };

}
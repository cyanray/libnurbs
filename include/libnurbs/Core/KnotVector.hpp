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

        struct KnotSpan;
    public:
        KnotVector() = default;

        explicit KnotVector(const vector <Numeric>& values);

        explicit KnotVector(const vector <KnotPair>& knots);

        [[nodiscard]] vector <KnotPair> GetKnotPairs() const;

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

        [[nodiscard]] KnotSpan FindSpan(int degree, Numeric u) const;

    };

    struct KnotVector::KnotPair
    {
        Numeric Value{INVALID_INDEX};
        int Multiplicity{INVALID_VALUE};

        KnotPair() = default;

        KnotPair(Numeric value, int multiplicity) : Value(value), Multiplicity(multiplicity) {}
    };

    struct KnotVector::KnotSpan
    {
        KnotPair Left{};
        KnotPair Right{};

        KnotSpan() = default;

        KnotSpan(const KnotPair& left, const KnotPair& right) : Left(left), Right(right) {}
    };

}
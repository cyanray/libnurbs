#pragma once

#include <iostream>
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

        static KnotVector Uniform(int degree, int knots_count);

        [[nodiscard]] vector <KnotPair> GetKnotPairs() const;

        [[nodiscard]] auto& Values() const
        {
            return m_Values;
        }

        [[nodiscard]] auto& Values()
        {
            return m_Values;
        }

        [[nodiscard]] int Count() const
        {
            return (int)m_Values.size();
        }

        Numeric& operator[](int index)
        {
            return m_Values[index];
        }

        const Numeric& operator()(int index) const
        {
            return m_Values[index];
        }

        [[nodiscard]] bool IsValid() const;

        [[nodiscard]] bool IsClamped() const;

        [[nodiscard]] bool IsUniform() const;

        /**
         * @brief Find the index of the span contains u into the knot vector.
         * @param degree
         * @param u
         * @return Span1, [Span1,Span2).
         */
        [[nodiscard]] int FindSpanIndex(int degree, Numeric u) const;

        [[nodiscard]] KnotSpan FindSpan(Numeric u) const;

        [[nodiscard]] int DetectDegree() const;

        [[nodiscard]] int GetMultiplicity(Numeric u) const;

        /**
         * \brief Insert a knot value into the knot vector.
         * \param value
         * \param times
         * \return Span index in original KnotVector that contains the new knot value.
         */
        int InsertKnot(Numeric value, int times = 1);

        void Reverse();
    };

    struct KnotVector::KnotPair
    {
        int Index{INVALID_INDEX};
        Numeric Value{INVALID_VALUE};
        int Multiplicity{INVALID_VALUE};

        KnotPair() = default;

        KnotPair(int index, Numeric value, int multiplicity)
                : Index(index), Value(value), Multiplicity(multiplicity) {}
    };

    struct KnotVector::KnotSpan
    {
        int Index{INVALID_INDEX};
        KnotPair Left{};
        KnotPair Right{};

        KnotSpan() = default;

        KnotSpan(int index, const KnotPair& left, const KnotPair& right)
                : Index(index), Left(left), Right(right) {}
    };

}
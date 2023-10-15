#include <libnurbs/Core/KnotVector.hpp>
#include <stdexcept>
#include <algorithm>

using namespace std;

namespace
{
    using namespace libnurbs;

    bool ValidateKnots(const vector<Numeric>& knots)
    {
        if (knots.empty()) return false;
        Numeric last_value = 0.0;
        for (Numeric val: knots)
        {
            if (val < last_value) return false;
            last_value = val;
        }
        if (last_value != 1.0) return false;
        return true;
    }
}

namespace libnurbs
{

    KnotVector::KnotVector(const vector<KnotPair>& knots)
    {
        for (const auto& knot: knots)
        {
            for (int i = 0; i < knot.Multiplicity; ++i)
            {
                m_Values.emplace_back(knot.Value);
            }
        }
        if (!ValidateKnots(m_Values)) throw invalid_argument("knots");
    }

    KnotVector::KnotVector(const vector<Numeric>& values)
    {
        if (!ValidateKnots(values)) throw invalid_argument("values");
        m_Values = values;
    }

    vector<KnotVector::KnotPair> KnotVector::GetKnotPairs() const
    {
        if (m_Values.empty()) return {};
        vector<KnotPair> result;
        int len = (int) m_Values.size();
        Numeric last_value = m_Values[0];
        int multiplicity{1};
        for (int i = 1; i < len; ++i)
        {
            Numeric val = m_Values[i];
            if (val != last_value)
            {
                result.emplace_back(last_value, multiplicity);
                last_value = val;
                multiplicity = 1;
            }
            else ++multiplicity;
        }
        result.emplace_back(last_value, multiplicity);
        return result;
    }

    bool KnotVector::IsNonperiodic() const
    {
        auto knot_pairs = GetKnotPairs();
        if (knot_pairs.size() < 2) return false;
        if (knot_pairs.front().Multiplicity != knot_pairs.back().Multiplicity) return false;
        return true;
    }

    bool KnotVector::IsUniform() const
    {
        // TODO: implement
        return true;
    }

    int KnotVector::FindSpanIndex(int degree, Numeric u) const
    {
        int index_first_span = degree;
        int index_last_span = static_cast<int>(m_Values.size()) - degree - 2;
        assert(index_last_span >= 0);
        if (u >= m_Values.back()) return index_last_span;
        if (u <= m_Values.front()) return index_first_span;
        // binary search
        auto it = std::lower_bound(m_Values.begin(), m_Values.end(), u);
        int index_result = (int) std::distance(m_Values.begin(), it);
        return (*it == u) ? index_result : index_result - 1;
    }
}
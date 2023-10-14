#include <libnurbs/Core/KnotVector.hpp>
#include <stdexcept>

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
}
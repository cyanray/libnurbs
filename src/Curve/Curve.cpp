#include "libnurbs/Curve/Curve.hpp"

#include <cassert>
#include <fstream>

#include "libnurbs/Algorithm/DegreeAlgo.hpp"
#include "libnurbs/Algorithm/KnotRemoval.hpp"
#include "libnurbs/Algorithm/MathUtils.hpp"
#include "libnurbs/Basis/BSplineBasis.hpp"
#include "libnurbs/Utils/Serialization.hpp"


using namespace std;
using namespace libnurbs;

void Curve::LoadFromFile(const string& filename)
{
    ifstream file(filename);
    if (!file.is_open())
    {
        throw std::runtime_error("Failed to open file: " + filename);
    }
    LoadFromFile(file);
}

void Curve::LoadFromFile(std::istream& is)
{
    // Read header line
    std::string header_line;
    std::getline(is, header_line);
    if (is.fail())
    {
        throw std::runtime_error("Failed to read header line.");
    }
    std::stringstream ss_header(header_line);
    std::string magic, mode, type;
    ss_header >> magic >> mode >> type;

    if (magic != LIBNURBS_MAGIC)
    {
        throw std::runtime_error("Invalid file magic: " + magic);
    }

    bool mode_binary;
    if (mode == "TXT")
    {
        mode_binary = false;
    }
    else if (mode == "BIN")
    {
        mode_binary = true;
    }
    else
    {
        throw std::runtime_error("Unknown serialization mode: " + mode);
    }

    if (type != "CURVE")
    {
        throw std::runtime_error("Object type mismatch. Expected CURVE, got: " + type);
    }

    if (mode_binary)
    {
        is.read(reinterpret_cast<char*>(&Degree), sizeof(int));
        if (is.fail())
        {
            throw std::runtime_error("Failed to read Degree.");
        }

        Utils::ReadKnotVectorFromStream(Knots, is);

        int cp_count;
        is.read(reinterpret_cast<char*>(&cp_count), sizeof(int));
        if (is.fail() || cp_count < 0)
        {
            throw std::runtime_error("Failed to read ControlPoints count.");
        }
        ControlPoints.resize(cp_count);
        for (int i = 0; i < cp_count; ++i)
        {
            Utils::ReadVec4FromStream(ControlPoints[i], is);
        }
        return;
    }

    // TXT mode
    std::string line;
    std::string current_key;
    std::vector<std::string> content_lines;

    auto process_key = [&](const std::string& key, const std::vector<std::string>& lines)
    {
        if (key == "Degree")
        {
            if (lines.size() != 1)
            {
                throw std::runtime_error("Degree should have exactly one value line.");
            }
            std::stringstream ss(lines[0]);
            ss >> Degree;
            if (ss.fail())
            {
                throw std::runtime_error("Failed to parse Degree value.");
            }
        }
        else if (key == "Knots")
        {
            for (const auto& l : lines)
            {
                std::stringstream ss(l);
                std::string token;
                while (std::getline(ss, token, ','))
                {
                    Utils::TrimString(token);
                    if (!token.empty())
                    {
                        Numeric knot;
                        std::stringstream num_stream(token);
                        num_stream >> knot;
                        if (num_stream.fail())
                        {
                            throw std::runtime_error("Failed to parse knot value: " + token);
                        }
                        Knots.Values().push_back(knot);
                    }
                }
            }
        }
        else if (key == "ControlPoints")
        {
            for (const auto& l : lines)
            {
                std::stringstream ss(l);
                std::string cp_str;
                while (std::getline(ss, cp_str, ','))
                {
                    Utils::TrimString(cp_str);
                    if (!cp_str.empty())
                    {
                        Vec4 vec;
                        Utils::ReadVec4FromText(vec, cp_str);
                        ControlPoints.push_back(vec);
                    }
                }
            }
        }
        else
        {
            throw std::runtime_error("Unknown key encountered: " + key);
        }
    };

    while (std::getline(is, line))
    {
        Utils::TrimString(line);
        if (line.empty()) continue;

        if (line.front() == '*' && line.find(':') != std::string::npos)
        {
            // Process previous key
            if (!current_key.empty())
            {
                process_key(current_key, content_lines);
                content_lines.clear();
            }

            // Extract key and possible inline value
            size_t colon_pos = line.find(':');
            current_key = line.substr(1, colon_pos - 1);
            Utils::TrimString(current_key);

            std::string value_part = line.substr(colon_pos + 1);
            Utils::TrimString(value_part);

            if (!value_part.empty())
            {
                content_lines.push_back(value_part);
            }
        }
        else
        {
            if (current_key.empty())
            {
                throw std::runtime_error("Content found before any key.");
            }
            content_lines.push_back(line);
        }
    }

    // Process the last key
    if (!current_key.empty())
    {
        process_key(current_key, content_lines);
    }
}

void Curve::SaveToFile(const std::string& filename, bool binary_mode) const
{
    std::ofstream file;
    if (binary_mode)
    {
        file.open(filename, std::ios::binary);
    }
    else
    {
        file.open(filename);
    }

    if (!file.is_open())
    {
        throw std::runtime_error("Failed to open file for writing: " + filename);
    }

    SaveToFile(file, binary_mode);
}

void Curve::SaveToFile(std::ostream& os, bool binary_mode) const
{
    // Write header
    os << LIBNURBS_MAGIC << " " << (binary_mode ? "BIN" : "TXT") << " CURVE\n";

    if (binary_mode)
    {
        // Write binary data
        os.write(reinterpret_cast<const char*>(&Degree), sizeof(int));
        if (os.fail())
        {
            throw std::runtime_error("Failed to write Degree to binary stream.");
        }

        Utils::WriteKnotVectorToStream(Knots, os);

        int cp_count = static_cast<int>(ControlPoints.size());
        os.write(reinterpret_cast<const char*>(&cp_count), sizeof(int));
        if (os.fail())
        {
            throw std::runtime_error("Failed to write ControlPoints count to binary stream.");
        }

        for (const auto& cp : ControlPoints)
        {
            Utils::WriteVec4ToStream(cp, os);
        }
    }
    else
    {
        // Write textual data
        os << "*Degree: " << Degree << "\n";

        os << "*Knots:\n";
        for (size_t i = 0; i < Knots.Values().size(); ++i)
        {
            os << Knots.Values()[i];
            if (i != Knots.Values().size() - 1)
            {
                os << ", ";
            }
        }
        os << "\n";

        os << "*ControlPoints:\n";
        for (const auto& cp : ControlPoints)
        {
            os << cp.x() << " " << cp.y() << " " << cp.z() << " " << cp.w() << ",\n";
        }
    }

    if (os.fail())
    {
        throw std::runtime_error("Failed to write data to stream.");
    }
}


Vec3 Curve::Evaluate(Numeric x) const
{
    assert(x >= 0 && x <= 1);
    int index_span = Knots.FindSpanIndex(Degree, x);
    VecX basis = BSplineBasis::Evaluate(Degree, Knots.Values(), index_span, x);
    assert(basis.size() == Degree + 1);
    Vec4 result = Vec4::Zero();
    for (int i = 0; i <= Degree; i++)
    {
        auto point = ToHomo(ControlPoints[index_span - Degree + i]);
        result += basis(i) * point;
    }
    return result.head<3>() / result.w();
}


std::vector<Vec4> Curve::HomogeneousDerivative(Numeric x, int order) const
{
    assert(x >= 0 && x <= 1);
    int index_span = Knots.FindSpanIndex(Degree, x);
    MatX basis = BSplineBasis::EvaluateAll(Degree, Knots.Values(), index_span, x, order);
    assert(basis.rows() == order + 1);
    assert(basis.cols() == Degree + 1);
    std::vector<Vec4> result(order + 1);
    for (int k = 0; k <= order; ++k)
    {
        Vec4 tmp = Vec4::Zero();
        for (int i = 0; i <= Degree; i++)
        {
            auto point = ToHomo(ControlPoints[index_span - Degree + i]);
            tmp.noalias() += point * basis(k, i);
        }
        result[k] = tmp;
    }
    return result;
}

Vec3 Curve::EvaluateDerivative(Numeric x, int order) const
{
    return EvaluateAll(x, order)[order];
}

vector<Vec3> Curve::EvaluateAll(Numeric x, int order) const
{
    auto homo_ders = HomogeneousDerivative(x, order);
    vector<Vec3> result(order + 1, Vec3::Zero());

    // Compute rational derivatives
    Numeric Wders0 = homo_ders[0].w();
    for (int k = 0; k <= order; k++)
    {
        Vec3 Aders = homo_ders[k].head<3>();
        for (int i = 1; i <= k; i++)
        {
            Numeric Wders = homo_ders[i].w();
            Aders.noalias() -= Binomial(k, i) * Wders * result[k - i];
        }
        result[k] = (Aders / Wders0);
    }
    return result;
}

bool Curve::IsRational() const
{
    const auto& cps = ControlPoints;
    if (cps.empty()) return false;
    Numeric w = cps.front().w();
    for (size_t i = 1; i < cps.size(); ++i)
    {
        if (cps[i].w() != w) return true;
    }
    return false;
}

Numeric Curve::SearchParameter(const Vec3& point, Numeric init, Numeric epsilon, Numeric max_iteration_count) const
{
    auto Ri = [&point, this](Numeric u) -> Vec3 { return Evaluate(u) - point; };

    auto fi = [this, &Ri](Numeric u) -> Numeric
    {
        auto ri = Ri(u);
        auto Cu = EvaluateDerivative(u, 1);
        return ri.dot(Cu);
    };

    Numeric u_last = init;
    Numeric current_residual = std::numeric_limits<Numeric>::max();
    Numeric Hk = 1.0;

    auto line_search = [&](Numeric gk) -> Numeric
    {
        Numeric alpha = 1.0;
        Numeric c1 = 1e-4;
        Numeric beta = 0.9;
        int max_line_search_iterations = 20;
        int ls_count = 0;

        while (ls_count++ < max_line_search_iterations)
        {
            Numeric u_trial = u_last - alpha * Hk * gk;
            u_trial = std::clamp(u_trial, Numeric(0), Numeric(1));

            if (Ri(u_trial).norm() <= current_residual - c1 * alpha * gk * Hk * gk)
            {
                return u_trial;
            }
            alpha *= beta;
        }
        Numeric u_final = u_last - alpha * Hk * gk;
        return std::clamp(u_final, Numeric(0), Numeric(1));
    };

    int count = 0;
    while (!((current_residual < epsilon) || (count++ >= max_iteration_count)))
    {
        Numeric gk = fi(u_last);
        Numeric u_new = line_search(gk);

        Numeric sk = u_new - u_last;
        Numeric gk_new = fi(u_new);
        Numeric yk = gk_new - gk;

        Numeric yk_dot_sk = yk * sk;
        if (std::abs(yk_dot_sk) < 1e-16)
        {
            Hk = 1.0;
            yk_dot_sk = 1e-10;
        }

        Numeric rho = 1.0 / yk_dot_sk;
        Hk = (1.0 - rho * yk * sk) * Hk * (1.0 - rho * sk * yk) + rho * sk * sk;

        u_last = u_new;
        current_residual = Ri(u_last).norm();
    }
    return u_last;
}

Numeric Curve::BinarySearchParameter(const Vec3& point, Numeric epsilon, int max_iteration_count) const
{
    Numeric low = 0.0;
    Numeric high = 1.0;

    auto Ri = [&point, this](Numeric u) -> Numeric
    {
        return (Evaluate(u) - point).norm();
    };

    int count = 0;
    while (count++ < max_iteration_count)
    {
        Numeric mid = (low + high) / 2.0;

        Numeric mid_residual = Ri(mid);

        Numeric left_mid = mid - epsilon;
        Numeric right_mid = mid + epsilon;

        if (left_mid < low) left_mid = low;
        if (right_mid > high) right_mid = high;

        Numeric left_residual = Ri(left_mid);
        Numeric right_residual = Ri(right_mid);

        if (left_residual < mid_residual)
        {
            high = mid;
        }
        else if (right_residual < mid_residual)
        {
            low = mid;
        }
        else
        {
            if (std::abs(high - low) < epsilon || mid_residual < epsilon) return mid;
            low = left_mid;
            high = right_mid;
        }
    }

    return (low + high) / 2.0;
}

Curve Curve::InsertKnot(Numeric knot_value) const
{
    Curve result{*this};
    int k = result.Knots.InsertKnot(knot_value);
    result.ControlPoints.insert(result.ControlPoints.begin() + k, Vec4::Zero());
    const auto& knots = this->Knots.Values();
    for (int i = k - Degree + 1; i <= k; i++)
    {
        Numeric alpha = (knot_value - knots[i]) / (knots[i + Degree] - knots[i]);
        auto point_i = ToHomo(this->ControlPoints[i]);
        auto point_im1 = ToHomo(this->ControlPoints[i - 1]);
        result.ControlPoints[i] = FromHomo(alpha * point_i + (1 - alpha) * point_im1);
    }
    return result;
}

Curve Curve::InsertKnot(Numeric knot_value, int times) const
{
    // TODO: better implementation
    std::vector list(times, knot_value);
    return InsertKnot(list);
}

Curve Curve::InsertKnot(std::span<Numeric> knots_to_insert) const
{
    Curve result = *this;
    auto& knots = result.Knots.Values();
    knots.reserve(knots.size() + knots_to_insert.size());
    auto& points = result.ControlPoints;
    points.reserve(points.size() + knots_to_insert.size());
    for (auto knot_value : knots_to_insert)
    {
        int k = result.Knots.FindSpanIndex(Degree, knot_value);
        points.insert(points.begin() + k, points[k]);
        Vec4 point_last = ToHomo(points[k - Degree]);
        for (int i = k - Degree + 1; i <= k; i++)
        {
            Numeric alpha = (knot_value - knots[i]) / (knots[i + Degree] - knots[i]);
            auto point_i = ToHomo(points[i]);
            Vec4 point_new = FromHomo(alpha * point_i + (1 - alpha) * point_last);
            point_last = point_i;
            points[i] = point_new;
        }
        knots.insert(knots.begin() + k + 1, knot_value);
    }
    return result;
}

std::tuple<Curve, int> Curve::RemoveKnot(Numeric knot_remove, int times, Numeric tolerance) const
{
    assert(knot_remove > 0 && knot_remove < 1);
    assert(times > 0);

    Curve result{*this};
    int degree = result.Degree;
    auto& points = result.ControlPoints;
    int t = KnotRemoval(result.Knots, points, degree, knot_remove, times, tolerance);
    return {result, t};
}

Curve Curve::ElevateDegree(int times) const
{
    assert(times >= 1);

    Curve result{*this};
    NurbsDegreeElevation(result.Knots, result.ControlPoints, result.Degree, times);
    return result;
}

Curve Curve::Transform(const Mat3x3& R) const
{
    Curve transformed_curve = *this;
    for (auto& point : transformed_curve.ControlPoints)
    {
        point.head<3>() = R * point.head<3>();
    }
    return transformed_curve;
}

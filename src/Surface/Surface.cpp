#include "libnurbs/Surface/Surface.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>

#include "libnurbs/Algorithm/DegreeAlgo.hpp"
#include "libnurbs/Algorithm/KnotRemoval.hpp"
#include "libnurbs/Algorithm/MathUtils.hpp"
#include "libnurbs/Basis/BSplineBasis.hpp"
#include "libnurbs/Utils/Serialization.hpp"

using namespace std;
using namespace libnurbs;

Surface& Surface::LoadFromFile(const std::string& filename)
{
    std::ifstream file(filename, std::ios::binary);
    if (!file.is_open())
    {
        throw std::runtime_error("Failed to open file: " + filename);
    }
    return LoadFromFile(file);
}

Surface& Surface::LoadFromFile(std::istream& is)
{
    // Read header line
    std::string header_line;
    std::getline(is, header_line);
    if (is.fail())
    {
        throw std::runtime_error("Failed to read header line from stream.");
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

    if (type != "SURFACE")
    {
        throw std::runtime_error("Object type mismatch. Expected SURFACE, got: " + type);
    }

    if (mode_binary)
    {
        // Binary mode
        // Read degrees
        is.read(reinterpret_cast<char*>(&DegreeU), sizeof(int));
        is.read(reinterpret_cast<char*>(&DegreeV), sizeof(int));
        if (is.fail())
        {
            throw std::runtime_error("Failed to read Degrees from binary stream.");
        }

        // Read knot vectors
        Utils::ReadKnotVectorFromStream(KnotsU, is);
        Utils::ReadKnotVectorFromStream(KnotsV, is);

        // Read control point grid dimensions
        is.read(reinterpret_cast<char*>(&ControlPoints.UCount), sizeof(int));
        is.read(reinterpret_cast<char*>(&ControlPoints.VCount), sizeof(int));
        if (is.fail() || ControlPoints.UCount <= 0 || ControlPoints.VCount <= 0)
        {
            throw std::runtime_error("Failed to read control point grid dimensions.");
        }

        // Read control points
        int cp_count = ControlPoints.UCount * ControlPoints.VCount;
        ControlPoints.Values.resize(cp_count);
        for (int i = 0; i < cp_count; ++i)
        {
            Utils::ReadVec4FromStream(ControlPoints.Values[i], is);
        }
    }
    else
    {
        // Text mode
        std::string line;
        std::string current_key;
        std::vector<std::string> content_lines;

        auto process_key = [&](const std::string& key, const std::vector<std::string>& lines)
        {
            if (key == "DegreeU")
            {
                if (lines.size() != 1)
                {
                    throw std::runtime_error("DegreeU should have exactly one value line.");
                }
                std::stringstream ss(lines[0]);
                ss >> DegreeU;
                if (ss.fail())
                {
                    throw std::runtime_error("Failed to parse DegreeU value.");
                }
            }
            else if (key == "DegreeV")
            {
                if (lines.size() != 1)
                {
                    throw std::runtime_error("DegreeV should have exactly one value line.");
                }
                std::stringstream ss(lines[0]);
                ss >> DegreeV;
                if (ss.fail())
                {
                    throw std::runtime_error("Failed to parse DegreeV value.");
                }
            }
            else if (key == "KnotsU")
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
                                throw std::runtime_error("Failed to parse KnotsU value: " + token);
                            }
                            KnotsU.Values().push_back(knot);
                        }
                    }
                }
            }
            else if (key == "KnotsV")
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
                                throw std::runtime_error("Failed to parse KnotsV value: " + token);
                            }
                            KnotsV.Values().push_back(knot);
                        }
                    }
                }
            }
            else if (key == "ControlPoints")
            {
                // Expecting first line to contain UCount and VCount
                if (lines.empty())
                {
                    throw std::runtime_error("ControlPoints section missing grid dimensions.");
                }
                std::stringstream ss_dims(lines[0]);
                ss_dims >> ControlPoints.UCount >> ControlPoints.VCount;
                if (ss_dims.fail() || ControlPoints.UCount <= 0 || ControlPoints.VCount <= 0)
                {
                    throw std::runtime_error("Invalid control point grid dimensions.");
                }

                // Read control points
                int expected_cp_count = ControlPoints.UCount * ControlPoints.VCount;
                if (static_cast<int>(lines.size() - 1) != expected_cp_count)
                {
                    throw std::runtime_error("ControlPoints count does not match grid dimensions.");
                }

                ControlPoints.Values.resize(expected_cp_count);
                for (int i = 1; i < static_cast<int>(lines.size()); ++i)
                {
                    std::string cp_line = lines[i];
                    Utils::TrimString(cp_line);
                    if (!cp_line.empty())
                    {
                        Vec4 vec;
                        Utils::ReadVec4FromText(vec, cp_line);
                        ControlPoints.Values[i - 1] = vec;
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
    return *this;
}

void Surface::SaveToFile(const std::string& filename, bool binary_mode) const
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


void Surface::SaveToFile(std::ostream& os, bool binary_mode) const
{
    // Write header
    os << LIBNURBS_MAGIC << " " << (binary_mode ? "BIN" : "TXT") << " SURFACE\n";

    if (binary_mode)
    {
        // Write binary data

        // Write degrees
        os.write(reinterpret_cast<const char*>(&DegreeU), sizeof(int));
        os.write(reinterpret_cast<const char*>(&DegreeV), sizeof(int));
        if (os.fail())
        {
            throw std::runtime_error("Failed to write Degrees to binary stream.");
        }

        // Write knot vectors
        Utils::WriteKnotVectorToStream(KnotsU, os);
        Utils::WriteKnotVectorToStream(KnotsV, os);

        // Write control point grid dimensions
        os.write(reinterpret_cast<const char*>(&ControlPoints.UCount), sizeof(int));
        os.write(reinterpret_cast<const char*>(&ControlPoints.VCount), sizeof(int));
        if (os.fail())
        {
            throw std::runtime_error("Failed to write control point grid dimensions.");
        }

        // Write control points
        int cp_count = ControlPoints.UCount * ControlPoints.VCount;
        for (int i = 0; i < cp_count; ++i)
        {
            Utils::WriteVec4ToStream(ControlPoints.Values[i], os);
        }
    }
    else
    {
        // Text mode

        // Write degrees
        os << "*DegreeU: " << DegreeU << "\n";
        os << "*DegreeV: " << DegreeV << "\n";

        os << std::fixed << std::setprecision(std::numeric_limits<double>::digits10 + 2);

        // Write knot vectors
        os << "*KnotsU:\n";
        for (size_t i = 0; i < KnotsU.Values().size(); ++i)
        {
            os << KnotsU.Values()[i];
            if (i != KnotsU.Values().size() - 1)
            {
                os << ", ";
            }
        }
        os << "\n";

        os << "*KnotsV:\n";
        for (size_t i = 0; i < KnotsV.Values().size(); ++i)
        {
            os << KnotsV.Values()[i];
            if (i != KnotsV.Values().size() - 1)
            {
                os << ", ";
            }
        }
        os << "\n";

        // Write control points
        os << "*ControlPoints:\n";
        os << ControlPoints.UCount << " " << ControlPoints.VCount << "\n";

        for (int v = 0; v < ControlPoints.VCount; ++v)
        {
            for (int u = 0; u < ControlPoints.UCount; ++u)
            {
                const auto& cp = ControlPoints.Get(u, v);
                os << cp.x() << " " << cp.y() << " " << cp.z() << " " << cp.w() << "\n";
            }
        }
    }

    if (os.fail())
    {
        throw std::runtime_error("Failed to write data to stream.");
    }
}


Vec3 Surface::Evaluate(Numeric u, Numeric v) const
{
    assert(u >= 0 && u <= 1);
    assert(v >= 0 && v <= 1);
    int index_span_u = KnotsU.FindSpanIndex(DegreeU, u);
    int index_span_v = KnotsV.FindSpanIndex(DegreeV, v);
    VecX basis_u = BSplineBasis::Evaluate(DegreeU, KnotsU.Values(), index_span_u, u);
    VecX basis_v = BSplineBasis::Evaluate(DegreeV, KnotsV.Values(), index_span_v, v);
    assert(basis_u.size() == DegreeU + 1);
    assert(basis_v.size() == DegreeV + 1);

    Vec4 result = Vec4::Zero();
    int index_pre_u = index_span_u - DegreeU;
    int index_pre_v = index_span_v - DegreeV;

    for (int j = 0; j <= DegreeV; j++)
    {
        int index_v = index_pre_v + j;
        Vec4 tmp = Vec4::Zero();
        for (int i = 0; i <= DegreeU; i++)
        {
            int index_u = index_pre_u + i;
            auto point = ToHomo(ControlPoints.Get(index_u, index_v));
            tmp.noalias() += basis_u(i) * point;
        }
        result.noalias() += basis_v(j) * tmp;
    }
    return result.head<3>() / result.w();
}

Vec3 Surface::EvaluateDerivative(Numeric u, Numeric v, int order_u, int order_v) const
{
    return EvaluateAll(u, v, order_u, order_v).Get(order_u, order_v);
}

Grid<Vec3> Surface::EvaluateAll(Numeric u, Numeric v, int order_u, int order_v) const
{
    auto homo_ders = HomogeneousDerivative(u, v, order_u, order_v);
    Grid<Vec3> result(order_u + 1, order_v + 1, Vec3::Zero());

    Numeric Wders00 = homo_ders.Get(0, 0).w();

    for (int ou = 0; ou <= order_u; ++ou)
    {
        for (int ov = 0; ov <= order_v; ++ov)
        {
            Vec3 Aders = homo_ders.Get(ou, ov).head<3>();
            for (int i = 1; i <= ou; ++i)
            {
                Numeric Wders = homo_ders.Get(i, 0).w();
                Aders.noalias() -= Binomial(ou, i) * Wders * result.Get(ou - i, ov);
            }

            for (int j = 1; j <= ov; ++j)
            {
                Numeric Wders = homo_ders.Get(0, j).w();
                Aders.noalias() -= Binomial(ov, j) * Wders * result.Get(ou, ov - j);
            }

            for (int i = 1; i <= ou; ++i)
            {
                const int bi = Binomial(ou, i);
                for (int j = 1; j <= ov; ++j)
                {
                    Numeric Wders = homo_ders.Get(i, j).w();
                    Aders.noalias() -= bi * Binomial(ov, j) * Wders * result.Get(ou - i, ov - j);
                }
            }
            result.Get(ou, ov) = Aders / Wders00;
        }
    }

    return result;
}

bool Surface::IsRational() const
{
    const auto& cps = ControlPoints.Values;
    if (cps.empty()) return false;
    Numeric w = cps.front().w();
    for (size_t i = 1; i < cps.size(); ++i)
    {
        if (cps[i].w() != w) return true;
    }
    return false;
}

auto Surface::SearchParameter(const Vec3& point, Numeric init_u, Numeric init_v,
                              Numeric epsilon, Numeric max_iteration_count) const
    -> std::pair<Numeric, Numeric>
{
    assert(init_u >= 0 && init_u <= 1);
    assert(init_v >= 0 && init_v <= 1);
    using Vec2 = Eigen::Vector2<Numeric>;
    using Mat2x2 = Eigen::Matrix<Numeric, 2, 2>;

    auto Ri = [&point, this](Numeric u, Numeric v) -> Vec3
    {
        return Evaluate(u, v) - point;
    };

    auto Ki = [this, &Ri](Numeric u, Numeric v) -> Vec2
    {
        auto ri = Ri(u, v);
        auto Su = EvaluateDerivative(u, v, 1, 0);
        auto Sv = EvaluateDerivative(u, v, 0, 1);
        return Vec2{ri.dot(Su), ri.dot(Sv)};
    };

    Numeric u_last = init_u, v_last = init_v;
    Numeric current_residual = std::numeric_limits<Numeric>::max();
    Mat2x2 Hk = Mat2x2::Identity();

    auto line_search = [&](const Vec2& gk) -> Vec2
    {
        Numeric alpha = 1.0;
        Numeric c1 = 1e-4;
        Numeric beta = 0.9;
        int max_line_search_iterations = 20;
        int ls_count = 0;

        while (ls_count++ < max_line_search_iterations)
        {
            Vec2 uv_trial = Vec2{u_last, v_last} - alpha * Hk * gk;
            uv_trial.x() = std::clamp(uv_trial.x(), Numeric(0), Numeric(1));
            uv_trial.y() = std::clamp(uv_trial.y(), Numeric(0), Numeric(1));

            if (Ri(uv_trial.x(), uv_trial.y()).norm() <= current_residual - c1 * alpha * gk.dot(Hk * gk))
            {
                return uv_trial;
            }
            alpha *= beta;
        }
        Vec2 uv_last = Vec2{u_last, v_last} - alpha * Hk * gk;
        uv_last.x() = std::clamp(uv_last.x(), Numeric(0), Numeric(1));
        uv_last.y() = std::clamp(uv_last.y(), Numeric(0), Numeric(1));
        return uv_last;
    };

    int count = 0;
    while (!((current_residual < epsilon) || (count++ >= max_iteration_count)))
    {
        auto gk = Ki(u_last, v_last);
        Vec2 uv_last = line_search(gk);

        Vec2 sk = uv_last - Vec2{u_last, v_last};
        Vec2 gk_new = Ki(uv_last.x(), uv_last.y());
        Vec2 yk = gk_new - gk;

        Numeric yk_dot_sk = yk.dot(sk);
        if (std::abs(yk_dot_sk) < 1e-16)
        {
            Hk = Mat2x2::Identity();
            yk_dot_sk = 1e-10;
        }

        Numeric rho = 1.0 / yk_dot_sk;
        Mat2x2 I = Mat2x2::Identity();
        Hk = (I - rho * sk * yk.transpose()) * Hk * (I - rho * yk * sk.transpose()) + rho * sk * sk.transpose();

        u_last = uv_last.x();
        v_last = uv_last.y();
        current_residual = Ri(u_last, v_last).norm();
    }
    return {u_last, v_last};
}

auto Surface::SearchParameterOn(const Vec3& point, int direction, Numeric constant,
                                Numeric init_value,
                                Numeric epsilon, Numeric max_iteration_count) const
    -> std::pair<Numeric, Numeric>
{
    auto Ri = [&point, this, direction, constant](Numeric val) -> Vec3
    {
        Numeric u = direction == 0 ? constant : val;
        Numeric v = direction == 1 ? constant : val;
        return Evaluate(u, v) - point;
    };

    auto fi = [this, &Ri, direction, constant](Numeric val) -> Numeric
    {
        Numeric u = direction == 0 ? constant : val;
        Numeric v = direction == 1 ? constant : val;
        auto ri = Ri(val);
        auto Cu = direction == 1
                      ? EvaluateDerivative(u, v, 1, 0)
                      : EvaluateDerivative(u, v, 0, 1);
        return ri.dot(Cu);
    };

    Numeric val_last = init_value;
    Numeric current_residual = std::numeric_limits<Numeric>::max();
    Numeric Hk = 1.0;

    auto line_search = [&](Numeric gk) -> Numeric
    {
        Numeric alpha = 1.0;
        Numeric c1 = 1e-4;
        Numeric beta = 0.9;
        int max_line_search_iterations = 0;
        int ls_count = 0;

        while (ls_count++ < max_line_search_iterations)
        {
            Numeric u_trial = val_last - alpha * Hk * gk;
            u_trial = std::clamp(u_trial, Numeric(0), Numeric(1));

            if (Ri(u_trial).norm() <= current_residual - c1 * alpha * gk * Hk * gk)
            {
                return u_trial;
            }
            alpha *= beta;
        }
        Numeric u_final = val_last - alpha * Hk * gk;
        return std::clamp(u_final, Numeric(0), Numeric(1));
    };

    int count = 0;
    while (!((current_residual < epsilon) || (count++ >= max_iteration_count)))
    {
        Numeric gk = fi(val_last);
        Numeric u_new = line_search(gk);

        Numeric sk = u_new - val_last;
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

        val_last = u_new;
        current_residual = Ri(val_last).norm();
    }
    return direction == 0
               ? std::make_pair(constant, val_last)
               : std::make_pair(val_last, constant);
}

auto Surface::BinarySearchParameterOn(const Vec3& point, int direction, Numeric constant,
                                      Numeric epsilon, Numeric max_iteration_count) const
    -> std::pair<Numeric, Numeric>
{
    Numeric low = 0.0;
    Numeric high = 1.0;

    auto Ri = [&point, this, direction, constant](Numeric val) -> Numeric
    {
        Numeric u = direction == 0 ? constant : val;
        Numeric v = direction == 1 ? constant : val;
        return (Evaluate(u, v) - point).norm();
    };

    int count = 0;
    Numeric mid_residual;
    while (count++ < max_iteration_count)
    {
        Numeric mid = (low + high) / 2.0;

        mid_residual = Ri(mid);

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
            if (std::abs(high - low) < epsilon || mid_residual < epsilon)
            {
                return direction == 0
                           ? std::make_pair(constant, mid)
                           : std::make_pair(mid, constant);
            }
            low = left_mid;
            high = right_mid;
        }
    }

    Numeric final_val = (low + high) / 2.0;
    return direction == 0
               ? std::make_pair(constant, final_val)
               : std::make_pair(final_val, constant);
}


Surface Surface::InsertKnotU(Numeric knot_value) const
{
    Surface result{*this};
    int k = result.KnotsU.InsertKnot(knot_value);
    result.ControlPoints.InsertU(k, Vec4::Zero());
    const auto& knots = this->KnotsU.Values();
    vector<Numeric> alpha_list(result.DegreeU);
    for (int i = k - DegreeU + 1; i <= k; i++)
    {
        int idx = i - k + DegreeU - 1;
        alpha_list[idx] = (knot_value - knots[i]) / (knots[i + DegreeU] - knots[i]);
    }

    for (int j = 0; j < result.ControlPoints.VCount; j++)
    {
        for (int i = k - DegreeU + 1; i <= k; i++)
        {
            Vec4 point_i = ToHomo(this->ControlPoints.Get(i, j));
            Vec4 point_im1 = ToHomo(this->ControlPoints.Get(i - 1, j));
            int idx = i - k + DegreeU - 1;
            Numeric alpha = alpha_list[idx];
            Vec4& point_new = result.ControlPoints.Get(i, j);
            point_new = FromHomo(alpha * point_i + (1 - alpha) * point_im1);
        }
    }
    return result;
}

Surface Surface::InsertKnotU(Numeric knot_value, int times) const
{
    std::vector list(times, knot_value);
    return InsertKnotU(list);
}

Surface Surface::InsertKnotU(std::span<Numeric> knots_to_insert) const
{
    Surface result{*this};
    for (auto knot : knots_to_insert)
    {
        result = result.InsertKnotU(knot);
    }
    return result;
}

Surface Surface::InsertKnotV(Numeric knot_value) const
{
    Surface result{*this};
    int k = result.KnotsV.InsertKnot(knot_value);
    result.ControlPoints.InsertV(k, Vec4::Zero());
    const auto& knots = this->KnotsV.Values();
    vector<Numeric> alpha_list(result.DegreeV);
    for (int i = k - DegreeV + 1; i <= k; i++)
    {
        int idx = i - k + DegreeV - 1;
        alpha_list[idx] = (knot_value - knots[i]) / (knots[i + DegreeV] - knots[i]);
    }

    for (int j = 0; j < result.ControlPoints.UCount; j++)
    {
        for (int i = k - DegreeV + 1; i <= k; i++)
        {
            Vec4 point_i = ToHomo(this->ControlPoints.Get(j, i));
            Vec4 point_im1 = ToHomo(this->ControlPoints.Get(j, i - 1));

            int idx = i - k + DegreeV - 1;
            Numeric alpha = alpha_list[idx];
            Vec4& point_new = result.ControlPoints.Get(j, i);
            point_new = FromHomo(alpha * point_i + (1 - alpha) * point_im1);
        }
    }
    return result;
}

Surface Surface::InsertKnotV(Numeric knot_value, int times) const
{
    std::vector list(times, knot_value);
    return InsertKnotV(list);
}

Surface Surface::InsertKnotV(std::span<Numeric> knots_to_insert) const
{
    Surface result{*this};
    for (auto knot : knots_to_insert)
    {
        result = result.InsertKnotV(knot);
    }
    return result;
}

std::tuple<Surface, int> Surface::RemoveKnotU(Numeric knot_remove, int times, Numeric tolerance) const
{
    Surface result{*this};
    auto& points_ref = result.ControlPoints;
    points_ref = {points_ref.UCount - times, points_ref.VCount};
    int t = 0;
    for (int it_v = 0; it_v < ControlPoints.VCount; ++it_v)
    {
        auto knots = result.KnotsU;
        auto points = this->ControlPoints.GetV(it_v);
        int t_tmp = KnotRemoval(knots, points, result.DegreeU, knot_remove, times, tolerance);
        if (t_tmp == 0 || it_v != 0 && t != t_tmp) return {*this, 0};
        if (it_v == ControlPoints.VCount - 1) result.KnotsU = knots;
        points_ref.SetV(it_v, points);
        t = t_tmp;
    }
    return {result, t};
}

std::tuple<Surface, int> Surface::RemoveKnotV(Numeric knot_remove, int times, Numeric tolerance) const
{
    Surface result{*this};
    auto& points_ref = result.ControlPoints;
    points_ref = {points_ref.UCount, points_ref.VCount - times};
    int t = 0;
    for (int it_u = 0; it_u < ControlPoints.UCount; ++it_u)
    {
        auto knots = result.KnotsV;
        auto points = this->ControlPoints.GetU(it_u);
        int t_tmp = KnotRemoval(knots, points, result.DegreeV, knot_remove, times, tolerance);
        if (t_tmp == 0 || it_u != 0 && t != t_tmp) return {*this, 0};
        if (it_u == ControlPoints.UCount - 1) result.KnotsV = knots;
        points_ref.SetU(it_u, points);
        t = t_tmp;
    }
    return {result, t};
}

Surface Surface::ElevateDegreeU(int times) const
{
    Surface result{*this};
    auto& points_ref = result.ControlPoints;
    points_ref = {points_ref.UCount + times, points_ref.VCount};
    for (int it_v = 0; it_v < ControlPoints.VCount; ++it_v)
    {
        auto degree = result.DegreeU;
        auto knots = result.KnotsU;
        auto points = this->ControlPoints.GetV(it_v);
        NurbsDegreeElevation(knots, points, degree, times);
        if (it_v == ControlPoints.VCount - 1)
        {
            result.KnotsU = knots;
            result.DegreeU = degree;
        }
        points_ref.SetV(it_v, points);
    }
    return result;
}

Surface Surface::ElevateDegreeV(int times) const
{
    Surface result{*this};
    auto& points_ref = result.ControlPoints;
    points_ref = {points_ref.UCount, points_ref.VCount + times};
    for (int it_u = 0; it_u < ControlPoints.UCount; ++it_u)
    {
        auto degree = result.DegreeV;
        auto knots = result.KnotsV;
        auto points = this->ControlPoints.GetU(it_u);
        NurbsDegreeElevation(knots, points, degree, times);
        if (it_u == ControlPoints.UCount - 1)
        {
            result.KnotsV = knots;
            result.DegreeV = degree;
        }
        points_ref.SetU(it_u, points);
    }
    return result;
}

Surface Surface::Transform(const Mat3x3& R) const
{
    Surface transformed_surface = *this;
    for (auto& point : transformed_surface.ControlPoints.Values)
    {
        point.head<3>() = R * point.head<3>();
    }
    return transformed_surface;
}

Surface Surface::AlignParameterDomain(AlignAxis u_axis, AlignAxis v_axis)
{
    // 创建当前 Surface 的副本，以便进行操作
    Surface result = *this;

    // 检查控制点网格是否至少有 2 个控制点，以计算参数方向
    if (ControlPoints.UCount < 2 || ControlPoints.VCount < 2)
    {
        // 控制点太少，无法计算参数方向，直接返回原始 Surface
        return result;
    }

    // 1. 获取四个角点的控制点
    int u_max = ControlPoints.UCount - 1;
    int v_max = ControlPoints.VCount - 1;

    Vec3 S00 = ControlPoints.Get(0, 0).head<3>();
    Vec3 S10 = ControlPoints.Get(u_max, 0).head<3>();
    Vec3 S01 = ControlPoints.Get(0, v_max).head<3>();
    Vec3 S11 = ControlPoints.Get(u_max, v_max).head<3>();

    // 2. 计算包围盒
    Numeric min_x = std::min({S00.x(), S10.x(), S01.x(), S11.x()});
    Numeric max_x = std::max({S00.x(), S10.x(), S01.x(), S11.x()});
    Numeric min_y = std::min({S00.y(), S10.y(), S01.y(), S11.y()});
    Numeric max_y = std::max({S00.y(), S10.y(), S01.y(), S11.y()});
    Numeric min_z = std::min({S00.z(), S10.z(), S01.z(), S11.z()});
    Numeric max_z = std::max({S00.z(), S10.z(), S01.z(), S11.z()});

    Numeric dx = max_x - min_x;
    Numeric dy = max_y - min_y;
    Numeric dz = max_z - min_z;

    // 3. 确定主要的两个坐标轴和法向量
    std::vector<std::pair<Numeric, AlignAxis>> extents = {
        {dx, AlignAxis::XAxis},
        {dy, AlignAxis::YAxis},
        {dz, AlignAxis::ZAxis}
    };

    // 按照尺寸从大到小排序
    std::sort(extents.begin(),
              extents.end(),
              [](const auto& a, const auto& b)
              {
                  return a.first > b.first;
              });

    // 主要的两个方向（参数方向）
    AlignAxis major_axis1 = extents[0].second;
    AlignAxis major_axis2 = extents[1].second;
    // 法向量方向
    AlignAxis normal_axis = extents[2].second;

    // 4. 计算控制点之间的矢量
    struct VectorInfo
    {
        Vec3 vector;
        // 角度
        Numeric angle1; // 与 major_axis1 的夹角
        Numeric angle2; // 与 major_axis2 的夹角
    };
    std::vector<VectorInfo> vectors;

    auto compute_angle = [](const Vec3& vec, AlignAxis axis) -> Numeric
    {
        Vec3 axis_vector;
        switch (axis)
        {
        case AlignAxis::XAxis:
            axis_vector = Vec3::UnitX();
            break;
        case AlignAxis::YAxis:
            axis_vector = Vec3::UnitY();
            break;
        case AlignAxis::ZAxis:
            axis_vector = Vec3::UnitZ();
            break;
        }
        // 计算夹角的余弦值
        Numeric cos_angle = vec.normalized().dot(axis_vector);
        // 为了比较方便，使用夹角的绝对值
        return std::acos(std::abs(cos_angle));
    };

    // 计算非零矢量并计算其与主要方向的夹角
    auto add_vector_info = [&](const Vec3& vec)
    {
        if (vec.norm() > 1e-6)
        {
            VectorInfo info;
            info.vector = vec;
            info.angle1 = compute_angle(vec, major_axis1);
            info.angle2 = compute_angle(vec, major_axis2);
            vectors.push_back(info);
        }
    };

    add_vector_info(S10 - S00);
    add_vector_info(S01 - S00);
    add_vector_info(S11 - S10);
    add_vector_info(S11 - S01);

    if (vectors.empty())
    {
        // 无法确定参数方向，返回原始 Surface
        return result;
    }

    // 5. 为每个矢量分配最接近的主要方向
    AlignAxis current_u_axis, current_v_axis;
    int u_sign = 0, v_sign = 0;

    if (vectors.size() < 4)
    {
        int count_axis1 = 0;
        int count_axis2 = 0;
        Numeric du = 0;
        Numeric dv = 0;

        auto GetAxisValue = [](const Vec3& point, AlignAxis axis) -> Numeric
        {
            switch (axis)
            {
            case AlignAxis::XAxis:
                return point.x();
            case AlignAxis::YAxis:
                return point.y();
            case AlignAxis::ZAxis:
                return point.z();
            }
            return 0.0;
        };

        for (const auto& info : vectors)
        {
            if (info.angle1 < info.angle2)
            {
                // 分配给 major_axis1
                ++count_axis1;
                // 累加在该轴上的增量
                du += GetAxisValue(info.vector, major_axis1);
            }
            else
            {
                // 分配给 major_axis2
                ++count_axis2;
                dv += GetAxisValue(info.vector, major_axis2);
            }
        }

        if (count_axis1 >= count_axis2)
        {
            current_u_axis = major_axis1;
            current_v_axis = major_axis2;
            u_sign = (du >= 0) ? 1 : -1;
            v_sign = (dv >= 0) ? 1 : -1;
        }
        else
        {
            current_u_axis = major_axis2;
            current_v_axis = major_axis1;
            u_sign = (dv >= 0) ? 1 : -1;
            v_sign = (du >= 0) ? 1 : -1;
        }
    }
    else
    {
        auto GetDominantAxis = [](const Vec3& vec) -> std::pair<AlignAxis, int>
        {
            Numeric abs_x = std::abs(vec.x());
            Numeric abs_y = std::abs(vec.y());
            Numeric abs_z = std::abs(vec.z());
            Numeric max_abs = std::max({abs_x, abs_y, abs_z});

            AlignAxis axis;
            if (max_abs == abs_x)
                axis = AlignAxis::XAxis;
            else if (max_abs == abs_y)
                axis = AlignAxis::YAxis;
            else
                axis = AlignAxis::ZAxis;

            int sign = 0;
            switch (axis)
            {
            case AlignAxis::XAxis:
                sign = (vec.x() >= 0) ? 1 : -1;
                break;
            case AlignAxis::YAxis:
                sign = (vec.y() >= 0) ? 1 : -1;
                break;
            case AlignAxis::ZAxis:
                sign = (vec.z() >= 0) ? 1 : -1;
                break;
            }

            return {axis, sign};
        };

        std::tie(current_u_axis, u_sign) = GetDominantAxis(vectors[0].vector);
        std::tie(current_v_axis, v_sign) = GetDominantAxis(vectors[1].vector);
    }


    // 7. 判断是否需要交换 u 和 v
    bool need_swap = false;
    if (current_u_axis != u_axis)
    {
        if (current_v_axis == u_axis && current_u_axis == v_axis)
        {
            // 需要交换 u 和 v
            need_swap = true;
            std::swap(current_u_axis, current_v_axis);
            std::swap(u_sign, v_sign);
        }
        else
        {
            // 无法通过交换实现对齐，返回原始 Surface
            return result;
        }
    }

    // 8. 判断是否需要反转 u 和 v 的方向
    bool reverse_u = false;
    if (current_u_axis == u_axis)
    {
        if (u_sign < 0)
            reverse_u = true;
    }

    bool reverse_v = false;
    if (current_v_axis == v_axis)
    {
        if (v_sign < 0)
            reverse_v = true;
    }

    // 9. 执行必要的操作：交换和反转
    if (need_swap)
    {
        // 交换次数
        std::swap(result.DegreeU, result.DegreeV);
        // 交换节点矢量
        std::swap(result.KnotsU, result.KnotsV);

        // 交换控制点网格
        auto& grid = result.ControlPoints;
        Grid<Vec4> transposed_grid(grid.VCount, grid.UCount);
        for (int u = 0; u < grid.UCount; ++u)
        {
            for (int v = 0; v < grid.VCount; ++v)
            {
                transposed_grid.Get(v, u) = grid.Get(u, v);
            }
        }
        result.ControlPoints = std::move(transposed_grid);
    }

    if (reverse_u)
    {
        // 反转 u 方向的节点矢量
        result.KnotsU.Reverse();

        // 反转 u 方向的控制点顺序
        auto& grid = result.ControlPoints;
        for (int v = 0; v < grid.VCount; ++v)
        {
            for (int u = 0; u < grid.UCount / 2; ++u)
            {
                std::swap(grid.Get(u, v), grid.Get(grid.UCount - 1 - u, v));
            }
        }
    }

    if (reverse_v)
    {
        // 反转 v 方向的节点矢量
        result.KnotsV.Reverse();

        // 反转 v 方向的控制点顺序
        auto& grid = result.ControlPoints;
        for (int u = 0; u < grid.UCount; ++u)
        {
            for (int v = 0; v < grid.VCount / 2; ++v)
            {
                std::swap(grid.Get(u, v), grid.Get(u, grid.VCount - 1 - v));
            }
        }
    }

    // 10. 返回调整后的 Surface
    return result;
}


Grid<Vec4> Surface::HomogeneousDerivative(Numeric u, Numeric v, int order_u, int order_v) const
{
    assert(u >= 0 && u <= 1);
    assert(v >= 0 && v <= 1);
    int index_span_u = KnotsU.FindSpanIndex(DegreeU, u);
    int index_span_v = KnotsV.FindSpanIndex(DegreeV, v);
    MatX basis_u = BSplineBasis::EvaluateAll(DegreeU, KnotsU.Values(), index_span_u, u, order_u);
    MatX basis_v = BSplineBasis::EvaluateAll(DegreeV, KnotsV.Values(), index_span_v, v, order_v);
    assert(basis_u.cols() == DegreeU + 1);
    assert(basis_v.cols() == DegreeV + 1);

    int index_pre_u = index_span_u - DegreeU;
    int index_pre_v = index_span_v - DegreeV;

    Grid<Vec4> result(order_u + 1, order_v + 1, Vec4::Zero());

    for (int l = 0; l <= order_v; ++l)
    {
        for (int k = 0; k <= order_u; ++k)
        {
            for (int j = 0; j <= DegreeV; j++)
            {
                int index_v = index_pre_v + j;
                Vec4 tmp = Vec4::Zero();
                for (int i = 0; i <= DegreeU; i++)
                {
                    int index_u = index_pre_u + i;
                    auto point = ToHomo(ControlPoints.Get(index_u, index_v));
                    tmp.noalias() += basis_u(k, i) * point;
                }
                result.Get(k, l).noalias() += basis_v(l, j) * tmp;
            }
        }
    }
    return result;
}

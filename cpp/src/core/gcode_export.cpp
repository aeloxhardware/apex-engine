// Voir apex/gcode_export.hpp : transcription terme a terme de
// file_io/gcode_export.py::export_gcode (corps), texte identique.

#include "apex/gcode_export.hpp"

#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <limits>

namespace apex {

namespace {

constexpr double kDegToRad = 3.14159265358979323846 / 180.0;   // np.radians
constexpr double kRadToDeg = 180.0 / 3.14159265358979323846;   // np.degrees

// Modulo flottant de Python (signe du diviseur).
inline double py_mod(double a, double b) {
    double m = std::fmod(a, b);
    if (m != 0.0) {
        if ((b < 0.0) != (m < 0.0)) m += b;
    } else {
        m = std::copysign(0.0, b);
    }
    return m;
}

// f"{v:.{dec}f}".rstrip("0").rstrip(".") ; "-0" et "" -> "0"
void fmt(double v, int dec, std::string& out) {
    char buf[64];
    auto r = std::to_chars(buf, buf + sizeof(buf), v, std::chars_format::fixed, dec);
    std::string s(buf, r.ptr);
    while (!s.empty() && s.back() == '0') s.pop_back();
    if (!s.empty() && s.back() == '.') s.pop_back();
    if (s.empty() || s == "-0") s = "0";
    out += s;
}

struct Fw {
    double x, y, z, b;
};

inline double norm3(double x, double y, double z) { return std::sqrt(x * x + y * y + z * z); }
inline double norm4(double a, double b, double c, double d) { return std::sqrt(a * a + b * b + c * c + d * d); }

}  // namespace

std::string export_gcode_body(const double* pos, const double* ang, const double* ext, const std::uint8_t* skin,
                              const double* slow, std::int64_t n, const std::vector<double>& zgrid,
                              const std::vector<double>& zvals, const GcodeExportParams& P,
                              GcodeExportStats& st) {
    std::string out;
    out.reserve(static_cast<std::size_t>(n) * 60);
    const double inf = std::numeric_limits<double>::infinity();
    st.tip_z_min = inf;
    st.head_z_min = inf;
    const bool has_table = !zgrid.empty();

    auto head_zlow = [&](double b) {
        std::size_t best = 0;
        double bd = std::fabs(zgrid[0] - b);
        for (std::size_t k = 1; k < zgrid.size(); ++k) {
            const double d = std::fabs(zgrid[k] - b);
            if (d < bd) { bd = d; best = k; }
        }
        return zvals[best];
    };

    auto to_fw = [&](const double* p, double b, const double* prev_y) {
        const double x = p[0], y = p[1], z = p[2];
        const double r = std::hypot(x, y);
        double fy;
        if (r < P.singularity_radius && prev_y != nullptr) {
            fy = *prev_y;
        } else {
            fy = std::atan2(y, x) * kRadToDeg;
            if (prev_y != nullptr) fy = *prev_y + (py_mod(fy - *prev_y + 180.0, 360.0) - 180.0);
        }
        const double br = b * kDegToRad;
        const double sr = std::sin(br), cz = std::cos(br);
        const double zoff = P.z_reference_tip ? (cz - 1.0) : cz;
        return Fw{r + P.tilt_sign * P.pivot_to_tip * sr, fy, z + P.pivot_to_tip * zoff, b};
    };

    bool have_prev = false;
    Fw prev{};
    double tip[3] = {0, 0, 0};
    double e = 0.0;

    auto move = [&](const double* point, double b, double tip_speed, double e_add, const std::string& comment) {
        const Fw fw = to_fw(point, b, have_prev ? &prev.y : nullptr);
        const double tz = point[2];
        st.tip_z_min = std::min(st.tip_z_min, tz);
        st.tip_below_plate += tz < 0.0 ? 1 : 0;
        if (has_table) {
            const double hz = tz + head_zlow(b);
            st.head_z_min = std::min(st.head_z_min, hz);
            st.head_below_margin += hz < P.platform_margin ? 1 : 0;
        }
        if (!st.any) {
            st.x_min = st.x_max = fw.x; st.y_min = st.y_max = fw.y; st.z_min = st.z_max = fw.z; st.any = true;
        } else {
            st.x_min = std::min(st.x_min, fw.x); st.x_max = std::max(st.x_max, fw.x);
            st.y_min = std::min(st.y_min, fw.y); st.y_max = std::max(st.y_max, fw.y);
            st.z_min = std::min(st.z_min, fw.z); st.z_max = std::max(st.z_max, fw.z);
        }
        const double tip_d = have_prev ? norm3(point[0] - tip[0], point[1] - tip[1], point[2] - tip[2]) : 0.0;
        const double fw_d = have_prev ? norm4(fw.x - prev.x, fw.y - prev.y, fw.z - prev.z, fw.b - prev.b) : 0.0;
        double feed = tip_d > 1e-9 ? tip_speed * 60.0 * fw_d / tip_d : P.rotation_feed;
        feed = std::min(std::max(feed, 1.0), P.max_feed);
        e += e_add;
        out += "G00 X"; fmt(fw.x, P.decimals, out);
        out += " Y"; fmt(fw.y, P.decimals, out);
        out += " Z"; fmt(fw.z, P.decimals, out);
        out += " B"; fmt(fw.b, P.decimals, out);
        if (e_add != 0.0) { out += " E"; fmt(e, P.decimals, out); }
        out += " F"; fmt(feed, P.decimals, out);
        if (!comment.empty()) { out += " ; "; out += comment; }
        out += '\n';
        prev = fw;
        have_prev = true;
        tip[0] = point[0]; tip[1] = point[1]; tip[2] = point[2];
        st.moves += 1;
    };

    auto retract = [&](double sign) {
        e += sign * P.retract_length;
        out += "G01 E"; fmt(e, P.decimals, out);
        out += " F"; fmt(P.retract_speed * 60.0, P.decimals, out);
        out += '\n';
    };

    // clearance_move : rend l'axe ('Z', 'R', 'n' pour none) et le point degage
    auto clearance = [&](const double* p, double b, std::array<double, 3>& cand) -> char {
        const double br = b * kDegToRad;
        const double sr = std::sin(br), cz = std::cos(br);
        const double x = p[0], y = p[1];
        const double r = std::hypot(x, y);
        double rad[3] = {1.0, 0.0, 0.0};
        if (r > 1e-9) { rad[0] = x / r; rad[1] = y / r; rad[2] = 0.0 / r; }
        const char order[2] = {std::fabs(cz) >= std::fabs(sr) ? 'Z' : 'R', std::fabs(cz) >= std::fabs(sr) ? 'R' : 'Z'};
        for (char axis : order) {
            double step[3];
            if (axis == 'Z') {
                step[0] = 0.0; step[1] = 0.0; step[2] = P.travel_clearance * (cz >= 0 ? 1.0 : -1.0);
            } else {
                if (r < 1e-9 && sr < 0) continue;
                const double s = sr >= 0 ? 1.0 : -1.0;
                for (int k = 0; k < 3; ++k) step[k] = rad[k] * P.travel_clearance * s;
            }
            for (int k = 0; k < 3; ++k) cand[k] = p[k] + step[k];
            const double lowest = cand[2] + (has_table ? head_zlow(b) : 0.0);
            if (lowest >= P.platform_margin || cand[2] >= p[2]) return axis;
        }
        cand = {p[0], p[1], p[2]};
        return 'n';
    };
    auto count_axis = [&](char a) {
        if (a == 'Z') st.clearance_z++; else if (a == 'R') st.clearance_r++; else st.clearance_none++;
    };
    auto axis_name = [](char a) { return a == 'Z' ? std::string("Z") : a == 'R' ? std::string("R") : std::string("none"); };

    int last_kind = -1;
    std::array<double, 3> up{}, up0{}, up1{};
    for (std::int64_t i = 0; i < n; ++i) {
        const double* p = pos + i * 3;
        const int kind = skin[i] ? 1 : 0;
        const double b = ang[i * 2 + 1];
        if (i > 0) {
            const Fw f = to_fw(p, b, nullptr);
            if (std::hypot(p[0], p[1]) > 1e-3) {
                st.c_mismatch_deg = std::max(
                    st.c_mismatch_deg, std::fabs(py_mod(f.y - (ang[i * 2] - 90.0) + 180.0, 360.0) - 180.0));
            }
        }
        if (kind != last_kind) {
            out += kind ? "; --- lissage ---\n" : "; --- corps ---\n";
            last_kind = kind;
        }
        if (i == 0 || ext[i] <= 0) {
            st.travels++;
            if (i == 0) {
                const char a = clearance(p, b, up);
                move(up.data(), b, P.travel_speed, 0.0, "approche");
                count_axis(a);
                move(p, b, P.travel_speed, 0.0, "");
                continue;
            }
            const double* q = pos + (i - 1) * 3;
            const double length = norm3(p[0] - q[0], p[1] - q[1], p[2] - q[2]);
            if (length >= P.retract_min_travel) {
                retract(-1.0);
                st.retractions++;
                const double b0 = ang[(i - 1) * 2 + 1];
                const char a0 = clearance(q, b0, up0);
                clearance(p, b, up1);
                count_axis(a0);
                move(up0.data(), b0, P.travel_speed, 0.0, "degagement " + axis_name(a0));
                move(up1.data(), b, P.travel_speed, 0.0, "");
                move(p, b, P.travel_speed, 0.0, "");
                retract(+1.0);
            } else {
                move(p, b, P.travel_speed, 0.0, "");
            }
            continue;
        }
        const double* q = pos + (i - 1) * 3;
        const double seg = norm3(p[0] - q[0], p[1] - q[1], p[2] - q[2]);
        const double h = skin[i] ? P.skin_thickness : P.layer_height;
        const double e_add = seg * h * P.line_width * P.flow / P.fil_area;
        st.filament_mm += e_add;
        const double speed = P.print_speed * slow[i - 1];
        move(p, b, speed, e_add, "");
    }
    if (n) {
        retract(-1.0);
        const double* p = pos + (n - 1) * 3;
        const double b = ang[(n - 1) * 2 + 1];
        const char a = clearance(p, b, up);
        move(up.data(), b, P.travel_speed, 0.0, "degagement final " + axis_name(a));
    }
    if (!out.empty() && out.back() == '\n') out.pop_back();
    return out;
}

namespace {
inline int np_sign(double v) { return v > 0.0 ? 1 : (v < 0.0 ? -1 : 0); }   // np.sign (NaN exclu)
}

void feed_plan_caps(const double* d, std::int64_t count, double vmax, double amax, double* cap) {
    for (std::int64_t i = 0; i < count; ++i) cap[i] = vmax;
    double prev = 0.0;
    for (std::int64_t i = 0; i < count; ++i) {
        if (i && np_sign(d[i]) != np_sign(d[i - 1])) prev = 0.0;
        cap[i] = std::min(cap[i], std::sqrt(std::pow(prev, 2.0) + 2 * amax * std::fabs(d[i])));
        prev = cap[i];
    }
    double nxt = 0.0;
    for (std::int64_t i = count - 1; i >= 0; --i) {
        if (i < count - 1 && np_sign(d[i]) != np_sign(d[i + 1])) nxt = 0.0;
        cap[i] = std::min(cap[i], std::sqrt(std::pow(nxt, 2.0) + 2 * amax * std::fabs(d[i])));
        nxt = cap[i];
    }
}

}  // namespace apex

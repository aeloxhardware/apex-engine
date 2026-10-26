#include "apex/step8_assembly.hpp"

#include <cmath>

namespace apex {

RawPointsOutput build_raw_points(std::int64_t layer_count, const double* frames,
                                 const std::uint8_t* reversed, const std::int64_t* layer_ids,
                                 const std::int64_t* poly_begin, const std::int64_t* point_begin,
                                 const double* points, std::int64_t columns, std::int64_t segment_count,
                                 const std::int64_t* seg_layer, const std::int64_t* seg_index,
                                 const std::int64_t* seg_pass) {
    RawPointsOutput out;
    // Normales de couche, une fois par couche (comme le cache du Python).
    std::vector<double> normals(static_cast<std::size_t>(layer_count) * 3);
    std::vector<double> frame_normals(static_cast<std::size_t>(layer_count) * 3);
    for (std::int64_t l = 0; l < layer_count; ++l) {
        const double* u = frames + 9 * l + 3;
        const double* v = frames + 9 * l + 6;
        const double c[3] = {u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]};
        for (int a = 0; a < 3; ++a) frame_normals[3 * static_cast<std::size_t>(l) + a] = c[a];
        const double n = std::sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2]);
        double d[3] = {0.0, 0.0, 1.0};
        if (!(n < 1e-12)) {
            d[0] = c[0] / n; d[1] = c[1] / n; d[2] = c[2] / n;
        }
        if (reversed[l]) {
            d[0] = -d[0]; d[1] = -d[1]; d[2] = -d[2];
        }
        for (int a = 0; a < 3; ++a) normals[3 * static_cast<std::size_t>(l) + a] = d[a];
    }

    bool has_previous = false;
    std::int64_t prev_layer = 0, prev_seg = 0;
    auto push = [&](const double p[3], std::int64_t table, std::int64_t seg, std::int64_t pass, bool travel) {
        out.position.insert(out.position.end(), p, p + 3);
        out.desired_direction.insert(out.desired_direction.end(), normals.begin() + 3 * table,
                                     normals.begin() + 3 * table + 3);
        out.layer_id.push_back(layer_ids[table]);
        out.segment_index.push_back(seg);
        out.pass_index.push_back(pass);
        out.is_travel.push_back(travel ? 1 : 0);
    };

    for (std::int64_t s = 0; s < segment_count; ++s) {
        const std::int64_t table = seg_layer[s];
        const std::int64_t seg = seg_index[s];
        const std::int64_t first = point_begin[poly_begin[table] + seg];
        const std::int64_t last = point_begin[poly_begin[table] + seg + 1];
        if (last == first) continue;
        const double* origin = frames + 9 * table;
        const double* u = frames + 9 * table + 3;
        const double* v = frames + 9 * table + 6;
        const double* nrm = &frame_normals[3 * static_cast<std::size_t>(table)];
        const bool new_key = has_previous && !(layer_ids[table] == prev_layer && seg == prev_seg);
        for (std::int64_t k = first; k < last; ++k) {
            const double* q = points + k * columns;
            double p[3];
            for (int a = 0; a < 3; ++a) {
                p[a] = origin[a] + q[0] * u[a];
                p[a] = p[a] + q[1] * v[a];
                if (columns >= 3) p[a] = p[a] + q[2] * nrm[a];
            }
            if (k == first && new_key) push(p, table, seg, seg_pass[s], true);
            push(p, table, seg, seg_pass[s], false);
        }
        has_previous = true;
        prev_layer = layer_ids[table];
        prev_seg = seg;
    }
    return out;
}

}  // namespace apex

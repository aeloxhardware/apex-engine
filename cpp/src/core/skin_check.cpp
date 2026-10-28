// Voir apex/skin_check.hpp.
#include "apex/skin_check.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <thread>

namespace apex {

namespace {

// Segment (a, b) de dense() sans son extremite b : a + (k/m)(b - a), k = 0..m-1.
template <class F>
void segment_points(const double* a, const double* b, double step, F&& emit) {
    const double dx = b[0] - a[0], dy = b[1] - a[1], dz = b[2] - a[2];
    const double seg = std::sqrt(dx * dx + dy * dy + dz * dz);
    const std::int64_t m = std::max<std::int64_t>(1, static_cast<std::int64_t>(std::ceil(seg / step)));
    for (std::int64_t k = 0; k < m; ++k) {
        const double t = static_cast<double>(k) / static_cast<double>(m);
        emit(a[0] + t * dx, a[1] + t * dy, a[2] + t * dz);
    }
}

}  // namespace

SkinCheckIndex::SkinCheckIndex(double cell) : cell_(cell) {}

std::int64_t SkinCheckIndex::cell_key(std::int64_t i, std::int64_t j, std::int64_t k) const {
    // Cles exactes (pas de collision de hachage) : indices de cellule bornes a +-2^20.
    const std::int64_t off = std::int64_t(1) << 20;
    return ((i + off) << 42) | ((j + off) << 21) | (k + off);
}

void SkinCheckIndex::add_point(Layer& layer, double x, double y, double z) const {
    const std::int64_t id = static_cast<std::int64_t>(layer.xyz.size() / 3);
    layer.xyz.push_back(x);
    layer.xyz.push_back(y);
    layer.xyz.push_back(z);
    const auto ci = static_cast<std::int64_t>(std::floor(x / cell_));
    const auto cj = static_cast<std::int64_t>(std::floor(y / cell_));
    const auto ck = static_cast<std::int64_t>(std::floor(z / cell_));
    layer.grid[cell_key(ci, cj, ck)].push_back(id);
}

void SkinCheckIndex::add_polylines(const double* pts, const std::int64_t* starts, std::int64_t n_lines, double step) {
    for (std::int64_t l = 0; l < n_lines; ++l) {
        const std::int64_t s = starts[l], e = starts[l + 1];
        if (e - s < 2) {
            for (std::int64_t i = s; i < e; ++i) add_point(perm_, pts[i * 3], pts[i * 3 + 1], pts[i * 3 + 2]);
            continue;
        }
        for (std::int64_t j = s; j + 1 < e; ++j)
            segment_points(pts + j * 3, pts + (j + 1) * 3, step,
                           [&](double x, double y, double z) { add_point(perm_, x, y, z); });
        const double* q = pts + (e - 1) * 3;
        add_point(perm_, q[0], q[1], q[2]);
    }
}

double SkinCheckIndex::nearest(const Layer& layer, double x, double y, double z, double bound) const {
    double best2 = std::numeric_limits<double>::infinity();
    const auto ci = static_cast<std::int64_t>(std::floor(x / cell_));
    const auto cj = static_cast<std::int64_t>(std::floor(y / cell_));
    const auto ck = static_cast<std::int64_t>(std::floor(z / cell_));
    for (std::int64_t di = -1; di <= 1; ++di)
        for (std::int64_t dj = -1; dj <= 1; ++dj)
            for (std::int64_t dk = -1; dk <= 1; ++dk) {
                const auto it = layer.grid.find(cell_key(ci + di, cj + dj, ck + dk));
                if (it == layer.grid.end()) continue;
                for (const std::int64_t id : it->second) {
                    const double dx = x - layer.xyz[id * 3];
                    const double dy = y - layer.xyz[id * 3 + 1];
                    const double dz = z - layer.xyz[id * 3 + 2];
                    const double d2 = dx * dx + dy * dy + dz * dz;
                    if (d2 < best2) best2 = d2;
                }
            }
    const double d = std::sqrt(best2);
    return d < bound ? d : std::numeric_limits<double>::infinity();
}

void SkinCheckIndex::check_zone(const double* head, std::int64_t n_head, const double* axes, std::int64_t n_axes,
                                const double* pts, const std::int64_t* starts, std::int64_t n_lines,
                                const double* cosv, const double* sinv, const std::int64_t* tested,
                                std::int64_t n_tested, double plate_radius, double margin, double reach, double step,
                                std::uint8_t* plate_hit, double* distance) const {
    Layer tmp;
    const double bound = reach + 1.0;
    std::vector<double> verts(static_cast<std::size_t>(n_head * 3));
    std::int64_t t = 0;
    for (std::int64_t l = 0; l < n_lines; ++l) {
        const std::int64_t s = starts[l], e = starts[l + 1];
        std::int64_t seg_done = s;   // segments (j, j+1) deja ajoutes pour j < seg_done
        auto add_segments_upto = [&](std::int64_t j_end) {
            for (; seg_done < j_end; ++seg_done)
                segment_points(pts + seg_done * 3, pts + (seg_done + 1) * 3, step,
                               [&](double x, double y, double z) { add_point(tmp, x, y, z); });
        };
        while (t < n_tested && tested[t] < e) {
            const std::int64_t i = tested[t];
            if (i > s) {
                // dense(p[:i]) = segments s..i-2 entiers + le point p[i-1]
                add_segments_upto(i - 1);
                const double* q = pts + (i - 1) * 3;
                add_point(tmp, q[0], q[1], q[2]);
                // Particularite de la recherche (blocs de 256 points) : quand le point teste
                // est EXACTEMENT un debut de bloc (>= 256), l'arbre du bloc est dense(p[:i+1]) :
                // il contient deja le segment (i-1, i) et p[i]. Reproduit tel quel.
                const std::int64_t local = i - s;
                if (local >= 256 && local % 256 == 0) {
                    add_segments_upto(i);
                    add_point(tmp, pts[i * 3], pts[i * 3 + 1], pts[i * 3 + 2]);
                }
            }
            for (std::int64_t h = 0; h < n_head * 3; ++h) verts[h] = head[h];
            for (std::int64_t ax = n_axes - 1; ax >= 0; --ax) {
                const double kx = axes[ax * 3], ky = axes[ax * 3 + 1], kz = axes[ax * 3 + 2];
                const double c = cosv[i * n_axes + ax], sn = sinv[i * n_axes + ax];
                for (std::int64_t h = 0; h < n_head; ++h) {
                    double* v = &verts[h * 3];
                    const double dot = v[0] * kx + v[1] * ky + v[2] * kz;
                    const double cx = ky * v[2] - kz * v[1];
                    const double cy = kz * v[0] - kx * v[2];
                    const double cz = kx * v[1] - ky * v[0];
                    const double nx = v[0] * c + cx * sn + kx * dot * (1.0 - c);
                    const double ny = v[1] * c + cy * sn + ky * dot * (1.0 - c);
                    const double nz = v[2] * c + cz * sn + kz * dot * (1.0 - c);
                    v[0] = nx;
                    v[1] = ny;
                    v[2] = nz;
                }
            }
            const double* p = pts + i * 3;
            bool plate = false;
            double d = std::numeric_limits<double>::infinity();
            for (std::int64_t h = 0; h < n_head; ++h) {
                const double x = verts[h * 3] + p[0], y = verts[h * 3 + 1] + p[1], z = verts[h * 3 + 2] + p[2];
                if (z < margin && std::sqrt(x * x + y * y) < plate_radius + margin) plate = true;
                d = std::min(d, nearest(perm_, x, y, z, bound));
                if (!tmp.xyz.empty()) d = std::min(d, nearest(tmp, x, y, z, bound));
            }
            plate_hit[t] = plate ? 1 : 0;
            distance[t] = d;
            ++t;
        }
        // ligne finie : toute la ligne densifiee rejoint les obstacles de la zone
        if (e - s >= 2) {
            add_segments_upto(e - 1);
            const double* q = pts + (e - 1) * 3;
            add_point(tmp, q[0], q[1], q[2]);
        } else if (e - s == 1) {
            add_point(tmp, pts[s * 3], pts[s * 3 + 1], pts[s * 3 + 2]);
        }
    }
}

void SkinCheckIndex::check_zone_parallel(const double* head, std::int64_t n_head, const double* axes,
                                         std::int64_t n_axes, const double* pts, const std::int64_t* starts,
                                         std::int64_t n_lines, const double* cosv, const double* sinv,
                                         const std::int64_t* tested, std::int64_t n_tested, double plate_radius,
                                         double margin, double reach, double step, std::uint8_t* plate_hit,
                                         double* distance, int threads) const {
    // 1. Couche de la zone dans l'ordre d'insertion final (lignes entieres, segment
    //    par segment, puis le dernier point), et limite de prefixe par point teste.
    Layer zone;
    std::vector<std::int64_t> limit(static_cast<std::size_t>(n_tested), 0);
    std::vector<std::int64_t> seg_first;
    std::int64_t t = 0;
    for (std::int64_t l = 0; l < n_lines; ++l) {
        const std::int64_t s = starts[l], e = starts[l + 1];
        const std::int64_t line_start = static_cast<std::int64_t>(zone.xyz.size() / 3);
        seg_first.assign(static_cast<std::size_t>(std::max<std::int64_t>(e - s, 0)), 0);
        if (e - s >= 2) {
            for (std::int64_t j = s; j + 1 < e; ++j) {
                seg_first[j - s] = static_cast<std::int64_t>(zone.xyz.size() / 3);
                segment_points(pts + j * 3, pts + (j + 1) * 3, step,
                               [&](double x, double y, double z) { add_point(zone, x, y, z); });
            }
            seg_first[e - 1 - s] = static_cast<std::int64_t>(zone.xyz.size() / 3);   // dernier point
            const double* q = pts + (e - 1) * 3;
            add_point(zone, q[0], q[1], q[2]);
        } else if (e - s == 1) {
            seg_first[0] = static_cast<std::int64_t>(zone.xyz.size() / 3);
            add_point(zone, pts[s * 3], pts[s * 3 + 1], pts[s * 3 + 2]);
        }
        while (t < n_tested && tested[t] < e) {
            const std::int64_t i = tested[t];
            if (i <= s) {
                limit[t] = line_start;                          // lignes precedentes seulement
            } else {
                const std::int64_t local = i - s;
                // dense(p[:i]) : segments s..i-2 + p[i-1] (premier point du segment i-1) ;
                // blocs de 256 : + segment i-1 + p[i] (premier point du segment i, ou le dernier point).
                const std::int64_t last = (local >= 256 && local % 256 == 0) ? i : i - 1;
                limit[t] = seg_first[last - s] + 1;
            }
            ++t;
        }
    }
    // 2. Points testes, independants : repartis sur les fils.
    const double bound = reach + 1.0;
    auto nearest_prefix = [&](double x, double y, double z, std::int64_t lim) {
        double best2 = std::numeric_limits<double>::infinity();
        if (lim <= 0) return std::numeric_limits<double>::infinity();
        const auto ci = static_cast<std::int64_t>(std::floor(x / cell_));
        const auto cj = static_cast<std::int64_t>(std::floor(y / cell_));
        const auto ck = static_cast<std::int64_t>(std::floor(z / cell_));
        for (std::int64_t di = -1; di <= 1; ++di)
            for (std::int64_t dj = -1; dj <= 1; ++dj)
                for (std::int64_t dk = -1; dk <= 1; ++dk) {
                    const auto it = zone.grid.find(cell_key(ci + di, cj + dj, ck + dk));
                    if (it == zone.grid.end()) continue;
                    for (const std::int64_t id : it->second) {
                        if (id >= lim) break;                   // rangs croissants dans chaque cellule
                        const double dx = x - zone.xyz[id * 3];
                        const double dy = y - zone.xyz[id * 3 + 1];
                        const double dz = z - zone.xyz[id * 3 + 2];
                        const double d2 = dx * dx + dy * dy + dz * dz;
                        if (d2 < best2) best2 = d2;
                    }
                }
        const double d = std::sqrt(best2);
        return d < bound ? d : std::numeric_limits<double>::infinity();
    };
    auto work = [&](std::int64_t begin, std::int64_t end) {
        std::vector<double> verts(static_cast<std::size_t>(n_head * 3));
        for (std::int64_t tt = begin; tt < end; ++tt) {
            const std::int64_t i = tested[tt];
            for (std::int64_t h = 0; h < n_head * 3; ++h) verts[h] = head[h];
            for (std::int64_t ax = n_axes - 1; ax >= 0; --ax) {
                const double kx = axes[ax * 3], ky = axes[ax * 3 + 1], kz = axes[ax * 3 + 2];
                const double c = cosv[i * n_axes + ax], sn = sinv[i * n_axes + ax];
                for (std::int64_t h = 0; h < n_head; ++h) {
                    double* v = &verts[h * 3];
                    const double dot = v[0] * kx + v[1] * ky + v[2] * kz;
                    const double cx = ky * v[2] - kz * v[1];
                    const double cy = kz * v[0] - kx * v[2];
                    const double cz = kx * v[1] - ky * v[0];
                    const double nx = v[0] * c + cx * sn + kx * dot * (1.0 - c);
                    const double ny = v[1] * c + cy * sn + ky * dot * (1.0 - c);
                    const double nz = v[2] * c + cz * sn + kz * dot * (1.0 - c);
                    v[0] = nx;
                    v[1] = ny;
                    v[2] = nz;
                }
            }
            const double* p = pts + i * 3;
            bool plate = false;
            double d = std::numeric_limits<double>::infinity();
            for (std::int64_t h = 0; h < n_head; ++h) {
                const double x = verts[h * 3] + p[0], y = verts[h * 3 + 1] + p[1], z = verts[h * 3 + 2] + p[2];
                if (z < margin && std::sqrt(x * x + y * y) < plate_radius + margin) plate = true;
                d = std::min(d, nearest(perm_, x, y, z, bound));
                d = std::min(d, nearest_prefix(x, y, z, limit[tt]));
            }
            plate_hit[tt] = plate ? 1 : 0;
            distance[tt] = d;
        }
    };
    unsigned count = threads > 0 ? static_cast<unsigned>(threads) : std::max(1u, std::thread::hardware_concurrency());
    count = static_cast<unsigned>(std::min<std::int64_t>(count, std::max<std::int64_t>(1, n_tested / 8)));
    if (count <= 1) {
        work(0, n_tested);
        return;
    }
    std::vector<std::thread> pool;
    const std::int64_t chunk = (n_tested + count - 1) / count;
    for (unsigned k = 0; k < count; ++k) {
        const std::int64_t b = k * chunk, e = std::min(n_tested, b + chunk);
        if (b >= e) break;
        pool.emplace_back(work, b, e);
    }
    for (auto& th : pool) th.join();
}

}  // namespace apex

// Voir apex/isosurface.hpp pour le contrat d'equivalence (au bit pres) et
// l'ordre de sortie.

#include "apex/isosurface.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>

namespace apex {

namespace {

struct P3 {
    double x, y, z;
};

// `p_from + t * (p_to - p_from)`, t = (iso - v_from) / (v_to - v_from),
// composante par composante comme numpy (_interpolate_edge_batch).
inline P3 interp(const P3& pf, double vf, const P3& pt, double vt, double iso) {
    const double t = (iso - vf) / (vt - vf);
    return P3{pf.x + t * (pt.x - pf.x), pf.y + t * (pt.y - pf.y), pf.z + t * (pt.z - pf.z)};
}

// _orient_triangle_batch : normal = cross(p1 - p0, p2 - p0) ; dot avec la
// reference ; si dot < 0, l'ordre devient (p0, p2, p1).
inline void emit_oriented(const P3& p0, const P3& p1, const P3& p2, const P3& ref,
                          std::vector<double>& out) {
    const double ax = p1.x - p0.x, ay = p1.y - p0.y, az = p1.z - p0.z;
    const double bx = p2.x - p0.x, by = p2.y - p0.y, bz = p2.z - p0.z;
    // np.cross : (a1*b2 - a2*b1, a2*b0 - a0*b2, a0*b1 - a1*b0)
    const double nx = ay * bz - az * by;
    const double ny = az * bx - ax * bz;
    const double nz = ax * by - ay * bx;
    const double dot = nx * ref.x + ny * ref.y + nz * ref.z;
    const P3& q1 = dot < 0.0 ? p2 : p1;
    const P3& q2 = dot < 0.0 ? p1 : p2;
    out.insert(out.end(), {p0.x, p0.y, p0.z, q1.x, q1.y, q1.z, q2.x, q2.y, q2.z});
}

}  // namespace

IsosurfaceExtractor::IsosurfaceExtractor(const double* tet_vertices, std::int64_t vertex_count,
                                         const std::int64_t* tetrahedra, std::int64_t tet_count,
                                         const double* values)
    : vertices_(tet_vertices, tet_vertices + vertex_count * 3),
      tets_(tetrahedra, tetrahedra + tet_count * 4),
      values_(values, values + vertex_count),
      tet_min_(static_cast<std::size_t>(tet_count)),
      tet_max_(static_cast<std::size_t>(tet_count)) {
    for (std::int64_t t = 0; t < tet_count; ++t) {
        const std::int64_t* q = &tets_[static_cast<std::size_t>(t * 4)];
        double lo = values_[q[0]], hi = values_[q[0]];
        for (int k = 1; k < 4; ++k) {
            lo = std::min(lo, values_[q[k]]);
            hi = std::max(hi, values_[q[k]]);
        }
        tet_min_[t] = lo;
        tet_max_[t] = hi;
    }
}

IsosurfaceResult IsosurfaceExtractor::extract(double iso, double weld_tolerance) const {
    const std::int64_t T = tet_count();
    std::vector<std::int64_t> caps, quads;
    for (std::int64_t t = 0; t < T; ++t) {
        if (!(tet_min_[t] <= iso && iso <= tet_max_[t])) {
            continue;
        }
        const std::int64_t* q = &tets_[static_cast<std::size_t>(t * 4)];
        int n_above = 0;
        for (int k = 0; k < 4; ++k) {
            n_above += values_[q[k]] >= iso ? 1 : 0;
        }
        // champ constant (hi == lo) ou aucun changement de cote : rien
        if (tet_max_[t] == tet_min_[t] || n_above == 0 || n_above == 4) {
            continue;
        }
        (n_above == 2 ? quads : caps).push_back(t);
    }

    IsosurfaceResult res;
    std::vector<double> raw;
    raw.reserve((caps.size() * 3 + quads.size() * 6) * 3);
    res.source_tets.reserve(caps.size() + quads.size() * 2);

    auto load = [&](std::int64_t t, std::array<P3, 4>& p, std::array<double, 4>& v) {
        const std::int64_t* q = &tets_[static_cast<std::size_t>(t * 4)];
        for (int k = 0; k < 4; ++k) {
            const double* c = &vertices_[static_cast<std::size_t>(q[k] * 3)];
            p[k] = P3{c[0], c[1], c[2]};
            v[k] = values_[q[k]];
        }
    };
    auto reference = [](const std::array<P3, 4>& p, const std::array<double, 4>& v) {
        int hi = 0, lo = 0;  // premier extremum, comme argmax/argmin
        for (int k = 1; k < 4; ++k) {
            if (v[k] > v[hi]) hi = k;
            if (v[k] < v[lo]) lo = k;
        }
        return P3{p[hi].x - p[lo].x, p[hi].y - p[lo].y, p[hi].z - p[lo].z};
    };

    std::array<P3, 4> p;
    std::array<double, 4> v;
    for (std::int64_t t : caps) {
        load(t, p, v);
        bool above[4];
        int n_above = 0;
        for (int k = 0; k < 4; ++k) {
            above[k] = v[k] >= iso;
            n_above += above[k] ? 1 : 0;
        }
        const bool lone_is_above = n_above == 1;
        int lone = 0;
        while (above[lone] != lone_is_above) ++lone;
        P3 e[3];
        int j = 0;
        for (int k = 0; k < 4; ++k) {
            if (k == lone) continue;
            e[j++] = lone_is_above ? interp(p[k], v[k], p[lone], v[lone], iso)
                                   : interp(p[lone], v[lone], p[k], v[k], iso);
        }
        emit_oriented(e[0], e[1], e[2], reference(p, v), raw);
        res.source_tets.push_back(t);
    }
    for (std::int64_t t : quads) {
        load(t, p, v);
        int a[2], b[2], na = 0, nb = 0;
        for (int k = 0; k < 4; ++k) {
            if (v[k] >= iso) a[na++] = k; else b[nb++] = k;
        }
        const P3 a1b1 = interp(p[b[0]], v[b[0]], p[a[0]], v[a[0]], iso);
        const P3 a1b2 = interp(p[b[1]], v[b[1]], p[a[0]], v[a[0]], iso);
        const P3 a2b1 = interp(p[b[0]], v[b[0]], p[a[1]], v[a[1]], iso);
        const P3 a2b2 = interp(p[b[1]], v[b[1]], p[a[1]], v[a[1]], iso);
        const P3 ref = reference(p, v);
        emit_oriented(a1b1, a1b2, a2b2, ref, raw);
        emit_oriented(a1b1, a2b2, a2b1, ref, raw);
        res.source_tets.push_back(t);
        res.source_tets.push_back(t);
    }

    const std::size_t nv = raw.size() / 3;
    if (weld_tolerance <= 0.0 || nv == 0) {
        res.vertices = std::move(raw);
        res.faces.resize(nv);
        std::iota(res.faces.begin(), res.faces.end(), std::int64_t{0});
        return res;
    }

    // _weld_vertices : cle = round(v / tol) (arrondi au pair), tri stable
    // lexicographique, premiere occurrence representative.
    std::vector<std::array<std::int64_t, 3>> keys(nv);
    for (std::size_t i = 0; i < nv; ++i) {
        for (int c = 0; c < 3; ++c) {
            keys[i][c] = static_cast<std::int64_t>(std::nearbyint(raw[i * 3 + c] / weld_tolerance));
        }
    }
    std::vector<std::size_t> order(nv);
    std::iota(order.begin(), order.end(), std::size_t{0});
    std::stable_sort(order.begin(), order.end(),
                     [&](std::size_t x, std::size_t y) { return keys[x] < keys[y]; });
    std::vector<std::int64_t> inverse(nv);
    std::int64_t group = -1;
    for (std::size_t r = 0; r < nv; ++r) {
        const std::size_t i = order[r];
        if (r == 0 || keys[i] != keys[order[r - 1]]) {
            ++group;
            res.vertices.insert(res.vertices.end(), {raw[i * 3], raw[i * 3 + 1], raw[i * 3 + 2]});
        }
        inverse[i] = group;
    }
    res.faces = std::move(inverse);  // faces avant soudure = 0..nv-1, donc faces = inverse
    return res;
}

BoundaryLoops boundary_loops(const std::int64_t* faces, std::int64_t face_count, std::int64_t vertex_count) {
    BoundaryLoops out;
    std::vector<std::uint64_t> edges;
    edges.reserve(static_cast<std::size_t>(face_count) * 3);
    for (std::int64_t f = 0; f < face_count; ++f) {
        const std::int64_t a = faces[f * 3], b = faces[f * 3 + 1], c = faces[f * 3 + 2];
        if (a == b || b == c || a == c) {
            continue;  // triangle degenere : aucun bord (limite #60)
        }
        const std::int64_t pairs[3][2] = {{a, b}, {b, c}, {c, a}};
        for (const auto& e : pairs) {
            const auto lo = static_cast<std::uint64_t>(std::min(e[0], e[1]));
            const auto hi = static_cast<std::uint64_t>(std::max(e[0], e[1]));
            edges.push_back((lo << 32) | hi);
        }
    }
    if (edges.empty()) {
        out.ok = true;
        return out;
    }
    std::sort(edges.begin(), edges.end());
    std::vector<std::vector<std::int64_t>> adj(static_cast<std::size_t>(vertex_count));
    for (std::size_t i = 0; i < edges.size();) {
        std::size_t j = i;
        while (j < edges.size() && edges[j] == edges[i]) ++j;
        if (j - i == 1) {
            const auto lo = static_cast<std::int64_t>(edges[i] >> 32);
            const auto hi = static_cast<std::int64_t>(edges[i] & 0xffffffffULL);
            adj[lo].push_back(hi);
            adj[hi].push_back(lo);
        }
        i = j;
    }
    std::vector<char> visited(static_cast<std::size_t>(vertex_count), 0);
    for (std::int64_t s = 0; s < vertex_count; ++s) {
        const auto& ns = adj[s];
        if (ns.empty()) continue;
        if (ns.size() != 2) {
            return out;  // non-variete ou ouvert : repli trimesh
        }
    }
    for (std::int64_t s = 0; s < vertex_count; ++s) {
        if (adj[s].empty() || visited[s]) continue;
        std::vector<std::int64_t> loop{s};
        visited[s] = 1;
        std::int64_t prev = s;
        std::int64_t cur = std::min(adj[s][0], adj[s][1]);
        while (cur != s) {
            if (visited[cur]) {
                return BoundaryLoops{};  // incoherence : repli
            }
            visited[cur] = 1;
            loop.push_back(cur);
            const auto& n = adj[cur];
            const std::int64_t next = n[0] == prev ? n[1] : n[0];
            prev = cur;
            cur = next;
        }
        loop.push_back(s);
        out.loops.push_back(std::move(loop));
    }
    out.ok = true;
    return out;
}

}  // namespace apex

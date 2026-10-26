// Qualité du maillage découpé (Labelle-Shewchuk) sur une sphère et un cube : angles dièdres,
// bord fermé, écart de peau, volume. Seuils choisis sous les valeurs mesurées (marge de 25 %).
#include "tetranet/tetranet.hpp"
#include <cmath>
#include <cstdio>
#include <array>
#include <map>
#include <vector>

static tetranet::SurfaceMesh icosphere(int subdiv, double radius) {
    tetranet::SurfaceMesh s;
    const double t = (1.0 + std::sqrt(5.0)) / 2.0;
    std::vector<std::array<double, 3>> v = {{-1, t, 0}, {1, t, 0}, {-1, -t, 0}, {1, -t, 0}, {0, -1, t}, {0, 1, t},
                                            {0, -1, -t}, {0, 1, -t}, {t, 0, -1}, {t, 0, 1}, {-t, 0, -1}, {-t, 0, 1}};
    std::vector<std::array<int, 3>> f = {{0, 11, 5}, {0, 5, 1}, {0, 1, 7}, {0, 7, 10}, {0, 10, 11}, {1, 5, 9}, {5, 11, 4}, {11, 10, 2},
                                         {10, 7, 6}, {7, 1, 8}, {3, 9, 4}, {3, 4, 2}, {3, 2, 6}, {3, 6, 8}, {3, 8, 9}, {4, 9, 5},
                                         {2, 4, 11}, {6, 2, 10}, {8, 6, 7}, {9, 8, 1}};
    for (auto& p : v) { double n = std::sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]); for (auto& c : p) c /= n; }
    for (int it = 0; it < subdiv; ++it) {
        std::map<std::pair<int, int>, int> mid;
        auto midpoint = [&](int a, int b) {
            auto key = std::make_pair(std::min(a, b), std::max(a, b));
            auto found = mid.find(key);
            if (found != mid.end()) return found->second;
            std::array<double, 3> p{(v[a][0] + v[b][0]) / 2, (v[a][1] + v[b][1]) / 2, (v[a][2] + v[b][2]) / 2};
            double n = std::sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]); for (auto& c : p) c /= n;
            v.push_back(p); return mid[key] = int(v.size()) - 1;
        };
        std::vector<std::array<int, 3>> nf;
        for (auto& t3 : f) {
            int a = midpoint(t3[0], t3[1]), b = midpoint(t3[1], t3[2]), c = midpoint(t3[2], t3[0]);
            nf.push_back({t3[0], a, c}); nf.push_back({t3[1], b, a}); nf.push_back({t3[2], c, b}); nf.push_back({a, b, c});
        }
        f = nf;
    }
    for (auto& p : v) { s.vertices.push_back(p[0] * radius); s.vertices.push_back(p[1] * radius); s.vertices.push_back(p[2] * radius); }
    for (auto& t3 : f) { s.faces.push_back(t3[0]); s.faces.push_back(t3[1]); s.faces.push_back(t3[2]); }
    return s;
}

int main() {
    int bad = 0;
    auto check = [&](bool ok, const char* what, double value) { std::printf("%-46s %10.4f : %s\n", what, value, ok ? "ok" : "ECHEC"); if (!ok) ++bad; };
    tetranet::Params p; p.cell_size = 1.5;
    tetranet::TetMesh m;
    auto r = tetranet::mesh(icosphere(4, 20.0), p, m);
    check(r.min_dihedral_deg > 12.0, "sphère : angle dièdre minimal (deg)", r.min_dihedral_deg);
    check(r.dropped_degenerate == 0, "sphère : tétraèdres écartés", double(r.dropped_degenerate));
    check(r.nonmanifold_edges == 0, "sphère : arêtes de bord non fermées", double(r.nonmanifold_edges));
    check(r.skin_max_error < 0.06, "sphère : écart de peau maximal (mm)", r.skin_max_error);
    check(std::fabs(r.volume_ratio - 1.0) < 0.01, "sphère : volume / volume de la surface", r.volume_ratio);
    p.use_cuts = false;
    auto r0 = tetranet::mesh(icosphere(4, 20.0), p, m);
    check(r0.min_dihedral_deg < r.min_dihedral_deg, "sans découpes, la qualité est pire (garde)", r0.min_dihedral_deg);
    return bad;
}

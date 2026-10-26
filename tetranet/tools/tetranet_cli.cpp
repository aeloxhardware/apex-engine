// tetranet_cli model.stl --cell 0.5 [--alpha 0.4] [--tol 0.05] [--threads N] [--out prefix]
#include "tetranet/tetranet.hpp"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <vector>

static bool load_stl(const char* path, tetranet::SurfaceMesh& s) {
    std::ifstream f(path, std::ios::binary);
    if (!f) return false;
    std::vector<char> d((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
    if (d.size() >= 84) {
        uint32_t n; std::memcpy(&n, &d[80], 4);
        if (84 + 50ull * n == d.size()) {
            for (uint32_t i = 0; i < n; ++i) {
                const char* r = &d[84 + 50ull * i];
                for (int v = 0; v < 3; ++v) {
                    float c[3]; std::memcpy(c, r + 12 + 12 * v, 12);
                    for (int k = 0; k < 3; ++k) s.vertices.push_back(c[k]);
                    s.faces.push_back(int32_t(3 * i + v));
                }
            }
            return true;
        }
    }
    std::istringstream is(std::string(d.begin(), d.end()));
    std::string w; int c = 0;
    while (is >> w) if (w == "vertex") { double x, y, z; is >> x >> y >> z; s.vertices.insert(s.vertices.end(), {x, y, z}); s.faces.push_back(c++); }
    return c > 0;
}

int main(int argc, char** argv) {
    if (argc < 2) { std::fprintf(stderr, "usage: tetranet_cli model.stl --cell h [--alpha a] [--tol t] [--threads n] [--out prefix]\n"); return 2; }
    tetranet::Params p; std::string out;
    for (int i = 2; i + 1 < argc; i += 2) {
        std::string k = argv[i];
        if (k == "--cell") p.cell_size = std::atof(argv[i + 1]);
        else if (k == "--alpha") p.snap_alpha = std::atof(argv[i + 1]);
        else if (k == "--tol") p.skin_tolerance = std::atof(argv[i + 1]);
        else if (k == "--threads") p.threads = std::atoi(argv[i + 1]);
        else if (k == "--out") out = argv[i + 1];
    }
    tetranet::SurfaceMesh s; tetranet::TetMesh m;
    if (!load_stl(argv[1], s)) { std::fprintf(stderr, "lecture impossible : %s\n", argv[1]); return 2; }
    auto r = tetranet::mesh(s, p, m);
    std::printf("{\"status\":%d,\"reason\":\"%s\",\"seconds\":%.3f,\"t_bvh\":%.3f,\"t_classify\":%.3f,\"t_tets\":%.3f,\"t_checks\":%.3f,"
                "\"lattice_vertices\":%lld,\"band_vertices\":%lld,\"snapped_vertices\":%lld,\"tets\":%lld,\"dropped_degenerate\":%lld,"
                "\"sign_conflicts\":%lld,\"skin_max_error\":%.5f,\"skin_p99_error\":%.5f,\"volume_ratio\":%.5f,\"min_dihedral_deg\":%.3f,"
                "\"min_volume\":%.3e,\"boundary_faces\":%lld,\"nonmanifold_edges\":%lld}\n",
                int(r.status), r.reason.c_str(), r.seconds, r.t_bvh, r.t_classify, r.t_tets, r.t_checks,
                (long long)r.lattice_vertices, (long long)r.band_vertices, (long long)r.snapped_vertices, (long long)r.tets,
                (long long)r.dropped_degenerate, (long long)r.sign_conflicts, r.skin_max_error, r.skin_p99_error, r.volume_ratio,
                r.min_dihedral_deg, r.min_volume, (long long)r.boundary_faces, (long long)r.nonmanifold_edges);
    if (!out.empty() && !m.tetrahedra.empty()) {
        std::ofstream(out + ".nodes", std::ios::binary).write((const char*)m.vertices.data(), m.vertices.size() * 8);
        std::ofstream(out + ".tets", std::ios::binary).write((const char*)m.tetrahedra.data(), m.tetrahedra.size() * 4);
    }
    return int(r.status);
}

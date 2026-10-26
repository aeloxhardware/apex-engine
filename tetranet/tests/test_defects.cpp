// Un cube de 10 mm, intact puis abîmé : le volume du maillage doit rester celui du cube.
#include "tetranet/tetranet.hpp"
#include <cmath>
#include <cstdio>

static tetranet::SurfaceMesh cube() {
    tetranet::SurfaceMesh s;
    const double v[8][3] = {{0,0,0},{10,0,0},{10,10,0},{0,10,0},{0,0,10},{10,0,10},{10,10,10},{0,10,10}};
    for (auto& p : v) s.vertices.insert(s.vertices.end(), {p[0], p[1], p[2]});
    const int f[12][3] = {{0,2,1},{0,3,2},{4,5,6},{4,6,7},{0,1,5},{0,5,4},{1,2,6},{1,6,5},{2,3,7},{2,7,6},{3,0,4},{3,4,7}};
    for (auto& t : f) s.faces.insert(s.faces.end(), {t[0], t[1], t[2]});
    return s;
}

static double volume(const tetranet::TetMesh& m) {
    double vol = 0;
    for (size_t t = 0; t + 3 < m.tetrahedra.size(); t += 4) {
        const double* p[4];
        for (int k = 0; k < 4; ++k) p[k] = &m.vertices[3 * size_t(m.tetrahedra[t + k])];
        double a[3], b[3], c[3];
        for (int i = 0; i < 3; ++i) { a[i] = p[1][i] - p[0][i]; b[i] = p[2][i] - p[0][i]; c[i] = p[3][i] - p[0][i]; }
        vol += (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0])) / 6.0;
    }
    return vol;
}

static int check(const char* name, tetranet::SurfaceMesh s, double lo, double hi) {
    tetranet::Params p; p.cell_size = 1.0; p.snap_alpha = 0.5;
    tetranet::TetMesh m;
    auto r = tetranet::mesh(s, p, m);
    double v = volume(m);
    bool ok = v > lo && v < hi;
    std::printf("%-22s volume %.1f (attendu %.0f-%.0f) statut %d : %s\n", name, v, lo, hi, int(r.status), ok ? "ok" : "ECHEC");
    return ok ? 0 : 1;
}

int main() {
    int bad = 0;
    bad += check("cube intact", cube(), 930, 1070);
    { auto s = cube(); s.faces.resize(s.faces.size() - 3); bad += check("un triangle manquant", s, 930, 1070); }       // trou : moitié d'une face
    { auto s = cube(); s.faces.resize(s.faces.size() - 6); bad += check("une face manquante", s, 900, 1070); }          // trou : une face entière
    { auto s = cube(); std::swap(s.faces[6], s.faces[7]); bad += check("un triangle retourné", s, 930, 1070); }       // normale inversée
    { auto s = cube(); for (size_t i = 0; i + 2 < s.faces.size(); i += 3) std::swap(s.faces[i + 1], s.faces[i + 2]); bad += check("cube entier à l'envers", s, 930, 1070); }
    return bad;
}

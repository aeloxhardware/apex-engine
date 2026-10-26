// TetraNet — API publique (prototype). Millimètres, tableaux plats row-major.
#pragma once
#include <cstdint>
#include <string>
#include <vector>

namespace tetranet {

struct SurfaceMesh {
    std::vector<double> vertices;    // 3*N
    std::vector<int32_t> faces;      // 3*M
};

struct Params {
    double cell_size = 1.5;          // arête du réseau BCC (mm)
    double skin_tolerance = 0.05;    // écart toléré entre peau du maillage et surface (mm)
    double snap_alpha = 0.4;         // un sommet à moins de alpha*cell_size de la peau y est ramené
    double min_volume_ratio = 0.01;  // volume minimal d'un tétraèdre / volume nominal (cell^3/12)
    bool has_platform = false;
    double platform_z = 0.0;
    double platform_cell_size = 0.0; // 0 = pas de raffinement local (phase 5)
    uint32_t threads = 0;            // 0 = tous les cœurs
};

enum class Status { Ok, OutOfTolerance, InvalidInput, TooLarge };

struct TetMesh {
    std::vector<double> vertices;    // 3*P
    std::vector<int32_t> tetrahedra; // 4*T, orientés avec un volume positif
};

struct Report {
    Status status = Status::InvalidInput;
    std::string reason;
    // temps (s)
    double seconds = 0.0, t_bvh = 0.0, t_classify = 0.0, t_tets = 0.0, t_checks = 0.0;
    // taille
    int64_t lattice_vertices = 0, band_vertices = 0, snapped_vertices = 0;
    int64_t tets = 0, dropped_degenerate = 0, sign_conflicts = 0;
    // contrôles
    double skin_max_error = 0.0, skin_p99_error = 0.0;  // mm, maillage -> surface
    double volume_ratio = 0.0;       // volume du maillage / volume de la surface
    double min_dihedral_deg = 0.0;
    double min_volume = 0.0;         // mm3
    int64_t boundary_faces = 0, nonmanifold_edges = 0;
};

const char* version();
Report mesh(const SurfaceMesh& surface, const Params& params, TetMesh& out);

}  // namespace tetranet

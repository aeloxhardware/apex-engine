// TetraNet — API publique. Prototype : seul l'état du contrat est déclaré.
#pragma once
#include <cstdint>
#include <string>
#include <vector>

namespace tetranet {

// Surface d'entrée : tableaux plats, row-major, millimètres.
struct SurfaceMesh {
    std::vector<double> vertices;    // 3*N
    std::vector<int32_t> faces;      // 3*M
};

struct Params {
    double cell_size = 1.5;          // arête du réseau BCC (mm)
    double skin_tolerance = 0.05;    // écart toléré entre peau du maillage et STL (mm)
    bool has_platform = false;
    double platform_z = 0.0;
    double platform_cell_size = 0.0; // 0 = pas de raffinement local
    uint32_t threads = 0;            // 0 = tous les cœurs
};

enum class Status { Ok, OutOfTolerance, InvalidInput, TooLarge };

// Même contenu que core.mesh_types.TetMesh.
struct TetMesh {
    std::vector<double> vertices;    // 3*P
    std::vector<int32_t> tetrahedra; // 4*T
};

struct Report {
    Status status = Status::InvalidInput;
    std::string reason;              // lisible par l'utilisateur
    double seconds = 0.0;
    double skin_max_error = 0.0;     // mm
    double volume_ratio = 0.0;       // volume du maillage / volume de la surface
    double min_dihedral_deg = 0.0;
};

const char* version();
Report mesh(const SurfaceMesh& surface, const Params& params, TetMesh& out);

}  // namespace tetranet

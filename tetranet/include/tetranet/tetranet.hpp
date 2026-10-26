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
    double warn_conflict_ratio = 0.001;   // conflits de signe / sommets de la bande : avertissement
    double refuse_conflict_ratio = 0.05;  // au-delà : InconsistentSurface (surface trop cassée)
};

enum class Status { Ok, OutOfTolerance, InvalidInput, TooLarge, InconsistentSurface };

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
    // diagnostic de la surface d'entrée (nombres seulement ; le texte est écrit par le moteur / l'interface)
    int64_t surface_faces = 0, surface_edges = 0;
    int64_t surface_open_edges = 0;          // arêtes d'une seule face
    int64_t surface_holes = 0;               // boucles de bord distinctes : nombre de trous
    int64_t surface_nonmanifold_edges = 0;   // arêtes de plus de deux faces
    int64_t surface_flipped_edges = 0;       // deux faces qui parcourent l'arête dans le même sens : normales retournées
    int64_t surface_degenerate_faces = 0;    // aires nulles ou sommets confondus, écartées
    int64_t surface_duplicate_faces = 0;
    double conflict_ratio = 0.0;             // sign_conflicts / band_vertices
    bool warning = false;                    // défauts présents mais sous le seuil de refus
};

const char* version();
Report mesh(const SurfaceMesh& surface, const Params& params, TetMesh& out);

}  // namespace tetranet

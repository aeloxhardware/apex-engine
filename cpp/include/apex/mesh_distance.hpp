// Distance exacte point -> maillage de triangles (T9 du 2026-09-30).
//
// Cible mesuree : dans le creusement de P3 (lissage non planaire, recherche),
// trimesh.proximity.closest_point fait 67 % du temps (9,3 s sur 14 s pour 12
// couches du Benchy). Requetes independantes par point : c'est aussi la
// forme qui passe telle quelle sur GPU.
//
// CONTRAT D'INTERFACE (demande de l'utilisateur) :
// - le coeur ne depend JAMAIS de pybind11 : il recoit des vues brutes
//   (pointeur + taille, `MeshView`) ; aujourd'hui le binding les remplit
//   depuis numpy (sorties du code Python existant), demain les etapes C++
//   amont les rempliront directement, et le binding sera retire comme le
//   reste de l'echafaudage (limite #56) ;
// - trois backends derriere un seul point d'entree : CpuSerial (ORACLE),
//   CpuParallel (threads), Vulkan (GPU, optionnel). L'absence de GPU n'est
//   jamais une erreur : `resolve_backend` retombe sur le CPU et le dit.
//
// Equivalence : trimesh calcule le point le plus proche par d'autres
// operations flottantes -- ecart de quelques ULP, tolerance declaree.

#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "apex/types.hpp"

namespace apex {

struct MeshView {
    const double* vertices = nullptr;      // (n_vertices, 3), ligne par sommet
    std::size_t n_vertices = 0;
    const std::int64_t* faces = nullptr;   // (n_faces, 3)
    std::size_t n_faces = 0;
};

enum class DistanceBackend { CpuSerial = 0, CpuParallel = 1, Vulkan = 2 };

class MeshDistanceIndex {
public:
    explicit MeshDistanceIndex(const MeshView& mesh);

    // Distance de chaque point au maillage (et, si `closest` non nul, le point
    // le plus proche, (n,3)). `threads` <= 0 : nombre de coeurs.
    void distances(const double* points, std::size_t n, double* out, double* closest,
                   DistanceBackend backend, int threads = 0) const;

    std::size_t face_count() const { return tris_.size(); }

    // Backend reellement utilise pour une demande donnee (Vulkan -> CPU tant
    // que le backend GPU n'est pas construit ou qu'aucun peripherique n'existe).
    static DistanceBackend resolve_backend(DistanceBackend requested, std::string* reason = nullptr);

private:
    struct Node {
        AABB box;
        std::int32_t left = -1, right = -1, first = 0, count = 0;
    };
    std::int32_t build(std::int32_t lo, std::int32_t hi);
    double query(const Vec3& p, Vec3* best_point) const;
    void run_range(const double* points, std::size_t begin, std::size_t end, double* out, double* closest) const;

    std::vector<Triangle3> tris_;
    std::vector<AABB> boxes_;
    std::vector<Vec3> centroids_;
    std::vector<std::int32_t> order_;
    std::vector<Node> nodes_;
    std::int32_t root_ = -1;
};

// Point du triangle le plus proche de p (Ericson, Real-Time Collision
// Detection, 5.1.5) -- expose pour les tests.
Vec3 closest_point_on_triangle(const Vec3& p, const Triangle3& t);

}  // namespace apex

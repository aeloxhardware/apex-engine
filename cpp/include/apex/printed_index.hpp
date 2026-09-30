// Index spatial de la matiere deja imprimee - portage de
// `pipeline/step9_collision_sim.py::_PrintedGeometryIndex` (lignes 231-281).
//
// POURQUOI CE MODULE EST LE PREMIER PORTE : le profilage de `death.stl`
// (resolution 0.8) attribue a `candidates_near` 119,0 s de temps propre, plus
// 93 s de `set.add` (**1,41 milliard d'appels**) et 8,3 s de `dict.get`
// (46,7 M d'appels) qu'elle declenche - soit ~46 % du temps total du pipeline.
// Un seul portage sert QUATRE etapes : step9, step9b, step8b, step10 (+ step5).
//
// DEUX DEFAUTS DISTINCTS, mesures separement (voir plan) :
//   1. cout d'interpretation : liste de tuples materialisee par requete,
//      un `dict.get` par cellule, un `set.add` par (cellule x capsule) ;
//   2. defaut ALGORITHMIQUE : une capsule est inseree dans toutes les
//      cellules que son AABB touche, donc revisitee une fois par cellule a
//      chaque requete. C'est la limite #27/#30 de CLAUDE.md, jamais chiffree
//      jusqu'ici.
//
// Les deux modes de requete coexistent volontairement (`QueryMode`) : c'est ce
// qui permet de mesurer le gain du LANGAGE et le gain de l'ALGORITHME en deux
// chiffres separes, avec un seul build, sans reecrire le module entre les deux.
//
// DEUX SORTIES, parce que les consommateurs se repartissent en deux motifs
// (verifie en lisant les cinq sites d'appel, tous insensibles a l'ordre) :
//   - `candidates_indices` : step9 `_check_printed_part` et step8b
//     `_path_collides` bouclent sur les capsules une par une (filtre AABB puis
//     distance) - ils ont besoin des objets Python, donc des indices ;
//   - `candidates_endpoints` : step9b `_has_support_below`,
//     `_has_support_in_any_direction` et step5 `_layer_support_mask`
//     reconstruisent IMMEDIATEMENT deux tableaux (C,3) par
//     `np.array([s.start for s in candidates])`. C'est la source des 41 s de
//     `numpy.array` mesurees (817 653 requetes x ~432 candidats en mediane).
//     Leur rendre directement les tableaux supprime a la fois l'ensemble de
//     dedoublonnage, la liste d'objets et les deux compréhensions.

#pragma once

#include <cstdint>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "apex/types.hpp"

namespace apex {

// Coordonnees entieres d'une cellule de la grille uniforme.
struct Cell {
    int64_t x = 0;
    int64_t y = 0;
    int64_t z = 0;

    bool operator==(const Cell& other) const {
        return x == other.x && y == other.y && z == other.z;
    }
};

struct CellHash {
    std::size_t operator()(const Cell& c) const noexcept {
        // hash_combine facon boost - suffisant pour des coordonnees de
        // cellules reelles (petites, denses), pas une fonction cryptographique.
        std::size_t h = static_cast<std::size_t>(c.x) * 0x9E3779B97F4A7C15ULL;
        h ^= static_cast<std::size_t>(c.y) + 0x9E3779B97F4A7C15ULL + (h << 6) + (h >> 2);
        h ^= static_cast<std::size_t>(c.z) + 0x9E3779B97F4A7C15ULL + (h << 6) + (h >> 2);
        return h;
    }
};

enum class QueryMode {
    // Transcription fidele du Python : un ensemble de hachage par requete.
    // Sert a mesurer le gain du seul changement de langage.
    Faithful = 0,
    // Dedoublonnage par estampilles : un compteur de requete et un tableau
    // `uint32_t` par capsule. O(1) sans hachage ni allocation - c'est ce qui
    // supprime les 1,41 milliard d'operations de `set.add`.
    Stamped = 1,
};

// Compteurs CUMULES sur toutes les instances de PrintedGrid du processus.
// Les compteurs par instance ne suffisent pas a la mesure : l'index est cree
// et detruit a l'interieur de `step9b.run`, donc inaccessible a l'appelant
// une fois l'etape terminee. Non proteges contre la concurrence — a revoir
// quand une etape deviendra multi-thread (aucune ne l'est aujourd'hui).
struct GlobalStats {
    std::uint64_t visits = 0;
    std::uint64_t uniques = 0;
    std::uint64_t queries = 0;
    std::uint64_t appends = 0;
};

GlobalStats& global_stats();
void reset_global_stats();

class PrintedGrid {
public:
    PrintedGrid(double cell_size, double radius, QueryMode mode);

    // Insere une capsule d'indice implicite `size()`, en calculant son AABB
    // exactement comme `_segment_aabb`.
    void append(const Vec3& start, const Vec3& end);

    // Indices des capsules dont l'AABB touche une cellule touchee par `box`.
    // Sur-ensemble exact des vrais candidats, comme la version Python : le
    // filtrage fin reste a la charge de l'appelant.
    void candidates(const AABB& box, std::vector<int32_t>& out) const;

    const Vec3& start_of(int32_t index) const { return starts_[static_cast<std::size_t>(index)]; }
    const Vec3& end_of(int32_t index) const { return ends_[static_cast<std::size_t>(index)]; }

    std::int64_t size() const { return count_; }
    double radius() const { return radius_; }
    QueryMode mode() const { return mode_; }
    void set_mode(QueryMode mode) { mode_ = mode; }

    // Nombre total de VISITES (cellule x capsule) et de candidats UNIQUES
    // rendus depuis la construction. Leur rapport est le facteur de
    // duplication de la grille - la mesure qui dit si le defaut restant est le
    // langage ou l'algorithme. Instrumentation explicite, jamais un silence
    // sur ce que fait le mecanisme.
    std::uint64_t visit_count() const { return visits_; }
    std::uint64_t unique_count() const { return uniques_; }
    std::uint64_t query_count() const { return queries_; }
    void reset_stats() const { visits_ = 0; uniques_ = 0; queries_ = 0; }

private:
    void cell_range(const AABB& box, Cell& lo, Cell& hi) const;

    double cell_size_;
    double radius_;
    QueryMode mode_;
    std::int64_t count_ = 0;
    std::vector<Vec3> starts_;
    std::vector<Vec3> ends_;
    std::unordered_map<Cell, std::vector<int32_t>, CellHash> grid_;

    // Etat du mode Stamped et compteurs. `mutable` car `candidates()` est
    // logiquement const (elle n'altere pas l'index) mais doit incrementer
    // l'estampille et les compteurs.
    mutable std::vector<std::uint32_t> stamp_;
    mutable std::uint32_t current_stamp_ = 0;
    mutable std::uint64_t visits_ = 0;
    mutable std::uint64_t uniques_ = 0;
    mutable std::uint64_t queries_ = 0;
};

}  // namespace apex

// Controle de collision du lissage non planaire de la voie Bird
// (pipeline/bird/skin/print_order.py::order_and_check, partie 3).
//
// Derive de S4 Slicer (Joshua Bird, GPL-3.0) via la recherche d'Apex.
//
// Pour chaque point teste d'une zone de lissage, les obstacles sont : le corps deja
// depose et les zones de lissage gardees (couche PERMANENTE), les lignes deja tracees
// de la zone et le debut de la ligne en cours (couche TEMPORAIRE). La recherche
// reconstruisait un arbre k-d pour le debut de ligne A CHAQUE point teste ; ici,
// une grille spatiale incrementale (pas = reach + 1). Meme distance que cKDTree
// avec distance_upper_bound = reach + 1 (racine de la somme x -> y -> z, +inf
// au-dela), meme densification que dense() (a + (k/m)(b - a)), memes regles de
// plateau. Le positionnement de la tete (Rodrigues, rotations en ordre inverse de la
// liste) recoit cos/sin calcules par numpy : seule la reduction des produits
// scalaires peut differer de quelques ULP.
#pragma once

#include <cstdint>
#include <unordered_map>
#include <vector>

namespace apex {

class SkinCheckIndex {
public:
    explicit SkinCheckIndex(double cell);
    // Polylignes (points concatenes, `starts` : n_lines + 1 bornes), densifiees au pas `step`.
    void add_polylines(const double* pts, const std::int64_t* starts, std::int64_t n_lines, double step);
    std::int64_t size() const { return static_cast<std::int64_t>(perm_.xyz.size() / 3); }

    // Une zone. `cosv`/`sinv` (N, A) des angles d'axes par point ; `tested` (T,)
    // indices globaux croissants. Sorties (T,) : contact plateau, distance (inf si
    // >= reach + 1).
    void check_zone(const double* head, std::int64_t n_head, const double* axes, std::int64_t n_axes,
                    const double* pts, const std::int64_t* starts, std::int64_t n_lines,
                    const double* cosv, const double* sinv, const std::int64_t* tested, std::int64_t n_tested,
                    double plate_radius, double margin, double reach, double step,
                    std::uint8_t* plate_hit, double* distance) const;

    // Meme resultat que check_zone, IDENTIQUE AU BIT, reparti sur `threads` fils
    // (<= 0 : nombre de coeurs). La couche temporaire de check_zone ne fait que
    // croitre dans l'ordre d'insertion ; l'ensemble vu par un point teste est donc
    // un PREFIXE de cet ordre (bizarrerie des blocs de 256 points comprise : elle
    // allonge le prefixe d'un segment). Chaque point de la zone porte son rang ; un
    // point teste ne lit que les rangs < sa limite. Le minimum d'un meme ensemble
    // de distances ne depend pas de l'ordre de parcours.
    void check_zone_parallel(const double* head, std::int64_t n_head, const double* axes, std::int64_t n_axes,
                             const double* pts, const std::int64_t* starts, std::int64_t n_lines,
                             const double* cosv, const double* sinv, const std::int64_t* tested,
                             std::int64_t n_tested, double plate_radius, double margin, double reach, double step,
                             std::uint8_t* plate_hit, double* distance, int threads = 0) const;

private:
    struct Layer {
        std::vector<double> xyz;
        std::unordered_map<std::int64_t, std::vector<std::int64_t>> grid;
    };
    std::int64_t cell_key(std::int64_t i, std::int64_t j, std::int64_t k) const;
    void add_point(Layer& layer, double x, double y, double z) const;
    double nearest(const Layer& layer, double x, double y, double z, double bound) const;

    double cell_;
    Layer perm_;
};

}  // namespace apex

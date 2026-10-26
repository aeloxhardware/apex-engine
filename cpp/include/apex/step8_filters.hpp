// Filtres de l'etape 8 portes en natif : normale par point et filtre median
// des axes libres. Python de reference : `_apply_per_point_normals` et
// `_median_filter_free_axes` (pipeline/step8_machine_toolpath.py).
#pragma once

#include <cstdint>

namespace apex {

// Remplace la direction visee par la normale locale `local` (n*3) quand elle
// est finie, orientee comme `current` (produit scalaire > 0) ; `reversed_flag`
// (n, 0/1) inverse la normale locale (deposition_reversed). `desired` (n*3) est
// ecrit en place ; `changed` (optionnel) recoit le nombre de points remplaces.
void apply_per_point_normals(const double* local, const double* current,
                             const std::uint8_t* reversed_flag, std::int64_t n, double* desired,
                             std::int64_t* changed);

// Filtre median (fenetre IMPAIRE, bord replique comme scipy mode="nearest") de
// chaque colonne d'angle `axes[k]`, par polyligne : une polyligne s'arrete a un
// changement de (couche, segment) et a chaque point de trajet a vide ; seules
// celles de plus de 2 points et hors trajet sont filtrees. `angles` : n*axis_count,
// ecrit en place.
void median_filter_free_axes(double* angles, std::int64_t n, std::int64_t axis_count,
                             const std::int64_t* free_axes, std::int64_t free_count,
                             const std::int64_t* layer_ids, const std::int64_t* segments,
                             const std::uint8_t* travels, std::int64_t window);

}  // namespace apex

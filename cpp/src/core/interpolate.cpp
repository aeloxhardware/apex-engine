#include "apex/interpolate.hpp"

#include <cmath>

namespace apex {
namespace {

// Transcription de `_COINCIDENT_EPS` (step6_toolpath2d.py).
constexpr double kCoincidentEps = 1e-9;

}  // namespace

void interpolate_out_of_plane(const double* query, std::int64_t query_count,
                              const double* sources, const double* offsets,
                              std::int64_t source_count, double power, double* out) {
    if (query_count <= 0) {
        return;
    }
    if (source_count <= 0) {
        for (std::int64_t k = 0; k < query_count; ++k) {
            out[k] = 0.0;  // `np.zeros(len(query_2d))` cote Python
        }
        return;
    }

    for (std::int64_t k = 0; k < query_count; ++k) {
        const double qx = query[2 * k];
        const double qy = query[2 * k + 1];

        // Une SEULE passe : plus proche voisin, detection de coincidence et
        // les deux accumulations, dans l'ordre des sources (comme le numpy).
        double nearest_distance = 0.0;
        std::int64_t nearest_index = 0;
        bool any_coincident = false;
        double total = 0.0;
        double weighted = 0.0;

        for (std::int64_t m = 0; m < source_count; ++m) {
            const double dx = qx - sources[2 * m];
            const double dy = qy - sources[2 * m + 1];
            // `np.linalg.norm(..., axis=2)` : hypot euclidien, meme ordre
            // d'operations (somme des carres puis racine).
            const double distance = std::sqrt(dx * dx + dy * dy);

            // `np.argmin` rend le PREMIER minimum : comparaison stricte.
            if (m == 0 || distance < nearest_distance) {
                nearest_distance = distance;
                nearest_index = m;
            }

            if (distance < kCoincidentEps) {
                any_coincident = true;
                continue;  // poids 0.0, comme `np.where(coincident, 0.0, ...)`
            }
            const double weight = std::pow(distance, -power);
            total += weight;
            weighted += weight * offsets[m];
        }

        if (any_coincident || !(total > 0.0)) {
            out[k] = offsets[nearest_index];
        } else {
            out[k] = weighted / total;
        }
    }
}

}  // namespace apex

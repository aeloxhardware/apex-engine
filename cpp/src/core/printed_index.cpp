#include "apex/printed_index.hpp"

#include <algorithm>
#include <cmath>

namespace apex {

GlobalStats& global_stats() {
    static GlobalStats stats;
    return stats;
}

void reset_global_stats() { global_stats() = GlobalStats{}; }

PrintedGrid::PrintedGrid(double cell_size, double radius, QueryMode mode)
    : cell_size_(cell_size), radius_(radius), mode_(mode) {}

void PrintedGrid::cell_range(const AABB& box, Cell& lo, Cell& hi) const {
    // Transcription EXACTE de `_cell_coords_for_aabb` :
    //     lo = np.floor(aabb.min_corner / cell_size).astype(int)
    //     hi = np.floor(aabb.max_corner / cell_size).astype(int)
    // `astype(int)` tronque vers zero, mais s'applique a une valeur deja
    // plancher : le resultat est identique a floor. Reproduire une troncature
    // directe donnerait des cellules differentes pour les coordonnees
    // negatives - d'ou floor explicite, jamais un cast seul.
    lo.x = static_cast<int64_t>(std::floor(box.lo.x / cell_size_));
    lo.y = static_cast<int64_t>(std::floor(box.lo.y / cell_size_));
    lo.z = static_cast<int64_t>(std::floor(box.lo.z / cell_size_));
    hi.x = static_cast<int64_t>(std::floor(box.hi.x / cell_size_));
    hi.y = static_cast<int64_t>(std::floor(box.hi.y / cell_size_));
    hi.z = static_cast<int64_t>(std::floor(box.hi.z / cell_size_));
}

void PrintedGrid::append(const Vec3& start, const Vec3& end) {
    const int32_t index = static_cast<int32_t>(count_);
    ++count_;
    ++global_stats().appends;
    starts_.push_back(start);
    ends_.push_back(end);

    Cell lo, hi;
    cell_range(segment_aabb(start, end, radius_), lo, hi);

    for (int64_t x = lo.x; x <= hi.x; ++x) {
        for (int64_t y = lo.y; y <= hi.y; ++y) {
            for (int64_t z = lo.z; z <= hi.z; ++z) {
                grid_[Cell{x, y, z}].push_back(index);
            }
        }
    }
}

void PrintedGrid::candidates(const AABB& box, std::vector<int32_t>& out) const {
    out.clear();
    ++queries_;
    ++global_stats().queries;

    Cell lo, hi;
    cell_range(box, lo, hi);

    if (mode_ == QueryMode::Stamped) {
        // Une estampille par requete : un indice est neuf si son estampille
        // n'est pas celle de la requete courante. Pas de hachage, pas
        // d'allocation, pas de remise a zero du tableau entre requetes.
        if (stamp_.size() < static_cast<std::size_t>(count_)) {
            stamp_.resize(static_cast<std::size_t>(count_), 0u);
        }
        ++current_stamp_;
        if (current_stamp_ == 0u) {
            // Debordement du compteur 32 bits : seul cas ou il faut remettre
            // le tableau a zero, sinon d'anciennes estampilles passeraient
            // pour courantes. Survient apres ~4 milliards de requetes.
            std::fill(stamp_.begin(), stamp_.end(), 0u);
            current_stamp_ = 1u;
        }

        for (int64_t x = lo.x; x <= hi.x; ++x) {
            for (int64_t y = lo.y; y <= hi.y; ++y) {
                for (int64_t z = lo.z; z <= hi.z; ++z) {
                    const auto it = grid_.find(Cell{x, y, z});
                    if (it == grid_.end()) {
                        continue;
                    }
                    for (const int32_t index : it->second) {
                        ++visits_;
                        ++global_stats().visits;
                        std::uint32_t& slot = stamp_[static_cast<std::size_t>(index)];
                        if (slot != current_stamp_) {
                            slot = current_stamp_;
                            out.push_back(index);
                        }
                    }
                }
            }
        }
    } else {
        // Mode fidele : meme structure que le Python (un ensemble par
        // requete), pour isoler le gain du seul changement de langage.
        std::unordered_set<int32_t> seen;
        for (int64_t x = lo.x; x <= hi.x; ++x) {
            for (int64_t y = lo.y; y <= hi.y; ++y) {
                for (int64_t z = lo.z; z <= hi.z; ++z) {
                    const auto it = grid_.find(Cell{x, y, z});
                    if (it == grid_.end()) {
                        continue;
                    }
                    for (const int32_t index : it->second) {
                        ++visits_;
                        ++global_stats().visits;
                        seen.insert(index);
                    }
                }
            }
        }
        out.assign(seen.begin(), seen.end());
    }

    uniques_ += static_cast<std::uint64_t>(out.size());
    global_stats().uniques += static_cast<std::uint64_t>(out.size());
}

}  // namespace apex

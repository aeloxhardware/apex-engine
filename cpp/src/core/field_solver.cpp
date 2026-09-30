#include "apex/field_solver.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>

#include "apex/kdtree.hpp"

namespace apex {
namespace {

constexpr double kInf = std::numeric_limits<double>::infinity();

// AABB des points d'une source, ELARGIE du rayon d'influence.
//
// Un tetraedre hors de cette boite est a plus de `radius` de TOUTE la
// source : s'il sort par un axe, l'ecart sur ce seul axe depasse deja le
// rayon, donc la distance euclidienne aussi. La borne est donc exacte, pas
// heuristique - c'est ce qui autorise a ne jamais calculer ces distances.
AABB expanded_bounds(const double* points, std::int64_t count, double radius) {
    AABB box;
    box.lo = Vec3{kInf, kInf, kInf};
    box.hi = Vec3{-kInf, -kInf, -kInf};
    for (std::int64_t i = 0; i < count; ++i) {
        for (int a = 0; a < 3; ++a) {
            const double v = points[3 * i + a];
            (&box.lo.x)[a] = std::min((&box.lo.x)[a], v);
            (&box.hi.x)[a] = std::max((&box.hi.x)[a], v);
        }
    }
    for (int a = 0; a < 3; ++a) {
        (&box.lo.x)[a] -= radius;
        (&box.hi.x)[a] += radius;
    }
    return box;
}

inline bool inside(const AABB& box, const double* p) {
    return p[0] >= box.lo.x && p[0] <= box.hi.x && p[1] >= box.lo.y && p[1] <= box.hi.y &&
           p[2] >= box.lo.z && p[2] <= box.hi.z;
}

}  // namespace

void blend_tet_directions(const double* tet_centers, std::int64_t tet_count,
                          const BlendSources& sources, const BlendParams& params, double* out) {
    const std::int64_t region_count = sources.region_count;
    const std::int64_t source_count = region_count + (sources.has_platform ? 1 : 0);

    if (source_count == 0 || tet_count == 0) {
        // Aucune source exploitable : verticale partout, meme repli que le
        // Python (jamais une direction indeterminee).
        for (std::int64_t t = 0; t < tet_count; ++t) {
            out[3 * t + 0] = 0.0;
            out[3 * t + 1] = 0.0;
            out[3 * t + 2] = 1.0;
        }
        return;
    }

    // Matrice des distances (T, S), initialisee a l'infini : une source
    // ecartee garde donc un poids strictement nul (1/inf^p = 0), exactement
    // comme si sa distance reelle avait ete calculee puis fenetree a zero.
    std::vector<double> distances(static_cast<std::size_t>(tet_count * source_count), kInf);

    std::vector<KDTree> trees;
    std::vector<AABB> boxes;
    trees.reserve(static_cast<std::size_t>(region_count));
    boxes.reserve(static_cast<std::size_t>(region_count));

    std::vector<double> query_buffer;
    std::vector<std::int64_t> query_tets;
    std::vector<std::int64_t> nn_indices;
    std::vector<double> nn_distances;

    for (std::int64_t s = 0; s < region_count; ++s) {
        const std::int32_t begin = sources.offsets[s];
        const std::int32_t end = sources.offsets[s + 1];
        const std::int64_t point_count = end - begin;
        const double* points = sources.points + static_cast<std::size_t>(begin) * 3;

        trees.emplace_back(points, point_count);
        boxes.push_back(expanded_bounds(points, point_count, params.radius));

        // Regroupe les tetraedres retenus, puis UNE seule requete par lot :
        // franchir la frontiere par point couterait plus cher que la requete
        // elle-meme (meme lecon qu'au pas 4 sur `append_segments`).
        query_buffer.clear();
        query_tets.clear();
        for (std::int64_t t = 0; t < tet_count; ++t) {
            const double* c = tet_centers + 3 * t;
            if (inside(boxes.back(), c)) {
                query_tets.push_back(t);
                query_buffer.insert(query_buffer.end(), c, c + 3);
            }
        }
        if (query_tets.empty()) {
            continue;
        }
        nn_indices.resize(query_tets.size());
        nn_distances.resize(query_tets.size());
        trees.back().query_nn(query_buffer.data(), static_cast<std::int64_t>(query_tets.size()),
                              nn_indices.data(), nn_distances.data());
        for (std::size_t i = 0; i < query_tets.size(); ++i) {
            distances[static_cast<std::size_t>(query_tets[i] * source_count + s)] = nn_distances[i];
        }
    }

    if (sources.has_platform) {
        const std::int64_t s = source_count - 1;
        for (std::int64_t t = 0; t < tet_count; ++t) {
            distances[static_cast<std::size_t>(t * source_count + s)] =
                std::abs(tet_centers[3 * t + 2] - sources.platform_level);
        }
    }

    // SECONDE PASSE, sans elagage, pour les seuls tetraedres qu'aucune source
    // n'atteint. Leur direction vient du repli `nearest_source`, qui exige le
    // vrai argmin GLOBAL - lequel peut tomber sur une source ecartee. Sans
    // cette passe, l'elagage cesserait d'etre exact.
    for (std::int64_t t = 0; t < tet_count; ++t) {
        double* row = distances.data() + static_cast<std::size_t>(t * source_count);
        bool any_in_radius = false;
        for (std::int64_t s = 0; s < source_count; ++s) {
            if (row[s] <= params.radius) {
                any_in_radius = true;
                break;
            }
        }
        if (any_in_radius) {
            continue;
        }
        const double* c = tet_centers + 3 * t;
        for (std::int64_t s = 0; s < region_count; ++s) {
            if (row[s] != kInf) {
                continue;  // deja calculee au premier passage
            }
            std::int64_t index = 0;
            double d = 0.0;
            trees[static_cast<std::size_t>(s)].query_nn(c, 1, &index, &d);
            row[s] = d;
        }
    }

    // --- Ponderation : transcription directe du numpy ------------------------
    const double pi = std::numbers::pi_v<double>;
    std::vector<double> weights(static_cast<std::size_t>(source_count));

    for (std::int64_t t = 0; t < tet_count; ++t) {
        const double* row = distances.data() + static_cast<std::size_t>(t * source_count);

        std::int64_t nearest = 0;
        double nearest_distance = row[0];
        for (std::int64_t s = 1; s < source_count; ++s) {
            if (row[s] < nearest_distance) {
                nearest_distance = row[s];
                nearest = s;
            }
        }

        // `effective_radius` : le rayon nominal, sauf quand meme la source la
        // plus proche est hors rayon - la fenetre serait alors nulle partout
        // et diviserait par zero. Elle est alors recentree sur cette distance.
        const double effective_radius = std::max(params.radius, nearest_distance + params.epsilon);

        double weight_sum = 0.0;
        for (std::int64_t s = 0; s < source_count; ++s) {
            const double d = row[s];
            // `within_radius`, avec la meme correction que le Python : la
            // source la plus proche est reactivee de force quand aucune
            // n'entre dans le rayon.
            const bool in_radius = (d <= params.radius) || (s == nearest && nearest_distance > params.radius);
            double w = 0.0;
            if (in_radius) {
                const double ratio = std::min(d / effective_radius, 1.0);
                const double window = 0.5 * (1.0 + std::cos(pi * ratio));
                w = window / std::pow(d + params.epsilon, params.power);
            }
            weights[static_cast<std::size_t>(s)] = w;
            weight_sum += w;
        }

        double bx = 0.0, by = 0.0, bz = 0.0;
        for (std::int64_t s = 0; s < source_count; ++s) {
            const double w = weights[static_cast<std::size_t>(s)] / weight_sum;
            const double* dir = sources.directions + 3 * s;
            bx += w * dir[0];
            by += w * dir[1];
            bz += w * dir[2];
        }

        const double norm = std::sqrt(bx * bx + by * by + bz * bz);
        if (norm < params.min_norm) {
            // Cas degenere (sources quasi opposees a poids proches) : replier
            // sur la source dominante plutot que normaliser un vecteur
            // instable.
            std::int64_t best = 0;
            double best_weight = weights[0];
            for (std::int64_t s = 1; s < source_count; ++s) {
                if (weights[static_cast<std::size_t>(s)] > best_weight) {
                    best_weight = weights[static_cast<std::size_t>(s)];
                    best = s;
                }
            }
            const double* dir = sources.directions + 3 * best;
            out[3 * t + 0] = dir[0];
            out[3 * t + 1] = dir[1];
            out[3 * t + 2] = dir[2];
        } else {
            out[3 * t + 0] = bx / norm;
            out[3 * t + 1] = by / norm;
            out[3 * t + 2] = bz / norm;
        }
    }
}

void compute_shape_gradients(const double* tet_vertices, const std::int64_t* tetrahedra,
                             std::int64_t tet_count, double* shape_gradients, double* volumes) {
    for (std::int64_t t = 0; t < tet_count; ++t) {
        const std::int64_t* nodes = tetrahedra + 4 * t;
        const double* p0 = tet_vertices + 3 * nodes[0];
        const double* p1 = tet_vertices + 3 * nodes[1];
        const double* p2 = tet_vertices + 3 * nodes[2];
        const double* p3 = tet_vertices + 3 * nodes[3];

        // Jacobien [p1-p0, p2-p0, p3-p0] en COLONNES, comme le numpy.
        const double j00 = p1[0] - p0[0], j01 = p2[0] - p0[0], j02 = p3[0] - p0[0];
        const double j10 = p1[1] - p0[1], j11 = p2[1] - p0[1], j12 = p3[1] - p0[1];
        const double j20 = p1[2] - p0[2], j21 = p2[2] - p0[2], j22 = p3[2] - p0[2];

        const double det = j00 * (j11 * j22 - j12 * j21) - j01 * (j10 * j22 - j12 * j20) +
                           j02 * (j10 * j21 - j11 * j20);

        const double signed_volume = det / 6.0;
        double* dest = shape_gradients + 12 * t;

        // Meme seuil de degenerescence que le Python (`|volume| > 1e-14`) et
        // non un test contre zero : un tetraedre plat garde un determinant
        // minuscule mais non nul, dont l'inverse exploserait. Gradient ET
        // volume restent alors nuls, si bien que sa contribution disparait
        // d'elle-meme des assemblages.
        if (!(std::abs(signed_volume) > 1e-14)) {
            for (int i = 0; i < 12; ++i) {
                dest[i] = 0.0;
            }
            volumes[t] = 0.0;
            continue;
        }

        // Volume en valeur ABSOLUE : tetgen ne garantit aucune orientation
        // coherente d'un tetraedre a l'autre (meme raison que la limite #10).
        volumes[t] = std::abs(signed_volume);

        // Inverse de J : les gradients des fonctions de forme 1..3 sont les
        // LIGNES de J^-1 ; celui de la fonction 0 est l'oppose de leur somme
        // (les quatre fonctions somment a 1, donc leurs gradients a 0).
        const double inv_det = 1.0 / det;
        const double c00 = (j11 * j22 - j12 * j21) * inv_det;
        const double c01 = (j02 * j21 - j01 * j22) * inv_det;
        const double c02 = (j01 * j12 - j02 * j11) * inv_det;
        const double c10 = (j12 * j20 - j10 * j22) * inv_det;
        const double c11 = (j00 * j22 - j02 * j20) * inv_det;
        const double c12 = (j02 * j10 - j00 * j12) * inv_det;
        const double c20 = (j10 * j21 - j11 * j20) * inv_det;
        const double c21 = (j01 * j20 - j00 * j21) * inv_det;
        const double c22 = (j00 * j11 - j01 * j10) * inv_det;

        // dest est indexe [axe][fonction de forme] : 3 lignes de 4 colonnes.
        dest[0 * 4 + 1] = c00;
        dest[1 * 4 + 1] = c01;
        dest[2 * 4 + 1] = c02;
        dest[0 * 4 + 2] = c10;
        dest[1 * 4 + 2] = c11;
        dest[2 * 4 + 2] = c12;
        dest[0 * 4 + 3] = c20;
        dest[1 * 4 + 3] = c21;
        dest[2 * 4 + 3] = c22;
        for (int axis = 0; axis < 3; ++axis) {
            dest[axis * 4 + 0] =
                -(dest[axis * 4 + 1] + dest[axis * 4 + 2] + dest[axis * 4 + 3]);
        }
    }
}



// ---------------------------------------------------------------------------
// Assemblage FEM direct en CSR
// ---------------------------------------------------------------------------

namespace {

// Position de `col` dans la ligne triee [begin, end) du motif.
//
// Toujours trouvee par construction : le motif contient, par definition,
// toute paire de noeuds partageant un tetraedre, et l'assemblage ne visite
// que ces paires-la. Une recherche infructueuse signalerait un motif
// desynchronise de la connectivite - d'ou l'absence de branche de repli
// silencieuse.
inline std::int32_t find_column(const std::int32_t* indices, std::int32_t begin, std::int32_t end,
                                std::int32_t col) {
    const std::int32_t* first = indices + begin;
    const std::int32_t* last = indices + end;
    const std::int32_t* it = std::lower_bound(first, last, col);
    return static_cast<std::int32_t>(it - indices);
}

}  // namespace

SparsityPattern build_fem_sparsity(const std::int64_t* tetrahedra, std::int64_t tet_count,
                                   std::int64_t node_count) {
    SparsityPattern pattern;
    pattern.indptr.assign(static_cast<std::size_t>(node_count) + 1, 0);
    if (tet_count <= 0 || node_count <= 0) {
        return pattern;
    }

    // Passe 1 : chaque tetraedre apporte 4 entrees a chacune de ses 4 lignes.
    // On compte AVANT d'allouer plutot que d'agrandir 87 836 vecteurs.
    std::vector<std::int32_t> counts(static_cast<std::size_t>(node_count), 0);
    for (std::int64_t t = 0; t < tet_count; ++t) {
        const std::int64_t* nodes = tetrahedra + 4 * t;
        for (int i = 0; i < 4; ++i) {
            counts[static_cast<std::size_t>(nodes[i])] += 4;
        }
    }

    std::vector<std::int32_t> offsets(static_cast<std::size_t>(node_count) + 1, 0);
    for (std::int64_t n = 0; n < node_count; ++n) {
        offsets[static_cast<std::size_t>(n) + 1] =
            offsets[static_cast<std::size_t>(n)] + counts[static_cast<std::size_t>(n)];
    }

    // Passe 2 : remplissage brut, avec doublons.
    std::vector<std::int32_t> raw(static_cast<std::size_t>(offsets.back()));
    std::vector<std::int32_t> cursor(offsets.begin(), offsets.end() - 1);
    for (std::int64_t t = 0; t < tet_count; ++t) {
        const std::int64_t* nodes = tetrahedra + 4 * t;
        for (int i = 0; i < 4; ++i) {
            std::int32_t& at = cursor[static_cast<std::size_t>(nodes[i])];
            for (int j = 0; j < 4; ++j) {
                raw[static_cast<std::size_t>(at++)] = static_cast<std::int32_t>(nodes[j]);
            }
        }
    }

    // Passe 3 : tri + deduplication LIGNE PAR LIGNE. Les lignes comptent
    // quelques dizaines d'entrees, si bien que N petits tris coutent bien
    // moins qu'un tri global des 6,2 M paires - c'est precisement ce que fait
    // `coo_matrix.tocsr()`, et ce que ce chemin evite.
    pattern.indices.reserve(raw.size() / 2);
    for (std::int64_t n = 0; n < node_count; ++n) {
        const std::int32_t begin = offsets[static_cast<std::size_t>(n)];
        const std::int32_t end = offsets[static_cast<std::size_t>(n) + 1];
        std::sort(raw.begin() + begin, raw.begin() + end);
        const auto last = std::unique(raw.begin() + begin, raw.begin() + end);
        pattern.indices.insert(pattern.indices.end(), raw.begin() + begin, last);
        pattern.indptr[static_cast<std::size_t>(n) + 1] =
            static_cast<std::int32_t>(pattern.indices.size());
    }
    return pattern;
}

void assemble_gram_matrix(const std::int64_t* tetrahedra, const double* shape_gradients,
                          const double* scale, const std::int32_t* active,
                          std::int64_t active_count, std::int64_t tet_count,
                          const SparsityPattern& pattern, double* data) {
    const std::size_t nnz = pattern.indices.size();
    std::fill(data, data + nnz, 0.0);

    const std::int64_t count = (active != nullptr) ? active_count : tet_count;
    for (std::int64_t k = 0; k < count; ++k) {
        const std::int64_t t = (active != nullptr) ? active[k] : k;
        const double w = scale[k];
        if (w == 0.0) {
            continue;
        }
        const std::int64_t* nodes = tetrahedra + 4 * t;
        // G est (3,4) : gradient de chaque fonction de forme selon chaque axe.
        const double* g = shape_gradients + 12 * (active != nullptr ? k : t);

        for (int i = 0; i < 4; ++i) {
            const std::int32_t row = static_cast<std::int32_t>(nodes[i]);
            const std::int32_t begin = pattern.indptr[static_cast<std::size_t>(row)];
            const std::int32_t end = pattern.indptr[static_cast<std::size_t>(row) + 1];
            for (int j = 0; j < 4; ++j) {
                // local[i][j] = somme sur les 3 axes de G[axe][i] * G[axe][j],
                // dans l'ordre des axes - meme expression que l'einsum
                // "tki,tkj->tij" du Python.
                const double v = g[0 * 4 + i] * g[0 * 4 + j] + g[1 * 4 + i] * g[1 * 4 + j] +
                                 g[2 * 4 + i] * g[2 * 4 + j];
                const std::int32_t slot =
                    find_column(pattern.indices.data(), begin, end,
                                static_cast<std::int32_t>(nodes[j]));
                data[static_cast<std::size_t>(slot)] += w * v;
            }
        }
    }
}

void assemble_gradient_rhs(const std::int64_t* tetrahedra, const double* shape_gradients,
                           const double* scale, const double* directions,
                           const std::int32_t* active, std::int64_t active_count,
                           std::int64_t node_count, double* b) {
    std::fill(b, b + node_count, 0.0);
    for (std::int64_t k = 0; k < active_count; ++k) {
        const std::int64_t t = active[k];
        const double w = scale[k];
        const std::int64_t* nodes = tetrahedra + 4 * t;
        const double* g = shape_gradients + 12 * k;
        const double* d = directions + 3 * k;
        for (int i = 0; i < 4; ++i) {
            // einsum "tki,tk->ti" : somme sur les axes de G[axe][i] * d[axe].
            const double v = g[0 * 4 + i] * d[0] + g[1 * 4 + i] * d[1] + g[2 * 4 + i] * d[2];
            b[nodes[i]] += w * v;
        }
    }
}

}  // namespace apex

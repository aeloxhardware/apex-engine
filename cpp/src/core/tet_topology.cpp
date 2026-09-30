#include "apex/tet_topology.hpp"

#include <algorithm>
#include <array>

namespace apex {
namespace {

// Les 4 faces locales d'un tetraedre, dans l'ORDRE du Python - il determine
// l'ordre d'insertion, donc celui des faces frontiere rendues.
constexpr int kLocalFaces[4][3] = {{0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {1, 2, 3}};

struct FaceRecord {
    int32_t n0, n1, n2;  // triplet trie
    int32_t tet;
    int32_t seq;  // rang d'apparition, pour restituer l'ordre d'insertion
};

inline void sort3(int32_t& a, int32_t& b, int32_t& c) {
    if (a > b) std::swap(a, b);
    if (b > c) std::swap(b, c);
    if (a > b) std::swap(a, b);
}

}  // namespace

TetFaceTopology compute_tet_face_topology(const int64_t* tetrahedra, int64_t tet_count) {
    TetFaceTopology result;
    if (tet_count <= 0) {
        return result;
    }

    std::vector<FaceRecord> faces;
    faces.reserve(static_cast<std::size_t>(tet_count) * 4);

    int32_t seq = 0;
    for (int64_t t = 0; t < tet_count; ++t) {
        const int64_t* tet = tetrahedra + 4 * t;
        for (const auto& local : kLocalFaces) {
            int32_t a = static_cast<int32_t>(tet[local[0]]);
            int32_t b = static_cast<int32_t>(tet[local[1]]);
            int32_t c = static_cast<int32_t>(tet[local[2]]);
            sort3(a, b, c);
            faces.push_back(FaceRecord{a, b, c, static_cast<int32_t>(t), seq++});
        }
    }

    // Tri par cle de face, `seq` croissant a cle egale : le premier element
    // d'un groupe porte donc le rang de PREMIERE apparition de cette face.
    std::sort(faces.begin(), faces.end(), [](const FaceRecord& x, const FaceRecord& y) {
        if (x.n0 != y.n0) return x.n0 < y.n0;
        if (x.n1 != y.n1) return x.n1 < y.n1;
        if (x.n2 != y.n2) return x.n2 < y.n2;
        return x.seq < y.seq;
    });

    // Faces frontiere collectees avec leur rang, puis remises dans l'ordre
    // d'insertion - c'est l'ordre d'iteration du dict Python, dont depend
    // l'ordre de la requete kdtree en aval.
    std::vector<std::pair<int32_t, std::size_t>> boundary_by_seq;

    std::size_t i = 0;
    while (i < faces.size()) {
        std::size_t j = i + 1;
        while (j < faces.size() && faces[j].n0 == faces[i].n0 && faces[j].n1 == faces[i].n1 &&
               faces[j].n2 == faces[i].n2) {
            ++j;
        }
        const std::size_t owners = j - i;
        if (owners == 1) {
            boundary_by_seq.emplace_back(faces[i].seq, i);
        } else if (owners == 2) {
            // Une face interne est partagee par exactement deux tetraedres.
            // Un groupe de taille > 2 signalerait un maillage non-manifold :
            // le Python l'ignorerait aussi (il ne teste que len(owners) == 2).
            result.adjacency_pairs.push_back(faces[i].tet);
            result.adjacency_pairs.push_back(faces[i + 1].tet);
        }
        i = j;
    }

    std::sort(boundary_by_seq.begin(), boundary_by_seq.end());

    result.boundary_face_nodes.reserve(boundary_by_seq.size() * 3);
    result.boundary_face_tets.reserve(boundary_by_seq.size());
    for (const auto& [face_seq, index] : boundary_by_seq) {
        (void)face_seq;
        const FaceRecord& record = faces[index];
        result.boundary_face_nodes.push_back(record.n0);
        result.boundary_face_nodes.push_back(record.n1);
        result.boundary_face_nodes.push_back(record.n2);
        result.boundary_face_tets.push_back(record.tet);
    }
    return result;
}

std::vector<int32_t> compute_unique_edges(const int64_t* tetrahedra, int64_t tet_count) {
    std::vector<int32_t> edges;
    if (tet_count <= 0) {
        return edges;
    }

    // Chaque arete est empaquetee dans un uint64 (min << 32 | max) : le tri
    // devient un tri d'entiers, et le dedoublonnage un `std::unique`. Evite
    // le hachage de tuples qui coutait 1,35 s de `set.add` cote Python.
    std::vector<uint64_t> packed;
    packed.reserve(static_cast<std::size_t>(tet_count) * 6);

    for (int64_t t = 0; t < tet_count; ++t) {
        const int64_t* tet = tetrahedra + 4 * t;
        for (int i = 0; i < 4; ++i) {
            for (int j = i + 1; j < 4; ++j) {
                const auto a = static_cast<uint32_t>(tet[i]);
                const auto b = static_cast<uint32_t>(tet[j]);
                const uint32_t lo = a < b ? a : b;
                const uint32_t hi = a < b ? b : a;
                packed.push_back((static_cast<uint64_t>(lo) << 32) | hi);
            }
        }
    }

    std::sort(packed.begin(), packed.end());
    packed.erase(std::unique(packed.begin(), packed.end()), packed.end());

    edges.reserve(packed.size() * 2);
    for (const uint64_t key : packed) {
        edges.push_back(static_cast<int32_t>(key >> 32));
        edges.push_back(static_cast<int32_t>(key & 0xFFFFFFFFull));
    }
    return edges;
}

}  // namespace apex

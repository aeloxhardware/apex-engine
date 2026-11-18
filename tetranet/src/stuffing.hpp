// TetraNet — réseau BCC rempli jusqu'à la peau (« isosurface stuffing », Labelle et Shewchuk 2007) : phases internes.
#pragma once
#include "internal.hpp"

namespace tetranet::detail {

struct Lattice {   // deux grilles décalées : sommets A (coins des cellules) et B (centres)
    double h = 0;
    V3 lo;
    int nx = 0, ny = 0, nz = 0;
    int64_t NA = 0, NB = 0, N = 0;
    int64_t idA(int i, int j, int k) const { return (int64_t(k) * (ny + 1) + j) * (nx + 1) + i; }
    int64_t idB(int i, int j, int k) const { return NA + (int64_t(k) * ny + j) * nx + i; }
};

struct VertexState {   // un élément par nœud du réseau
    std::vector<V3> pos, pos0, qpt;   // position (ramenée), position d'origine, point de peau le plus proche
    std::vector<int8_t> sgn, isband;  // sgn : +1 dehors, -1 dedans, 0 sur la peau
    std::vector<float> dist, wval;    // distance à la peau (bande), nombre d'enroulement
};

struct CutVertices {   // sommets créés sur les arêtes qui traversent la peau
    std::vector<V3> pos, ref;                     // position réelle (sur la peau) et position linéaire de référence
    std::unordered_map<uint64_t, int64_t> id;     // arête (u,v) -> indice global (N + rang)
};

struct Stuffing {   // état partagé des phases
    const Params& prm;
    unsigned threads;
    Progress& prog;
    Report& rep;
    const Surface& surf;
    const Bvh& bvh;
    Winding wn;
    std::unique_ptr<FeatureEdges> features;
    double feature_radius = 0;
    Lattice lat;
    VertexState st;
    CutVertices cuts;
};

bool make_lattice(Stuffing& c);                    // faux (rapport TooLarge) si le réseau est trop grand
void classify_vertices(Stuffing& c);               // bande, point de peau, dedans/dehors
void snap_to_skin(Stuffing& c);                    // ramène les sommets proches ; coupe les arêtes restantes
void emit_tetrahedra(Stuffing& c, TetMesh& out);   // tétraèdres (découpés) puis sommets compactés

}  // namespace tetranet::detail

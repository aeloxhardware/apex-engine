// TetraNet — contrôles du maillage produit : volume, angles, bord, écart de peau, statut.
#pragma once
#include "internal.hpp"

namespace tetranet::detail {

struct CheckContext {
    const Params& prm;
    unsigned threads;
    Progress& prog;
    Report& rep;
    const Surface& surf;
    const Bvh& bvh;
    double h;
    const Winding* wn;   // nul sans nombre d'enroulement : toute la surface est alors comptée comme peau extérieure
};

void check_mesh(const CheckContext& c, const TetMesh& out);   // remplit les contrôles du rapport (phase Verify)
void final_status(const Params& prm, Report& rep);            // statut et raison, à partir du rapport rempli

}  // namespace tetranet::detail

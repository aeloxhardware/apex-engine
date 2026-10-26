// Methode de la chaleur (Crane et al.) : distance lisse au plateau, champ
// generateur des couches du chemin par defaut (`use_heat_geodesic_field`,
// limite #79). Port de `step4_field_solver._heat_method_distance`.
//
// Trois temps, comme le Python : diffusion `(M + t L) u = source`, flux
// normalise `X = -grad(u)/|grad(u)|` par tetraedre (projete eventuellement
// sur le demi-plan radial atteignable, limite #64), puis Poisson
// `(L + ancrage) phi = div(X)`.
//
// EQUIVALENCE : tolerance, pas le bit (meme regime que le reste de step4 :
// LDLT d'Eigen contre SuperLU sous le seuil direct, CG diagonal contre scipy
// au-dessus). Le repli « 10x plus d'iterations » du Python est transcrit.

#pragma once

#include <cstdint>
#include <vector>

namespace apex {

struct HeatParams {
    double time_factor = 100.0;         // heat_method_time_factor
    double anchor_weight = 1e6;         // platform_anchor_weight
    double tolerance = 1e-8;            // convergence_tolerance
    int max_iterations = 2000;          // max_solver_iterations
    std::int64_t direct_threshold = 30000;
    bool project_reachable = false;
    double origin[3] = {0, 0, 0};       // radial_tilt_origin
    double up[3] = {0, 0, 1};
};

struct HeatResult {
    std::vector<double> phi;  // une valeur par noeud
    bool converged = false;   // faux => l'appelant leve FieldSolverError
};

HeatResult heat_method_distance(const double* tet_vertices, const std::int64_t* tetrahedra,
                                std::int64_t tet_count, std::int64_t node_count,
                                const double* shape_gradients, const double* volumes,
                                const std::int64_t* source_nodes, std::int64_t source_count,
                                const HeatParams& params);

}  // namespace apex

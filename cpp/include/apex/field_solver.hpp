// Noyau numerique de step4 (champ scalaire FEM) - assemblage et resolution.
//
// POURQUOI CE MODULE DIFFERE DE TOUS LES PRECEDENTS : c'est le premier
// portage dont l'equivalence N'EST PAS bit-a-bit, et c'est une decision
// explicite de l'utilisateur, pas un accident.
//
// Les pas 1-6 reproduisaient exactement les flottants du Python, ce qui
// garantissait que la metrique step9b ne bougeait jamais. Ici, deux choses
// l'interdisent :
//   - le solveur : `scipy.sparse.linalg.cg` et le gradient conjugue d'Eigen
//     resolvent le meme systeme par le meme algorithme, mais leurs produits
//     scalaires n'accumulent pas dans le meme ordre ; a tolerance 1e-8, les
//     deux champs concordent a ~1e-8 relatif, jamais au bit ;
//   - les produits matrice-vecteur denses (`weights @ directions`) passent
//     par BLAS cote numpy, dont l'ordre de sommation n'est pas specifie.
//
// L'oracle de step4 compare donc avec une TOLERANCE explicite, documentee
// test par test, au lieu de `array_equal`. Ce n'est pas un relachement de
// discipline : c'est la reconnaissance qu'exiger le bit ici reviendrait a
// reimplementer BLAS et scipy plutot qu'a porter le moteur.
//
// PERIMETRE : le CHEMIN PAR DEFAUT du pipeline (celui que `main.py`
// exerce). Les mecanismes opt-in de step4 - directions geodesiques
// (limite #52), retroaction de surplomb (#53), lissage vectoriel iteratif
// (#55), modulation du lissage, diffusion des directions - restent
// integralement en Python et ne sont pas touches : ils sont desactives par
// defaut, et les porter d'avance figerait un contrat pour du code que
// personne n'execute.

#pragma once

#include <cstdint>
#include <vector>

#include "apex/types.hpp"

namespace apex {

// ---------------------------------------------------------------------------
// Melange spatial des directions de gradient (`_blend_tet_directions`)
// ---------------------------------------------------------------------------
//
// POSTE LE PLUS COUTEUX DE STEP4, ET POUR UNE RAISON ALGORITHMIQUE, PAS
// LINGUISTIQUE. Le Python interroge un k-d tree par REGION avec la totalite
// des centres de tetraedres : 22 regions x 387 978 tetraedres = 8,5 millions
// de distances sur death.stl @ 0.8, soit 2,31 s des 2,68 s de la fonction -
// et ce k-d tree est DEJA natif depuis le pas 5.
//
// Or la fenetre de Hann annule toute source au-dela de `radius` : mesure
// directe, seules **21,4 %** de ces paires sont dans le rayon, et **69,2 %**
// sont ecartables par un simple test d'AABB elargie. Transcrire la fonction
// telle quelle n'aurait donc rien donne - c'est exactement la lecon du pas 3
// (« avant de porter une fonction DEJA vectorisee, demander ce qu'elle
// calcule EN TROP »).
//
// L'elagage est fait ici, dans le portage, parce qu'il est trivial en C++
// (un test de boite avant de descendre l'arbre) et couteux a exprimer en
// numpy vectorise sans reintroduire une boucle Python.
//
// EXACTITUDE DE L'ELAGAGE : une source ecartee a forcement une distance
// superieure a `radius`, donc une fenetre nulle, donc un poids nul - elle
// n'aurait rien apporte a la moyenne ponderee. Le seul usage restant d'une
// distance hors rayon est le REPLI `nearest_source`, emprunte par les
// tetraedres qu'AUCUNE source n'atteint ; ceux-la sont donc repris dans une
// seconde passe SANS elagage, qui retrouve l'argmin exact. Rien n'est
// approxime : on ne calcule que ce dont le resultat depend.

struct BlendSources {
    // Points de chaque source-region, concatenes ; `offsets` decoupe le
    // tableau en S plages (CSR), une par region.
    const double* points = nullptr;     // P * 3
    const std::int32_t* offsets = nullptr;  // S + 1
    std::int64_t region_count = 0;
    // Direction cible de chaque source, regions PUIS plateau si present.
    const double* directions = nullptr;  // (region_count + has_platform) * 3
    bool has_platform = false;
    double platform_level = 0.0;
};

struct BlendParams {
    double radius = 0.0;
    double power = 2.0;
    double min_norm = 0.1;
    double epsilon = 0.0;  // 1e-6 * etendue max du maillage, calcule cote Python
};

// Ecrit `tet_count * 3` doubles dans `out`.
void blend_tet_directions(const double* tet_centers, std::int64_t tet_count,
                          const BlendSources& sources, const BlendParams& params, double* out);

// ---------------------------------------------------------------------------
// Gradients des fonctions de forme P1 et volumes (`_compute_shape_gradients`)
// ---------------------------------------------------------------------------
//
// `shape_gradients` : T * 3 * 4, gradient constant de chaque fonction de
// forme sur chaque tetraedre. `volumes` : T, valeur ABSOLUE du volume signe
// (tetgen ne garantit aucune orientation coherente - meme raison que la
// limite #10 pour marching tetrahedra).
void compute_shape_gradients(const double* tet_vertices, const std::int64_t* tetrahedra,
                             std::int64_t tet_count, double* shape_gradients, double* volumes);

// ---------------------------------------------------------------------------
// Assemblage FEM direct en CSR
// ---------------------------------------------------------------------------
//
// POURQUOI CE N'EST PAS UNE SIMPLE TRANSCRIPTION. Le Python assemble par
// `coo_matrix((data, (rows, cols))).tocsr()` : 16 entrees par tetraedre,
// soit 6,2 M sur death.stl @ 0.8, que scipy doit TRIER puis sommer par
// doublons - et il le refait INTEGRALEMENT pour chacune des deux matrices
// (rigidite et terme de direction), alors qu'elles ont exactement le MEME
// motif de sparsite : celui de l'adjacence noeud-noeud du maillage, qui ne
// depend que de la topologie, jamais des poids.
//
// Ici le motif est construit UNE FOIS, puis chaque matrice n'a plus qu'a
// accumuler ses valeurs dedans par recherche dichotomique. Plus de tri de
// 6,2 M entrees, plus de passe de deduplication, et un seul motif pour deux
// matrices.
//
// EQUIVALENCE : l'ordre d'accumulation differe de celui de
// `sum_duplicates`, et l'addition flottante n'est pas associative - les
// coefficients concordent donc a ~1e-15 relatif, pas au bit. C'est
// coherent avec la tolerance deja actee pour step4 (voir l'en-tete).

struct SparsityPattern {
    std::vector<std::int32_t> indptr;   // node_count + 1
    std::vector<std::int32_t> indices;  // nnz, TRIES a l'interieur de chaque ligne
};

// Motif CSR de l'assemblage FEM : une entree (i, j) pour chaque paire de
// noeuds partageant un tetraedre, diagonale comprise.
SparsityPattern build_fem_sparsity(const std::int64_t* tetrahedra, std::int64_t tet_count,
                                   std::int64_t node_count);

// Accumule `sum_e scale[e] * G_e^T G_e` dans `data` (deja dimensionne a
// `pattern.indices.size()`, remis a zero par la fonction).
// `active` a nullptr = tous les tetraedres.
void assemble_gram_matrix(const std::int64_t* tetrahedra, const double* shape_gradients,
                          const double* scale, const std::int32_t* active,
                          std::int64_t active_count, std::int64_t tet_count,
                          const SparsityPattern& pattern, double* data);

// Accumule `sum_e scale[e] * G_e^T d_e` dans `b` (remis a zero par la
// fonction). `directions` est indexe COMME `active` (une direction par
// tetraedre actif), deja normalisee par l'appelant.
void assemble_gradient_rhs(const std::int64_t* tetrahedra, const double* shape_gradients,
                           const double* scale, const double* directions,
                           const std::int32_t* active, std::int64_t active_count,
                           std::int64_t node_count, double* b);

// ---------------------------------------------------------------------------
// Resolution du systeme lineaire (Eigen)
// ---------------------------------------------------------------------------

enum class SolverKind {
    Direct = 0,     // SimplicialLDLT - equivalent de scipy spsolve pour un SPD
    Iterative = 1,  // ConjugateGradient preconditionne par la diagonale
};

struct SolveResult {
    std::vector<double> values;
    bool converged = false;
    int iterations = 0;
    double error = 0.0;
};

// `A` en CSR (indptr/indices/data), symetrique definie positive par
// construction. `initial_guess` a nullptr = depart a zero, comme scipy.
SolveResult solve_linear_system(const std::int32_t* indptr, const std::int32_t* indices,
                                const double* data, std::int64_t node_count, const double* b,
                                SolverKind kind, double tolerance, int max_iterations,
                                const double* initial_guess);

}  // namespace apex

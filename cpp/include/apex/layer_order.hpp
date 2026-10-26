// Ordres de couches (etapes 5 et 5b) - noyaux combinatoires portes.
//
// `stable_topological_order` : tri topologique STABLE de step5b (Kahn, priorite
// = indice d'origine ; cycle rompu par la couche de plus petit degre entrant
// restant, a egalite la mieux classee). Resultat identique au Python (aucune
// arithmetique flottante).
#pragma once

#include <cstdint>
#include <vector>

namespace apex {

struct TopologicalOrder {
    std::vector<std::int64_t> order;
    std::int64_t broken = 0;   // contraintes violees faute d'ordre compatible
};

// Aretes (a[i] -> b[i]) : `a` doit preceder `b`. Les doublons sont ignores,
// comme dans le `set` du Python.
TopologicalOrder stable_topological_order(std::int64_t n, const std::int64_t* a,
                                          const std::int64_t* b, std::int64_t edge_count);

}  // namespace apex

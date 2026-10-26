#include "apex/layer_order.hpp"

#include <algorithm>
#include <functional>
#include <queue>

namespace apex {

TopologicalOrder stable_topological_order(std::int64_t n, const std::int64_t* a,
                                          const std::int64_t* b, std::int64_t edge_count) {
    std::vector<std::pair<std::int64_t, std::int64_t>> edges(static_cast<std::size_t>(edge_count));
    for (std::int64_t i = 0; i < edge_count; ++i) edges[static_cast<std::size_t>(i)] = {a[i], b[i]};
    std::sort(edges.begin(), edges.end());
    edges.erase(std::unique(edges.begin(), edges.end()), edges.end());

    std::vector<std::vector<std::int64_t>> succ(static_cast<std::size_t>(n));
    std::vector<std::int64_t> indeg(static_cast<std::size_t>(n), 0);
    for (const auto& [from, to] : edges) {
        succ[static_cast<std::size_t>(from)].push_back(to);
        ++indeg[static_cast<std::size_t>(to)];
    }

    std::priority_queue<std::int64_t, std::vector<std::int64_t>, std::greater<>> ready;
    for (std::int64_t i = 0; i < n; ++i) {
        if (indeg[static_cast<std::size_t>(i)] == 0) ready.push(i);
    }
    std::vector<char> placed(static_cast<std::size_t>(n), 0);
    TopologicalOrder result;
    result.order.reserve(static_cast<std::size_t>(n));
    std::int64_t remaining = n;
    while (remaining > 0) {
        if (ready.empty()) {
            // cycle : plus petit (degre entrant restant, indice) parmi les non places
            std::int64_t forced = -1;
            for (std::int64_t r = 0; r < n; ++r) {
                if (placed[static_cast<std::size_t>(r)]) continue;
                if (forced < 0 || indeg[static_cast<std::size_t>(r)] < indeg[static_cast<std::size_t>(forced)]) forced = r;
            }
            ready.push(forced);
        }
        const std::int64_t i = ready.top();
        ready.pop();
        if (placed[static_cast<std::size_t>(i)]) continue;
        placed[static_cast<std::size_t>(i)] = 1;
        --remaining;
        result.order.push_back(i);
        for (const std::int64_t j : succ[static_cast<std::size_t>(i)]) {
            if (--indeg[static_cast<std::size_t>(j)] == 0 && !placed[static_cast<std::size_t>(j)]) ready.push(j);
        }
    }
    std::vector<std::int64_t> position(static_cast<std::size_t>(n));
    for (std::int64_t k = 0; k < n; ++k) position[static_cast<std::size_t>(result.order[static_cast<std::size_t>(k)])] = k;
    for (const auto& [from, to] : edges) {
        if (position[static_cast<std::size_t>(from)] > position[static_cast<std::size_t>(to)]) ++result.broken;
    }
    return result;
}

}  // namespace apex

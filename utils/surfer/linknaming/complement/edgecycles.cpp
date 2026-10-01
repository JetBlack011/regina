//
//  edgecycles.cpp
//

#include "linknaming/complement/edgecycles.h"

#include <algorithm>
#include <functional>
#include <unordered_map>
#include <unordered_set>
#include <utility>

namespace edgecycles {

std::vector<EdgeEnds> endsOf(const std::vector<const regina::Edge<3> *> &edges) {
    std::vector<EdgeEnds> out;
    out.reserve(edges.size());
    for (const regina::Edge<3> *e : edges)
        out.push_back({e->index(), e->vertex(0)->index(), e->vertex(1)->index()});
    return out;
}

std::optional<std::vector<std::vector<size_t>>>
chainDirected(const std::vector<EdgeEnds> &edges, OpenChain open) {
    std::unordered_map<size_t, size_t> outFrom; // tail -> position (last wins)
    for (size_t i = 0; i < edges.size(); ++i)
        outFrom[edges[i].v0] = i;

    std::vector<std::vector<size_t>> curves;
    std::unordered_set<size_t> visited; // ids
    for (const EdgeEnds &e : edges) {
        if (visited.contains(e.id))
            continue;
        std::vector<size_t> curve;
        const size_t start = e.v0;
        size_t curr = start;
        for (size_t step = 0; step <= edges.size(); ++step) {
            auto it = outFrom.find(curr);
            if (it == outFrom.end()) {
                if (open == OpenChain::refuse)
                    return std::nullopt;
                break;
            }
            const size_t next = it->second;
            visited.insert(edges[next].id);
            curve.push_back(next);
            curr = edges[next].v1;
            if (curr == start)
                break;
        }
        curves.push_back(std::move(curve));
    }
    return curves;
}

std::optional<size_t> countDirectedCycles(const std::vector<EdgeEnds> &edges) {
    std::unordered_map<size_t, size_t> outOf, inCount;
    for (const EdgeEnds &e : edges) {
        if (!outOf.emplace(e.v0, e.v1).second)
            return std::nullopt; // two edges leave one vertex
        ++inCount[e.v1];
    }
    for (const auto &[v, n] : inCount)
        if (n != 1 || !outOf.contains(v))
            return std::nullopt;
    if (outOf.size() != inCount.size())
        return std::nullopt;
    std::unordered_map<size_t, bool> seen;
    size_t cycles = 0;
    for (const auto &[start, next] : outOf) {
        if (seen[start])
            continue;
        ++cycles;
        for (size_t v = start; !seen[v]; v = outOf.at(v))
            seen[v] = true;
    }
    return cycles;
}

namespace {

// (vertex, position) for both ends of every edge, sorted: each vertex's
// incident positions, in input order, are one contiguous run.
std::vector<std::pair<size_t, size_t>> incidences(const std::vector<EdgeEnds> &edges) {
    std::vector<std::pair<size_t, size_t>> at;
    at.reserve(2 * edges.size());
    for (size_t i = 0; i < edges.size(); ++i) {
        at.emplace_back(edges[i].v0, i);
        at.emplace_back(edges[i].v1, i);
    }
    std::sort(at.begin(), at.end());
    return at;
}

// The first position not yet used among those at vertex `v`, or npos.
size_t firstUnusedAt(const std::vector<std::pair<size_t, size_t>> &at, size_t v,
                     const std::vector<char> &used) {
    auto it = std::lower_bound(at.begin(), at.end(), std::make_pair(v, size_t(0)));
    for (; it != at.end() && it->first == v; ++it)
        if (!used[it->second])
            return it->second;
    return static_cast<size_t>(-1);
}

} // namespace

std::vector<std::vector<Step>> walkCurves(const std::vector<EdgeEnds> &edges) {
    const auto at = incidences(edges);
    std::vector<char> used(edges.size(), 0);
    std::vector<std::vector<Step>> curves;
    for (size_t s = 0; s < edges.size(); ++s) {
        if (used[s])
            continue;
        used[s] = 1;
        std::vector<Step> curve{{s, false}};
        size_t curr = edges[s].v1;
        for (;;) {
            const size_t j = firstUnusedAt(at, curr, used);
            if (j == static_cast<size_t>(-1))
                break;
            used[j] = 1;
            const bool reversed = edges[j].v0 != curr;
            curve.push_back({j, reversed});
            curr = reversed ? edges[j].v0 : edges[j].v1;
        }
        curves.push_back(std::move(curve));
    }
    return curves;
}

std::optional<std::vector<Step>> walkClosedCurve(const std::vector<EdgeEnds> &edges) {
    if (edges.empty())
        return std::nullopt;
    const auto at = incidences(edges);
    for (size_t k = 0; k < at.size();) {
        size_t m = k;
        while (m < at.size() && at[m].first == at[k].first)
            ++m;
        if (m - k != 2)
            return std::nullopt; // not a simple closed curve
        k = m;
    }
    std::vector<char> used(edges.size(), 0);
    std::vector<Step> steps;
    steps.reserve(edges.size());
    size_t slot = 0;
    const size_t start = edges[0].v0;
    size_t from = start;
    for (size_t n = 0; n < edges.size(); ++n) {
        used[slot] = 1;
        const bool reversed = edges[slot].v0 != from;
        steps.push_back({slot, reversed});
        from = reversed ? edges[slot].v0 : edges[slot].v1;
        if (n + 1 == edges.size())
            break;
        slot = firstUnusedAt(at, from, used);
        if (slot == static_cast<size_t>(-1))
            return std::nullopt; // more than one curve
    }
    if (from != start)
        return std::nullopt;
    return steps;
}

std::optional<size_t> countClosedCurves(const std::vector<EdgeEnds> &edges) {
    std::unordered_set<size_t> seen;
    std::unordered_map<size_t, int> degree;
    std::unordered_map<size_t, size_t> parent;
    std::function<size_t(size_t)> find = [&](size_t x) {
        while (parent[x] != x)
            x = parent[x] = parent[parent[x]];
        return x;
    };
    for (const EdgeEnds &e : edges) {
        if (!seen.insert(e.id).second)
            return std::nullopt;
        ++degree[e.v0];
        ++degree[e.v1];
        parent.try_emplace(e.v0, e.v0);
        parent.try_emplace(e.v1, e.v1);
        parent[find(e.v0)] = find(e.v1);
    }
    size_t roots = 0;
    for (const auto &[v, d] : degree) {
        if (d != 2)
            return std::nullopt;
        if (find(v) == v)
            ++roots;
    }
    return roots;
}

} // namespace edgecycles

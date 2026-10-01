//
//  gaussdiagram.cpp
//

#include "linknaming/diagrams/gaussdiagram.h"

#include <algorithm>
#include <array>
#include <map>
#include <numeric>

namespace exactnaming {

regina::Link GaussDiagram::link() const {
    return regina::Link::fromData(signs.begin(), signs.end(), comps.begin(), comps.end());
}

GaussDiagram GaussDiagram::of(const regina::Link &l, std::vector<size_t> origin) {
    GaussDiagram g;
    g.origin = std::move(origin);
    g.signs.reserve(l.size());
    for (size_t k = 0; k < l.size(); ++k) g.signs.push_back(l.crossing(k)->sign());
    g.comps.resize(l.countComponents());
    for (size_t c = 0; c < l.countComponents(); ++c) {
        const regina::StrandRef start = l.component(c);
        if (!start) continue;
        regina::StrandRef s = start;
        do {
            const long k = static_cast<long>(s.crossing()->index()) + 1;
            g.comps[c].push_back(s.strand() == 1 ? k : -k);
            s = s.next();
        } while (s != start);
    }
    return g;
}

namespace {

// The sub-diagram on the crossings in `keep` (old indices), made of the given
// component sequences (already restricted to those crossings), renumbered.
GaussDiagram restrict(const GaussDiagram &d, const std::vector<bool> &keep,
                      const std::vector<std::vector<long>> &comps,
                      const std::vector<size_t> &origin) {
    GaussDiagram out;
    std::vector<long> renum(d.crossings(), -1);
    for (size_t k = 0; k < d.crossings(); ++k)
        if (keep[k]) {
            renum[k] = static_cast<long>(out.signs.size());
            out.signs.push_back(d.signs[k]);
        }
    for (const auto &seq : comps) {
        std::vector<long> &o = out.comps.emplace_back();
        for (long v : seq) {
            const long k = renum[static_cast<size_t>(std::abs(v) - 1)] + 1;
            o.push_back(v > 0 ? k : -k);
        }
    }
    out.origin = origin;
    return out;
}

} // namespace

std::vector<GaussDiagram> splitPieces(const GaussDiagram &d) {
    const size_t n = d.components();
    std::vector<size_t> parent(n);
    std::iota(parent.begin(), parent.end(), 0);
    auto find = [&](size_t x) {
        while (parent[x] != x) x = parent[x] = parent[parent[x]];
        return x;
    };
    std::vector<long> firstComp(d.crossings(), -1);
    for (size_t c = 0; c < n; ++c)
        for (long v : d.comps[c]) {
            long &f = firstComp[static_cast<size_t>(std::abs(v) - 1)];
            if (f < 0) f = static_cast<long>(c);
            else parent[find(c)] = find(static_cast<size_t>(f));
        }
    std::map<size_t, std::vector<size_t>> groups; // root -> components, in order
    for (size_t c = 0; c < n; ++c) groups[find(c)].push_back(c);
    std::vector<GaussDiagram> out;
    for (const auto &[root, members] : groups) {
        std::vector<bool> keep(d.crossings(), false);
        std::vector<std::vector<long>> comps;
        std::vector<size_t> origin;
        for (size_t c : members) {
            for (long v : d.comps[c]) keep[static_cast<size_t>(std::abs(v) - 1)] = true;
            comps.push_back(d.comps[c]);
            origin.push_back(d.origin[c]);
        }
        out.push_back(restrict(d, keep, comps, origin));
    }
    return out;
}

std::optional<std::pair<GaussDiagram, GaussDiagram>> visibleSum(const GaussDiagram &d) {
    // Half-edge (crossing, position): positions run counterclockwise from the
    // incoming under-strand -- [under-in, over-in, under-out, over-out] at a
    // left-handed (-1) crossing, [under-in, over-out, under-out, over-in] at
    // a right-handed one (knotbuilder::DiagramDrawer's and the KnotTheory
    // convention).
    auto posIn = [&](size_t k, bool over) { return over ? (d.signs[k] < 0 ? 1 : 3) : 0; };
    auto posOut = [&](size_t k, bool over) { return over ? (d.signs[k] < 0 ? 3 : 1) : 2; };

    // Arcs: one per visit, from that visit to the next along its component.
    struct Arc { size_t comp, from; std::array<size_t, 2> half; };
    std::vector<Arc> arcs;
    std::vector<long> arcAt(4 * d.crossings(), -1); // half-edge -> arc
    std::vector<size_t> otherEnd(4 * d.crossings(), 0);
    for (size_t c = 0; c < d.components(); ++c) {
        const auto &seq = d.comps[c];
        for (size_t p = 0; p < seq.size(); ++p) {
            const long v = seq[p], w = seq[(p + 1) % seq.size()];
            const size_t x = static_cast<size_t>(std::abs(v) - 1), y = static_cast<size_t>(std::abs(w) - 1);
            const size_t h0 = 4 * x + posOut(x, v > 0), h1 = 4 * y + posIn(y, w > 0);
            arcAt[h0] = arcAt[h1] = static_cast<long>(arcs.size());
            otherEnd[h0] = h1;
            otherEnd[h1] = h0;
            arcs.push_back({c, p, {h0, h1}});
        }
    }
    for (long a : arcAt)
        if (a < 0) return std::nullopt; // not a closed 4-valent diagram; show nothing

    // Faces: from a half-edge, along its arc, then counterclockwise.
    std::vector<long> face(4 * d.crossings(), -1);
    long faces = 0;
    for (size_t h = 0; h < face.size(); ++h) {
        if (face[h] >= 0) continue;
        for (size_t g = h; face[g] < 0;) {
            face[g] = faces;
            const size_t e = otherEnd[g];
            g = 4 * (e / 4) + (e % 4 + 1) % 4;
        }
        ++faces;
    }

    // Arcs by the (unordered) pair of faces on their two sides.
    std::map<std::pair<long, long>, std::vector<size_t>> bySides;
    for (size_t a = 0; a < arcs.size(); ++a) {
        long f0 = face[arcs[a].half[0]], f1 = face[arcs[a].half[1]];
        if (f0 == f1) continue;
        bySides[{std::min(f0, f1), std::max(f0, f1)}].push_back(a);
    }
    for (const auto &[sides, as] : bySides)
        for (size_t i = 0; i < as.size(); ++i)
            for (size_t j = i + 1; j < as.size(); ++j) {
                const Arc &e = arcs[as[i]], &f = arcs[as[j]];
                if (e.comp != f.comp) continue; // cannot happen for a planar diagram
                // Crossings joined by every arc but e and f.
                std::vector<std::vector<size_t>> adj(d.crossings());
                for (size_t a = 0; a < arcs.size(); ++a) {
                    if (a == as[i] || a == as[j]) continue;
                    size_t x = arcs[a].half[0] / 4, y = arcs[a].half[1] / 4;
                    adj[x].push_back(y);
                    adj[y].push_back(x);
                }
                std::vector<bool> sideA(d.crossings(), false);
                std::vector<size_t> stack{e.half[1] / 4};
                sideA[stack.back()] = true;
                while (!stack.empty()) {
                    size_t x = stack.back();
                    stack.pop_back();
                    for (size_t y : adj[x])
                        if (!sideA[y]) { sideA[y] = true; stack.push_back(y); }
                }
                if (std::all_of(sideA.begin(), sideA.end(), [](bool b) { return b; }))
                    continue; // not a cut

                // The summed component: its visits after e up to f on one
                // side (e.half[1] is side A's), the rest on the other.
                const auto &seq = d.comps[e.comp];
                const size_t m = seq.size();
                std::vector<long> runA, runB;
                for (size_t k = 1; k <= m; ++k) {
                    const size_t p = (e.from + k) % m;
                    (sideA[static_cast<size_t>(std::abs(seq[p]) - 1)] ? runA : runB).push_back(seq[p]);
                }
                std::vector<bool> sideB(d.crossings());
                for (size_t k = 0; k < d.crossings(); ++k) sideB[k] = !sideA[k];
                std::vector<std::vector<long>> compsA, compsB;
                std::vector<size_t> originA, originB;
                for (size_t c = 0; c < d.components(); ++c) {
                    if (c == e.comp) {
                        compsA.push_back(runA); originA.push_back(d.origin[c]);
                        compsB.push_back(runB); originB.push_back(d.origin[c]);
                    } else if (!d.comps[c].empty() &&
                               sideA[static_cast<size_t>(std::abs(d.comps[c][0]) - 1)]) {
                        compsA.push_back(d.comps[c]); originA.push_back(d.origin[c]);
                    } else {
                        compsB.push_back(d.comps[c]); originB.push_back(d.origin[c]);
                    }
                }
                return std::make_pair(restrict(d, sideA, compsA, originA),
                                      restrict(d, sideB, compsB, originB));
            }
    return std::nullopt;
}

} // namespace exactnaming

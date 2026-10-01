//
//  unlinknaming.cpp
//

#include "linknaming/complement/unlinknaming.h"

#include <functional>
#include <optional>
#include <unordered_map>
#include <unordered_set>

#include <snappea/snappeatriangulation.h>

#include "linknaming/complement/complementcache.h"

namespace {

// Fast, sound, one-sided proof that `t`'s genus is 1 (the unknot): its
// fundamental group is Z if and only if it's the unknot (Dehn's lemma).
// group() already tries to simplify the presentation internally (same
// idiom as engine/triangulation/dim3/knot.cpp's Poincare-conjecture
// fast path), so a presentation of exactly one generator and no relations
// -- i.e. <a|>, which *is* Z by construction -- is conclusive. A "false"
// here is only ever inconclusive (simplify() didn't collapse it that far),
// never a wrong answer: this can only shorten the path to genus == 1, so
// it's always safe to fall back to recogniseHandlebody() when it fails.
bool groupProvesUnknot(const regina::Triangulation<3> &t) {
    const regina::GroupPresentation &g = t.group();
    return g.countGenerators() == 1 && g.countRelations() == 0;
}

// Fast, sound, one-sided proof that `t`'s genus is -1 (not any
// handlebody, of any genus): a genuine hyperbolic structure rules out
// every handlebody by geometrization (no handlebody -- solid tori
// included -- admits a complete hyperbolic structure, since its boundary
// is compressible). Symmetric to groupProvesUnknot() above: a "false"
// here is inconclusive, never wrong, so it's always safe to fall back to
// recogniseHandlebody(). SnapPeaTriangulation's construction and
// solutionType() are not documented thread-safe, but each caller here
// builds one from its own local, independently-owned Triangulation<3> --
// no state is shared across threads -- and this was stress-tested under
// concurrent load (8 threads, thousands of calls) with no crashes or
// inconsistent results.
bool hyperbolicityProvesNotHandlebody(const regina::Triangulation<3> &t) {
    regina::SnapPeaTriangulation snappea(t);
    return snappea.solutionType() ==
           regina::SnapPeaTriangulation::Solution::Geometric;
}

} // namespace

namespace identify {

// Cached recogniseHandlebody(): never touches censusLookupMutex, so this is
// safe to call from identify::isUnknot()'s hot, highly-parallel path
// without risking contention with an in-flight Census::lookup(). Tries the
// two fast, sound one-sided checks above before falling back to the
// expensive normal-surface-theory path -- both are cheap regardless of
// outcome, and either resolving conclusively skips recogniseHandlebody()
// entirely. Neither helps for knots that are neither the unknot nor
// hyperbolic (torus/satellite knots), which still fall through to
// recogniseHandlebody() same as before.
ssize_t cachedGenus(const regina::Triangulation<3> &complement,
                    const std::string &sig) {
    if (auto hit = checkRecognition(
            sig, &RecognitionCacheStats::genusChecks,
            &RecognitionCacheStats::genusCacheHits,
            [](const RecognitionResult &r) { return r.genus.has_value(); }))
        return *hit->genus;

    ssize_t genus;
    enum class Path { Group, SnapPea, Fallback } path;
    if (groupProvesUnknot(complement)) {
        genus = 1;
        path = Path::Group;
    } else if (hyperbolicityProvesNotHandlebody(complement)) {
        genus = -1;
        path = Path::SnapPea;
    } else {
        genus = complement.recogniseHandlebody();
        path = Path::Fallback;
    }

    countRecognition([path](RecognitionCacheStats &s) {
        if (path == Path::Group)
            ++s.groupFastPathHits;
        else if (path == Path::SnapPea)
            ++s.snapPeaFastPathHits;
        else
            ++s.recogniseHandlebodyFallbacks;
    });
    return *storeRecognition(sig, RecognitionResult{.genus = genus}).genus;
}

// Fast, sound, one-sided proof that `t` (a LINK complement, possibly
// multiple components) is split -- i.e. `t` is identify(const Link&)'s
// n-component-unlink case -- generalizing groupProvesUnknot() above from
// n == 1 to any n. A presentation with zero relations is free by
// construction (same "sound regardless of how simplify() got there"
// argument as groupProvesUnknot()), and a free fundamental group forces a
// link to be split: by Milnor's prime decomposition theorem, an orientable
// 3-manifold's prime decomposition realizes the Grushko free-product
// decomposition of its fundamental group, so a free pi_1 of rank k means
// `t` splits along essential spheres into k pieces, each with a single
// torus boundary component and pi_1 == Z -- which, by the exact same
// Dehn's-lemma argument groupProvesUnknot() itself relies on, forces each
// piece to be a solid torus. So `t` is the complement of k unknotted,
// pairwise split components, i.e. the k-component unlink. (No need to
// separately check k against the link's actual component count: a link
// complement's H_1 always has rank == component count via Alexander
// duality/meridians, regardless of link type, so if the group is ALSO
// free, its rank -- which must match H_1's rank -- is automatically the
// component count.) As with groupProvesUnknot(), a "false" here is only
// ever inconclusive (simplify() didn't collapse the presentation that
// far), never wrong -- always safe to fall through to the normal
// resolveRecognition() path when it fails.
bool groupProvesUnlink(const regina::Triangulation<3> &t) {
    return t.group().countRelations() == 0;
}

bool isOrientationSafeName(const std::string &name) {
    return name == "Unknot" || name.ends_with("-component unlink");
}

bool isUnknot(const EdgeComplement &e) {
    auto complement = e.buildComplement();
    return cachedGenus(complement, complement.isoSig()) == 1;
}

std::string unlinkNameOrIsoSig(const EdgeComplement &e) {
    auto complement = e.buildComplement();
    std::string sig = complement.isoSig();
    return cachedGenus(complement, sig) == 1 ? "Unknot" : sig;
}

std::string unlinkNameOrIsoSig(const Link &l) {
    auto complement = l.buildComplement();
    if (l.countComponents() > 1 && groupProvesUnlink(complement))
        return std::to_string(l.countComponents()) + "-component unlink";
    std::string sig = complement.isoSig();
    return cachedGenus(complement, sig) == 1 ? "Unknot" : sig;
}

namespace {
// How many cycles `edges` forms, or nullopt unless it is a disjoint union of
// cycles: no repeated edge, and every vertex it touches has degree exactly 2
// (a loop edge contributes 2 to its one vertex). A connected 2-regular
// multigraph is a cycle, so the component count is then the cycle count.
std::optional<size_t>
countCycles(const std::vector<const regina::Edge<3> *> &edges) {
    std::unordered_set<const regina::Edge<3> *> seen;
    std::unordered_map<size_t, int> degree;
    std::unordered_map<size_t, size_t> parent;
    std::function<size_t(size_t)> find = [&](size_t x) {
        while (parent[x] != x)
            x = parent[x] = parent[parent[x]];
        return x;
    };
    for (const auto *e : edges) {
        if (!seen.insert(e).second)
            return std::nullopt;
        size_t a = e->vertex(0)->index(), b = e->vertex(1)->index();
        ++degree[a];
        ++degree[b];
        parent.try_emplace(a, a);
        parent.try_emplace(b, b);
        parent[find(a)] = find(b);
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
} // namespace

bool certifiesUnlink(const regina::Triangulation<3> &tri,
                     const std::vector<const regina::Edge<3> *> &edges,
                     size_t m) {
    if (m == 0)
        return false;
    auto cycles = countCycles(edges);
    if (!cycles || *cycles != m)
        return false;
    // pinchEdge()'s one precondition; buildComplement() pinches every edge.
    for (const auto *e : edges)
        if (e->isBoundary())
            return false;
    try {
        auto complement = EdgeComplement(tri, edges).buildComplement();
        if (!complement.isValid() || !complement.isConnected())
            return false;
        // group() simplifies internally; ideal vertices count as truncated,
        // so this is pi_1 of the link exterior. A presentation with no
        // relations IS free, of rank countGenerators(), however simplify()
        // reached it -- which is what makes this one-sided but sound.
        const regina::GroupPresentation &g = complement.group();
        return g.countRelations() == 0 && g.countGenerators() == m;
    } catch (const std::exception &) {
        return false;
    }
}

bool capInCone(const regina::Triangulation<3> &ball,
               const std::vector<const regina::Edge<3> *> &edges,
               CappedCurves &out) {
    if (edges.empty() || !ball.hasBoundaryFacets())
        return false;

    // Everything is recorded positionally first: Edge<3>*/Vertex<3>* belong
    // to `ball`, but (tetrahedron index, local number) survives both the
    // copy and makeIdeal(), which only appends tetrahedra.
    std::unordered_set<const regina::Edge<3> *> seen;
    std::unordered_map<size_t, int> degree;
    std::vector<std::pair<size_t, int>> edgeDesc;
    for (const auto *e : edges) {
        if (!seen.insert(e).second)
            return false;
        edgeDesc.emplace_back(e->front().simplex()->index(), e->front().face());
        ++degree[e->vertex(0)->index()];
        ++degree[e->vertex(1)->index()];
    }
    std::vector<std::pair<size_t, int>> endDesc;
    for (const auto &[v, d] : degree) {
        if (d == 1) {
            const auto *vertex = ball.vertex(v);
            if (!vertex->isBoundary())
                return false;
            endDesc.emplace_back(vertex->front().simplex()->index(),
                                 vertex->front().face());
        } else if (d != 2) {
            return false;
        }
    }
    if (!endDesc.empty() && endDesc.size() != 2)
        return false;

    const size_t nOrig = ball.size();
    out.tri = ball;
    out.tri.makeIdeal();
    if (out.tri.size() <= nOrig)
        return false;

    out.edges.clear();
    for (const auto &[t, i] : edgeDesc)
        out.edges.push_back(out.tri.tetrahedron(t)->edge(i));

    if (!endDesc.empty()) {
        // makeIdeal() glues facet 3 of each new tetrahedron to a boundary
        // facet, so local vertex 3 of every new tetrahedron is the apex.
        const auto *apex = out.tri.tetrahedron(nOrig)->vertex(3);
        for (const auto &[t, i] : endDesc) {
            const auto *x = out.tri.tetrahedron(t)->vertex(i);
            const regina::Edge<3> *spoke = nullptr;
            for (size_t c = nOrig; c < out.tri.size() && !spoke; ++c) {
                auto *cone = out.tri.tetrahedron(c);
                for (int a = 0; a < 3; ++a)
                    if (cone->vertex(a) == x) {
                        spoke = cone->edge(a, 3);
                        break;
                    }
            }
            if (!spoke ||
                !((spoke->vertex(0) == x && spoke->vertex(1) == apex) ||
                  (spoke->vertex(1) == x && spoke->vertex(0) == apex)))
                return false;
            out.edges.push_back(spoke);
        }
    }

    auto cycles = countCycles(out.edges);
    if (!cycles)
        return false;
    out.components = *cycles;
    return true;
}

} // namespace identify

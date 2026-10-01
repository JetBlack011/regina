// predicate_order_test.cpp
//
// KnottedSurface::addFace() prunes on P_1, P_flat and P_transverse (paper
// Table tab:predicates). Each is anti-monotonic (Lemmas p1-anti-monotonic and
// smoothness-anti-monotonic), so whether a set S of triangles can be added
// is a property of S alone: it must not depend on the order the faces are
// added in, and it must agree with the definitions evaluated on S from
// scratch. EmbeddingSearch relies on both: the enumerator reaches a set
// along one parent chain, and the drain rebuilds it in another order
// (surfacesearch.cpp, processEntry_).
//
// The reference below computes, from Regina's own face mappings and nothing
// of EmbeddedSubmanifold's bookkeeping:
//   - P_1: every edge has valence <= 2, counting (triangle, edge) pairs;
//   - petals: corners at a vertex, identified across each valence-2 edge;
//   - a petal is closed when both spokes of each of its corners have
//     valence 2, and its trace is the link edge of every one of its corners;
//   - at every interior vertex, every closed petal's trace is an unknot and
//     every two closed petals' traces have linking number 0.
// Unknot recognition and linking numbers are the production ones
// (identify::isUnknot, EdgeComplement::linkingNumberWith): this test is about
// which corners make up a petal and when it is checked, not about knot
// recognition.
//
// The one-vertex closed triangulations matter most: every triangle there has
// all three corners at one interior vertex, so addFace() must register every
// corner before checking any petal. Checking corner by corner once pruned
// valid sets there, depended on the order, and threw from the linking check.

#include <algorithm>
#include <functional>
#include <iostream>
#include <map>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include <triangulation/dim3.h>
#include <triangulation/dim4.h>
#include <triangulation/example4.h>

#include "surfer/submanifold/submanifold.h"
#include "linknaming/census/identifycomplement.h"
#include "linknaming/complement/linkcomplement.h"
#include "surfer/submanifold/skeleton.h"

namespace {

int failures = 0;

using Corner = std::pair<int, int>; // (skeleton face index, local vertex)

// The edge of Lk(v) traced by corner (f, i): the edge opposite vertex i of
// triangle f, seen from one pentachoron containing f.
const regina::Edge<3> *linkEdge(const regina::Triangle<4> *tri, int i) {
    const regina::Vertex<4> *v = tri->vertex(i);
    const auto &emb = tri->front();
    const regina::Pentachoron<4> *pent = emb.simplex();
    const regina::Perm<5> perm = emb.vertices();
    const int vLocal = perm[i];
    const int w1 = perm[(i + 1) % 3];
    const int w2 = perm[(i + 2) % 3];

    size_t idx = 0;
    bool found = false;
    for (const auto &ve : v->embeddings()) {
        if (ve.simplex() == pent && ve.face() == vLocal) {
            found = true;
            break;
        }
        ++idx;
    }
    if (!found)
        throw std::logic_error("linkEdge: vertex embedding not found");

    const regina::Tetrahedron<3> *tet = v->buildLink().tetrahedron(idx);
    const regina::Perm<5> inv = pent->tetrahedronMapping(vLocal).inverse();
    return tet->edge(inv[w1], inv[w2]);
}

// Whether S satisfies P_1 and P_flat and P_transverse, from the definitions.
bool referenceAccepts(const regina::Triangulation<4> &tri,
                      const Skeleton<4, 2> &skel, const std::vector<int> &S) {
    const auto &nodes = skel.getNodes();

    // P_1, with the (triangle, edge) pairs over each ambient edge.
    std::map<size_t, std::vector<std::pair<int, int>>> pairs; // edge -> (f, j)
    for (int f : S) {
        const regina::Triangle<4> *t = nodes[f].face;
        for (int j = 0; j < 3; ++j)
            pairs[t->edge(j)->index()].push_back({f, j});
    }
    for (const auto &[e, ps] : pairs)
        if (ps.size() > 2)
            return false;

    // Corners, and the union-find identifying them across valence-2 edges.
    std::vector<Corner> corners;
    std::map<Corner, int> cornerId;
    for (int f : S)
        for (int i = 0; i < 3; ++i) {
            cornerId[{f, i}] = static_cast<int>(corners.size());
            corners.push_back({f, i});
        }
    std::vector<int> parent(corners.size());
    std::iota(parent.begin(), parent.end(), 0);
    std::function<int(int)> find = [&](int x) {
        return parent[x] == x ? x : parent[x] = find(parent[x]);
    };
    for (const auto &[e, ps] : pairs) {
        if (ps.size() != 2)
            continue;
        // Each end of the ambient edge carries one corner of each side.
        for (int end = 0; end < 2; ++end) {
            const auto [fa, ja] = ps[0];
            const auto [fb, jb] = ps[1];
            const int ia = nodes[fa].face->edgeMapping(ja)[end];
            const int ib = nodes[fb].face->edgeMapping(jb)[end];
            parent[find(cornerId[{fa, ia}])] = find(cornerId[{fb, ib}]);
        }
    }

    auto valence = [&](size_t edge) {
        auto it = pairs.find(edge);
        return it == pairs.end() ? 0 : static_cast<int>(it->second.size());
    };

    // Petals per interior vertex.
    std::map<size_t, std::map<int, std::vector<int>>> petals; // v -> root -> corners
    for (size_t c = 0; c < corners.size(); ++c) {
        const auto [f, i] = corners[c];
        const regina::Vertex<4> *v = nodes[f].face->vertex(i);
        if (v->isBoundary())
            continue;
        petals[v->index()][find(static_cast<int>(c))].push_back(
            static_cast<int>(c));
    }

    for (const auto &[vIdx, byRoot] : petals) {
        const regina::Vertex<4> *v = tri.vertex(vIdx);
        std::vector<std::vector<const regina::Edge<3> *>> closedTraces;
        for (const auto &[root, members] : byRoot) {
            bool closed = true;
            std::vector<const regina::Edge<3> *> trace;
            for (int c : members) {
                const auto [f, i] = corners[c];
                const regina::Triangle<4> *t = nodes[f].face;
                for (int j = 0; j < 3; ++j)
                    if (j != i && valence(t->edge(j)->index()) != 2)
                        closed = false;
                trace.push_back(linkEdge(t, i));
            }
            if (!closed)
                continue;
            if (!identify::isUnknot(Knot(v->buildLink(), trace)))
                return false;
            closedTraces.push_back(std::move(trace));
        }
        for (size_t a = 0; a < closedTraces.size(); ++a)
            for (size_t b = a + 1; b < closedTraces.size(); ++b)
                if (Knot(v->buildLink(), closedTraces[a])
                        .linkingNumberWith(
                            Knot(v->buildLink(), closedTraces[b])) != 0)
                    return false;
    }
    return true;
}

// A throw out of addFace() is a failure in its own right: in the search it
// escapes a worker thread and terminates the process.
long throws = 0;

bool knottedAccepts(const Skeleton<4, 2> &skel, const std::vector<int> &order) {
    KnottedSurface s(skel);
    try {
        return s.addFaces(order);
    } catch (const std::exception &e) {
        if (++throws <= 3) {
            std::cout << "  addFaces(";
            for (size_t i = 0; i < order.size(); ++i)
                std::cout << (i ? "," : "") << order[i];
            std::cout << ") THREW: " << e.what() << std::endl;
        }
        return false;
    }
}

std::string show(const std::vector<int> &v) {
    std::string out = "{";
    for (size_t i = 0; i < v.size(); ++i)
        out += (i ? "," : "") + std::to_string(v[i]);
    return out + "}";
}

void checkTriangulation(const std::string &name, const regina::Triangulation<4> &tri,
                        size_t maxSize) {
    Skeleton<4, 2> skel(tri);
    std::vector<int> candidates;
    for (size_t i = 0; i < skel.numFaces(); ++i)
        if (!EmbeddedSubmanifold<4, 2>::hasIrreparableSelfGluing(
                skel.getNodes()[i].gluings))
            candidates.push_back(static_cast<int>(i));

    std::cout << "  " << name << " ..." << std::endl;
    std::mt19937 rng(20260926);
    const long throwsBefore = throws;
    long sets = 0, accepted = 0, orderDisagree = 0, referenceDisagree = 0;
    const size_t n = candidates.size();
    // Every subset of candidates up to maxSize, by bitmask when small enough,
    // otherwise by combinations.
    std::vector<int> S;
    std::function<void(size_t)> rec = [&](size_t start) {
        if (!S.empty()) {
            ++sets;
            bool ref = false;
            try {
                ref = referenceAccepts(tri, skel, S);
            } catch (const std::exception &e) {
                std::cout << "  " << name << ": the REFERENCE threw on "
                          << show(S) << ": " << e.what() << std::endl;
                ++failures;
            }
            std::vector<int> fwd = S, rev(S.rbegin(), S.rend()), shuf = S;
            std::shuffle(shuf.begin(), shuf.end(), rng);
            const bool a = knottedAccepts(skel, fwd);
            const bool b = knottedAccepts(skel, rev);
            const bool c = knottedAccepts(skel, shuf);
            accepted += ref;
            if (a != b || a != c) {
                if (++orderDisagree <= 3)
                    std::cout << "  " << name << ": order matters for "
                              << show(S) << " (forward " << a << ", reverse "
                              << b << ", shuffled " << show(shuf) << " " << c
                              << ")\n";
            }
            if (a != ref || b != ref || c != ref) {
                if (++referenceDisagree <= 3)
                    std::cout << "  " << name << ": " << show(S)
                              << " reference " << ref << ", KnottedSurface "
                              << a << "/" << b << "/" << c << "\n";
            }
        }
        if (S.size() == maxSize)
            return;
        for (size_t k = start; k < n; ++k) {
            S.push_back(candidates[k]);
            rec(k + 1);
            S.pop_back();
        }
    };
    rec(0);

    const long threw = throws - throwsBefore;
    const bool ok = orderDisagree == 0 && referenceDisagree == 0 && threw == 0;
    std::cout << (ok ? "  PASS " : "  FAIL ") << name << ": " << sets
              << " sets of at most " << maxSize << " of " << n
              << " triangles, " << accepted << " accepted by the reference, "
              << orderDisagree << " order-dependent, " << referenceDisagree
              << " disagreeing with the reference, " << threw
              << " addFaces() calls threw" << std::endl;
    if (!ok)
        ++failures;
}

} // namespace

int main() {
    std::cout << "KnottedSurface acceptance is a property of the set:\n";
    // Controls: simplicial, so no triangle has a repeated vertex.
    checkTriangulation("fourSphere", regina::Example<4>::fourSphere(), 10);
    checkTriangulation("simplicialFourSphere",
                       regina::Example<4>::simplicialFourSphere(), 5);
    // One-vertex (or few-vertex) closed triangulations: triangles with
    // repeated interior vertices.
    checkTriangulation("s3xs1", regina::Example<4>::s3xs1(), 6);
    checkTriangulation("s3xs1Twisted", regina::Example<4>::s3xs1Twisted(), 6);
    checkTriangulation("cp2", regina::Example<4>::cp2(), 10);
    checkTriangulation("rp4", regina::Example<4>::rp4(), 12);

    if (failures) {
        std::cout << failures << " triangulation(s) FAILED\n";
        return 1;
    }
    std::cout << "all passed\n";
    return 0;
}

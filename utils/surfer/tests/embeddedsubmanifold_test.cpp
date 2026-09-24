// embeddedsubmanifold_test.cpp
//
// Tests whether EmbeddedSubmanifold<4,2>::isEmbedded() (inherited verbatim
// by KnottedSurface) accurately tracks whether the realization map from
// subtri_'s cells into the ambient Triangulation<4>'s cells is injective on
// every open cell. Since Phase 2 (addFace()'s codimension->=2 gating check)
// was disabled, addFace() legitimately accepts faces that leave the growing
// subcomplex NOT genuinely embedded -- that's expected, not a bug; what must
// hold is that isEmbedded() correctly reports false whenever that's true,
// and true only when the subcomplex is genuinely embedded.
//
// The check used here (isGenuinelyEmbedded(), below) is deliberately built
// from scratch against only the PUBLIC API of EmbeddedSubmanifold/
// KnottedSurface (addFace(), removeFace(), triangulation(), boundaryLinks())
// -- it never touches faceCount_ (private) or even faces_ (protected), so it
// shares none of addFace()'s/isEmbedded()'s own bookkeeping and can catch
// bugs in that bookkeeping rather than just re-deriving them. The audit
// (EmbeddednessAuditor, below) flags a violation exactly when isEmbedded()
// disagrees with this ground truth at any point along the DFS.
//
// Two independent verification strategies are used:
//   1. isGenuinelyEmbedded(): a brute-force pointwise-injectivity check,
//      driven through every connected subset of a triangulation's 2-skeleton
//      via the same ConnectedInducedSubgraphEnumerator machinery
//      EmbeddingSearch uses in production (see EmbeddednessAuditor below).
//   2. A boundary-link homology check: for triangulations capped by
//      CobordismBuilder<3>::cone() (so their boundary is S^3), any proper
//      connected embedded surface's boundary is a knot/link in that S^3, and
//      an n-component link complement always has H_1 = Z^n -- an algebraic
//      invariant with nothing to do with addFace()'s internals. Only
//      evaluated at states where isEmbedded() holds, since it presumes a
//      genuinely embedded surface.

#include <algorithm>
#include <chrono>
#include <iostream>
#include <optional>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <unistd.h>
#include <vector>

#include <maths/perm.h>
#include <triangulation/dim2.h>
#include <triangulation/dim3.h>
#include <triangulation/dim4.h>
#include <triangulation/example2.h>
#include <triangulation/example4.h>

#include "cobordismbuilder.h"
#include "embeddedsubmanifold.h"
#include "enumerate_cis.h"
#include "knotbuilder.h"
#include "linkcomplement.h"
#include "skeleton.h"

static int passed = 0, failed_count = 0;

namespace {
bool colorEnabled() {
    static bool enabled = isatty(fileno(stdout));
    return enabled;
}
std::ostream &green(std::ostream &os) {
    return colorEnabled() ? os << "\033[32m" : os;
}
std::ostream &red(std::ostream &os) {
    return colorEnabled() ? os << "\033[31m" : os;
}
std::ostream &bold(std::ostream &os) {
    return colorEnabled() ? os << "\033[1m" : os;
}
std::ostream &resetColor(std::ostream &os) {
    return colorEnabled() ? os << "\033[0m" : os;
}
} // namespace

#define EXPECT_EQ(actual, expected, desc)                                      \
    do {                                                                       \
        auto _a = (actual);                                                    \
        auto _e = (expected);                                                  \
        if (_a == _e) {                                                        \
            std::cout << green << "  PASS: " << resetColor << (desc) << "\n";  \
            ++passed;                                                          \
        } else {                                                               \
            std::cout << red << "  FAIL: " << (desc) << "\n"                   \
                      << "        expected " << _e << ", got " << _a           \
                      << resetColor << "\n";                                   \
            ++failed_count;                                                    \
        }                                                                      \
    } while (0)

// ─────────────────────────────────────────────────────────────────────────────
// Section A: the independent brute-force checker.
// ─────────────────────────────────────────────────────────────────────────────

// Maps each abstract simplex of a KnottedSurface's subtri_ to the ambient
// Face<4,2> (triangle) it represents. newSimplex()+join() (addFace()'s Phase
// 3) always attach the fresh abstract simplex using the ambient face's own
// local facet numbering verbatim, so abstract-simplex local slot j and
// ambient-face local slot j genuinely correspond -- this is a structural
// fact about how addFace() wires up subtri_, not part of the Condition-1/2
// logic under test, so relying on it doesn't compromise independence.
using Correspondence =
    std::unordered_map<const regina::Simplex<2> *, const regina::Face<4, 2> *>;

// Returns true iff the realization map subtri_ -> ambient triangulation
// implied by `correspondence` is injective on vertices, edges, and
// triangles. On failure, *why (if non-null) explains the collision found.
bool isGenuinelyEmbedded(const Correspondence &correspondence,
                         std::string *why = nullptr) {
    // Vertices (k = 0).
    {
        std::unordered_map<size_t, const regina::Vertex<2> *> seen;
        for (const auto &[simplex, face] : correspondence) {
            for (int j = 0; j < 3; ++j) {
                const regina::Vertex<2> *av = simplex->vertex(j);
                const regina::Vertex<4> *gv = face->vertex(j);
                auto [it, inserted] = seen.try_emplace(gv->index(), av);
                if (!inserted && it->second != av) {
                    if (why)
                        *why = "vertex collision at ambient vertex " +
                               std::to_string(gv->index());
                    return false;
                }
            }
        }
    }
    // Edges (k = 1).
    {
        std::unordered_map<size_t, const regina::Edge<2> *> seen;
        for (const auto &[simplex, face] : correspondence) {
            for (int j = 0; j < 3; ++j) {
                const regina::Edge<2> *ae = simplex->edge(j);
                const regina::Edge<4> *ge = face->edge(j);
                auto [it, inserted] = seen.try_emplace(ge->index(), ae);
                if (!inserted && it->second != ae) {
                    if (why)
                        *why = "edge collision at ambient edge " +
                               std::to_string(ge->index());
                    return false;
                }
            }
        }
    }
    // Triangles (k = 2): trivially injective by construction of
    // `correspondence` (an unordered_map keyed by distinct abstract
    // Simplex<2>*) -- verify no ambient triangle was claimed twice, as a
    // self-consistency check on the auditor rather than on addFace().
    {
        std::unordered_set<const regina::Face<4, 2> *> vals;
        for (const auto &[simplex, face] : correspondence) {
            if (!vals.insert(face).second) {
                if (why)
                    *why = "same ambient triangle claimed twice";
                return false;
            }
        }
    }
    return true;
}

// Independent ground truth for isProper(), which is now tracked
// incrementally (see the class comments on facetIsAmbientBoundary_/
// badProperCount_ in embeddedsubmanifold.h) rather than recomputed from
// scratch on every call. This is literally the O(size-of-ambient-
// triangulation) sweep the class used to perform, rewritten against
// `correspondence` instead of the class's private faces_/skeleton_ -- an
// oracle independent of the incremental bookkeeping under test, not a
// re-derivation of it. (boundaryComponentsMapInjectively() is unchanged --
// still the original O(n) sweep -- so it needs no such oracle.)
bool bruteForceIsProper(const Correspondence &correspondence) {
    for (const auto &[simplex, face] : correspondence) {
        for (int i = 0; i < 3; ++i) {
            if (simplex->adjacentSimplex(i) != nullptr)
                continue; // internal facet of subtri_
            if (!face->template face<1>(i)->isBoundary())
                return false;
        }
    }
    return true;
}

// ─────────────────────────────────────────────────────────────────────────────
// Section B: the auditor, driving addFace()/removeFace() through every
// connected subset of the 2-face gluing graph via
// ConnectedInducedSubgraphEnumerator -- the same mechanism
// EmbeddingSearch<4,2>::EmbeddednessPredicate uses in production.
// ─────────────────────────────────────────────────────────────────────────────

struct Violation {
    std::vector<int> faceSequence; // skeleton face indices, in the order added
    std::string reason;
};

class EmbeddednessAuditor : public ConditionalPredicate {
    KnottedSurface &embedding_;
    const Skeleton<4, 2> &skeleton_;
    const std::vector<int> &graphToSkel_;
    int maxDepth_;      // 0 = unbounded
    int violationCap_;  // 0 = unbounded; once reached, stop growing U further

    Correspondence correspondence_;
    std::unordered_map<int, const regina::Simplex<2> *> skelIndexToSimplex_;
    std::vector<int> path_; // running U, in skeleton-index terms

    std::vector<Violation> *violations_;

  public:
    EmbeddednessAuditor(KnottedSurface &embedding,
                        const Skeleton<4, 2> &skeleton,
                        const std::vector<int> &graphToSkel, int maxDepth,
                        int violationCap, std::vector<Violation> &violations)
        : embedding_(embedding), skeleton_(skeleton),
          graphToSkel_(graphToSkel), maxDepth_(maxDepth),
          violationCap_(violationCap), violations_(&violations) {}

    // Once violationCap_ violations have been recorded, addFace() has
    // already been shown to accept bad states repeatedly -- further
    // exploration just re-derives the same finding at more sites without
    // adding diagnostic value, so this stops growing U at all (checked
    // BEFORE calling addFace(), so no state is committed and there's
    // nothing to roll back -- satisfies the transactional tryAdd contract).
    // This lets the enumerator unwind quickly instead of exhaustively
    // cataloguing every consequence of a bug already confirmed.
    bool tryAdd(int v) override {
        if (violationCap_ > 0 &&
            static_cast<int>(violations_->size()) >= violationCap_)
            return false;
        if (maxDepth_ > 0 && static_cast<int>(path_.size()) >= maxDepth_)
            return false;

        int f = graphToSkel_[v - 1];
        if (!embedding_.addFace(f))
            return false;

        const auto &tri = embedding_.triangulation();
        const auto *simplex = tri.simplex(tri.size() - 1);
        correspondence_[simplex] = skeleton_.getNodes()[f].face;
        skelIndexToSimplex_[f] = simplex;
        path_.push_back(f);

        // The real thing under test: does isEmbedded() -- addFace()'s O(1)
        // incremental tracking -- agree with the brute-force ground truth?
        // A disagreement in EITHER direction is a bug: isEmbedded()==true
        // while genuinely not embedded is unsound (callers would wrongly
        // trust a singular subcomplex); isEmbedded()==false while genuinely
        // embedded is incomplete (a real result would be silently dropped
        // by EmbeddingSearch's output gating).
        std::string why;
        bool genuinelyEmbedded = isGenuinelyEmbedded(correspondence_, &why);
        if (embedding_.isEmbedded() != genuinelyEmbedded) {
            std::ostringstream reason;
            reason << "isEmbedded() reports "
                   << (embedding_.isEmbedded() ? "true" : "false")
                   << " but ground truth is "
                   << (genuinelyEmbedded ? "true" : "false");
            if (!genuinelyEmbedded)
                reason << " (" << why << ")";
            violations_->push_back({path_, reason.str()});
        }

        // Likewise for isProper(), now tracked incrementally (see
        // bruteForceIsProper() above) instead of recomputed on every call.
        // Unlike the isEmbedded() check above, this one is unconditional:
        // badProperCount_ is driven purely by per-facet faceCount_
        // transitions and facetIsAmbientBoundary_, neither of which depends
        // on codimension>=2 injectivity, so it should hold even in
        // !isEmbedded() states.
        {
            bool properOK = bruteForceIsProper(correspondence_);
            if (embedding_.isProper() != properOK) {
                std::ostringstream reason;
                reason << "isProper() reports " << embedding_.isProper()
                       << " but brute-force ground truth is " << properOK;
                violations_->push_back({path_, reason.str()});
            }
        }

        return true;
    }

    void undo(int v) override {
        int f = graphToSkel_[v - 1];
        const auto *simplex = skelIndexToSimplex_.at(f);
        correspondence_.erase(simplex);
        skelIndexToSimplex_.erase(f);
        path_.pop_back();
        embedding_.removeFace(f);
    }

    const std::vector<int> &path() const { return path_; }
};

// Builds the 1-indexed DFS adjacency list directly from skeleton.getNodes(),
// mirroring EmbeddingSearch<4,2>::buildGraph_ (private, so not reusable
// directly) -- this is test-infra, not the thing under test, so a small
// reimplementation here doesn't compromise independence.
struct Graph {
    std::vector<std::vector<int>> adj; // 1-indexed; adj[0] unused
    std::vector<int> graphToSkel;      // graph vertex index -> skeleton index
};

Graph buildTestGraph(const Skeleton<4, 2> &skeleton) {
    const auto &nodes = skeleton.getNodes();

    std::vector<int> graphToSkel;
    std::vector<int> skelToGraph(nodes.size(), -1);
    for (size_t i = 0; i < nodes.size(); ++i) {
        if (EmbeddedSubmanifold<4, 2>::hasIrreparableSelfGluing(
                nodes[i].gluings)
            // hasUnexplainedSelfCollision() filters codimension >= 2
            // (vertex-level) self-collisions -- disabled along with Phase 2
            // in addFace() (embeddedsubmanifold.cpp), so it no longer
            // exists to call. Mirrors buildGraph_() in embeddingsearch.cpp,
            // which drops the same disjunct.
            //
            // || EmbeddedSubmanifold<4, 2>::hasUnexplainedSelfCollision(
            //        nodes[i].face, nodes[i].gluings)
        )
            continue;
        skelToGraph[i] = static_cast<int>(graphToSkel.size());
        graphToSkel.push_back(static_cast<int>(i));
    }

    int n = static_cast<int>(graphToSkel.size());
    std::vector<std::vector<int>> adj(n + 1);
    for (int graphIdx = 0; graphIdx < n; ++graphIdx) {
        int i = graphToSkel[graphIdx];
        std::set<int> neighbors;
        for (const auto &g : nodes[i].gluings) {
            int u = static_cast<int>(g.dstIndex);
            if (u == i)
                continue;
            int graphUdx = skelToGraph[u];
            if (graphUdx != -1)
                neighbors.insert(graphUdx);
        }
        for (int denseU : neighbors)
            adj[graphIdx + 1].push_back(denseU + 1);
    }
    return {std::move(adj), std::move(graphToSkel)};
}

struct AuditResult {
    std::vector<Violation> violations;
    long long subgraphsVisited = 0;
};

// Runs the audit over every connected induced subgraph of tri's 2-face
// gluing graph, up to maxDepth faces (0 = unbounded), stopping early once
// violationCap violations have been recorded (0 = unbounded -- only safe for
// the small hand-rolled cases). If checkBoundaryHomology is true,
// additionally checks -- at every connected state that is non-closed and
// properly/connectedly embedded -- that each boundary component's Link
// complement has the homology of a genuine n-component link complement
// (H_1 = Z^n). If checkOrientability is true, additionally cross-checks
// -- at every embedded state -- that embedding.isOrientable() (the
// incremental union-find-with-parity tracking in EmbeddedSubmanifold)
// agrees with embedding.triangulation().isOrientable() (Regina's own,
// independently implemented ground truth computed from scratch on the
// assembled subtri_), catching any discrepancy across every connected
// subcomplex actually reachable, not just one hand-picked example.
AuditResult auditAllEmbeddings(const regina::Triangulation<4> &tri,
                               int maxDepth, bool checkBoundaryHomology = false,
                               int violationCap = 25,
                               bool checkOrientability = false,
                               bool checkOrientedBoundary = false) {
    Skeleton<4, 2> skeleton(tri);
    Graph graph = buildTestGraph(skeleton);
    KnottedSurface embedding(skeleton);

    AuditResult result;
    EmbeddednessAuditor auditor(embedding, skeleton, graph.graphToSkel,
                                maxDepth, violationCap, result.violations);

    auto visit = [&](const std::vector<int> &) {
        ++result.subgraphsVisited;
        if (checkOrientability && embedding.isEmbedded()) {
            bool tracked = embedding.isOrientable();
            bool actual = embedding.triangulation().isOrientable();
            if (tracked != actual) {
                std::ostringstream reason;
                reason << "isOrientable() tracking mismatch: tracked="
                      << tracked << " actual (Triangulation<2>::isOrientable())="
                      << actual;
                result.violations.push_back({auditor.path(), reason.str()});
            }
        }
        if (!checkBoundaryHomology && !checkOrientedBoundary)
            return;
        if (!embedding.isEmbedded())
            return; // not genuinely embedded; boundaryLinks() presumes it is
        if (embedding.isClosed())
            return;
        if (!embedding.isProper() || !embedding.boundaryComponentsMapInjectively())
            return;

        if (checkBoundaryHomology) {
            for (const auto &[component, link] : embedding.boundaryLinks()) {
                regina::Triangulation<3> complement = link.buildComplement();
                regina::AbelianGroup h1 = complement.homology();
                if (!h1.isFree(link.countComponents())) {
                    std::ostringstream reason;
                    reason << "boundary component " << component
                          << ": link complement H_1 is not Z^"
                          << link.countComponents() << " (got " << h1.str()
                          << ")";
                    result.violations.push_back({auditor.path(), reason.str()});
                }
            }
        }

        // orientedBoundaryLinks()'s own precondition: only meaningful on an
        // orientable surface.
        if (checkOrientedBoundary && embedding.isOrientable()) {
            for (const auto &[component, curves] :
                embedding.orientedBoundaryLinks()) {
                for (const OrientedCurve &curve : curves) {
                    if (curve.empty()) {
                        result.violations.push_back(
                            {auditor.path(),
                            "orientedBoundaryLinks(): component " +
                                std::to_string(component) +
                                " produced an empty curve"});
                        continue;
                    }
                    // Independently re-verify chainIntoCurves()'s own
                    // invariant: each edge's head must be the next edge's
                    // tail, and the last edge's head must close back onto
                    // the first edge's own tail.
                    bool headToTailOk = true;
                    for (size_t i = 0; i < curve.size(); ++i) {
                        const OrientedEdge &e = curve[i];
                        const OrientedEdge &n = curve[(i + 1) % curve.size()];
                        const regina::Vertex<3> *head =
                            e.reversed ? e.edge->vertex(0) : e.edge->vertex(1);
                        const regina::Vertex<3> *nextTail =
                            n.reversed ? n.edge->vertex(1) : n.edge->vertex(0);
                        if (head != nextTail)
                            headToTailOk = false;
                    }
                    if (!headToTailOk) {
                        result.violations.push_back(
                            {auditor.path(),
                            "orientedBoundaryLinks(): component " +
                                std::to_string(component) +
                                " produced a curve that doesn't close "
                                "head-to-tail"});
                    }
                }

                // Cross-check against boundaryLinks()'s own (undirected)
                // edge count for the same component: orientedBoundaryLinks()
                // must partition the exact same edge set, just with
                // direction added, never dropping or duplicating an edge.
                size_t directedEdgeCount = 0;
                for (const OrientedCurve &curve : curves)
                    directedEdgeCount += curve.size();
                size_t undirectedEdgeCount = 0;
                for (const auto &[c2, link] : embedding.boundaryLinks()) {
                    if (c2 != component)
                        continue;
                    for (const auto &knot : link.comps_)
                        undirectedEdgeCount += knot.edgeIndices().size();
                }
                if (directedEdgeCount != undirectedEdgeCount) {
                    std::ostringstream reason;
                    reason << "orientedBoundaryLinks(): component "
                          << component << " has " << directedEdgeCount
                          << " directed edges but boundaryLinks() has "
                          << undirectedEdgeCount << " undirected edges";
                    result.violations.push_back({auditor.path(), reason.str()});
                }
            }
        }
    };

    ConnectedInducedSubgraphEnumerator enumerator(
        static_cast<int>(graph.graphToSkel.size()), graph.adj);
    enumerator.enumerateFiltered(visit, auditor);

    return result;
}

void reportAudit(const std::string &label, const AuditResult &result) {
    if (result.violations.empty()) {
        std::cout << green << "  PASS: " << resetColor << label
                  << " -- no embeddedness violations in " << result.subgraphsVisited
                  << " connected subsets\n";
        ++passed;
        return;
    }

    std::cout << red << "  FAIL: " << label << " -- " << result.violations.size()
              << " violation(s) found (showing up to 10)\n" << resetColor;
    size_t shown = std::min<size_t>(10, result.violations.size());
    for (size_t vi = 0; vi < shown; ++vi) {
        const auto &v = result.violations[vi];
        std::cout << "        faces [";
        for (size_t i = 0; i < v.faceSequence.size(); ++i) {
            std::cout << v.faceSequence[i];
            if (i + 1 < v.faceSequence.size())
                std::cout << ", ";
        }
        std::cout << "]: " << v.reason << "\n";
    }
    if (result.violations.size() > shown)
        std::cout << "        ... and " << (result.violations.size() - shown)
                  << " more\n";
    ++failed_count;
}

// ─────────────────────────────────────────────────────────────────────────────
// Section C: hand-rolled pentachora examples. A pentachoron has C(5,3) = 10
// triangular 2-faces (not 5 -- that's the facet/tetrahedron count), and
// these are naturally, densely interconnected via shared ambient edges even
// with NO gluings at all (they're all faces of one simplex) -- so the
// single-pentachoron gluing graph is already nontrivial on its own, before
// any join() calls. Two-pentachora cases get a modest depth cap: when
// addFace() is buggy, MORE states become reachable than a correct
// implementation would allow (the bug itself removes pruning), so an
// uncapped exhaustive search over ~20 faces can blow up combinatorially --
// observed directly during development (one case hung past 60s uncapped).
// The interesting violations here all show up within a handful of faces, so
// a small cap loses little coverage.
// ─────────────────────────────────────────────────────────────────────────────

constexpr int kHandRolledMaxDepth = 6;

// Single pentachoron, zero gluings: trivial baseline -- no bug should be
// findable here (the full pentachoron obviously embeds in itself).
void test_single_pentachoron_no_gluing() {
    std::cout << "\n--- Single pentachoron, no gluings ---\n";

    regina::Triangulation<4> tri;
    tri.newPentachoron();

    auto result = auditAllEmbeddings(tri, /*maxDepth=*/0);
    reportAudit("single pentachoron, no gluings", result);
}

// Single pentachoron, self-glued facet 0 <-> facet 1 via the transposition
// (1 0 2 3 4): this identifies ONLY local vertex 0 with local vertex 1 (vertex
// 1 lies on facet 0 and the gluing sends it to gluing[1] = 0 on facet 1),
// leaving vertices 2, 3, 4 completely untouched -- a minimal, surgical
// self-gluing. That vertex collision is NOT confined to facets 0/1 -- it
// also shows up on every OTHER triangular face containing both vertex 0 and
// vertex 1 (the faces {0,1,2}, {0,1,3}, {0,1,4}, each a subface of facets
// disjoint from the self-gluing itself). See
// test_single_face_internal_vertex_collision() below, which reuses this
// exact triangulation to target one of those faces directly.
regina::Triangulation<4> buildSelfGluedPentachoron() {
    regina::Triangulation<4> tri;
    auto *p = tri.newPentachoron();
    p->join(0, p, regina::Perm<5>(1, 0, 2, 3, 4));
    return tri;
}

void test_single_pentachoron_self_gluing() {
    std::cout << "\n--- Single pentachoron, one nontrivial self-gluing ---\n";

    regina::Triangulation<4> tri = buildSelfGluedPentachoron();
    auto result = auditAllEmbeddings(tri, /*maxDepth=*/0);
    reportAudit("single pentachoron, self-gluing", result);
}

// Two pentachora glued along one shared facet: plain simplicial baseline, no
// self-identifications anywhere.
void test_two_pentachora_one_shared_facet() {
    std::cout << "\n--- Two pentachora, one shared facet ---\n";

    regina::Triangulation<4> tri;
    auto *p = tri.newPentachoron();
    auto *q = tri.newPentachoron();
    p->join(4, q, regina::Perm<5>());

    auto result = auditAllEmbeddings(tri, kHandRolledMaxDepth);
    reportAudit("two pentachora, one shared facet", result);
}

// Two pentachora glued along one shared facet, PLUS a self-gluing on the
// other pentachoron's remaining facets -- combines the plain-simplicial and
// self-folded patterns in one ambient triangulation.
void test_two_pentachora_shared_facet_and_self_gluing() {
    std::cout << "\n--- Two pentachora: shared facet + self-gluing ---\n";

    regina::Triangulation<4> tri;
    auto *p = tri.newPentachoron();
    auto *q = tri.newPentachoron();
    p->join(4, q, regina::Perm<5>());
    q->join(0, q, regina::Perm<5>(1, 0, 2, 3, 4));

    auto result = auditAllEmbeddings(tri, kHandRolledMaxDepth);
    reportAudit("two pentachora, shared facet + self-gluing", result);
}

// Two pentachora glued to each other along TWO different facets
// simultaneously: the sharpest hand-rolled stress case for Condition 2 --
// two locally-disconnected attachment points between the same pair of
// simplices, giving multiple independent chances for a lower-dimensional
// face to be shared "accidentally" (i.e. without a single facet gluing
// directly identifying it).
void test_two_pentachora_two_shared_facets() {
    std::cout << "\n--- Two pentachora, two shared facets ---\n";

    regina::Triangulation<4> tri;
    auto *p = tri.newPentachoron();
    auto *q = tri.newPentachoron();
    // facet 3 of p <-> facet 4 of q
    p->join(3, q, regina::Perm<5>(0, 1, 2, 4, 3));
    // facet 4 of p <-> facet 3 of q (different destination facet on q, so
    // this doesn't collide with the gluing above)
    p->join(4, q, regina::Perm<5>(0, 1, 4, 2, 3));

    auto result = auditAllEmbeddings(tri, kHandRolledMaxDepth);
    reportAudit("two pentachora, two shared facets", result);
}

// regina::Example<4>::fourSphere(): two pentachora with ALL 5 facets
// identity-glued to each other -- forces essentially every corresponding
// vertex/edge/triangle of the two pentachora together. Free stress case,
// no hand-tuned Perm<5> needed.
void test_foursphere_doubled_simplex() {
    std::cout << "\n--- Example<4>::fourSphere() ---\n";

    auto tri = regina::Example<4>::fourSphere();
    EXPECT_EQ((int)tri.size(), 2, "fourSphere() is two pentachora");

    auto result = auditAllEmbeddings(tri, kHandRolledMaxDepth);
    reportAudit("fourSphere() doubled simplex", result);
}

// The Condition-2 blind spot: addFace()'s Condition 2 only ever fires for a
// k-face that was already touched by a PREVIOUSLY added face (faceCount_[k]
// != 0). It never checks whether a face's OWN vertices collide with each
// other.
//
// Two vertices of a triangle colliding usually comes bundled with an
// edge-level self-fold too (the edge between them becomes a loop, which
// hasIrreparableSelfGluing()/Condition 1 *does* catch) -- e.g. every
// degenerate face in buildSelfGluedPentachoron() above is excluded from the
// search graph for exactly that reason. But it doesn't have to: in the
// two-pentachora-shared-facet-and-self-gluing triangulation from
// test_two_pentachora_shared_facet_and_self_gluing() above, ambient face 6
// has vertices [0, 0, 3] (two of its own corners collide) yet its three
// edges remain three genuinely distinct Edge<4> objects -- verified
// directly by inspecting Skeleton<4,2>'s nodes during development. Nothing
// about this face looks locally self-folded, so it sails through both
// Condition 1 and Condition 2 unchallenged, and addFace() accepts it as a
// singleton -- even though it's degenerate on its own (its own two corners
// map to the same point). This function locates that face by scanning for
// the signature (colliding vertices, non-colliding edges) rather than
// hardcoding index 6, so it stays correct if the construction ever changes.
void test_single_face_internal_vertex_collision() {
    std::cout << "\n--- Single face with internally-colliding vertices, "
                 "non-colliding edges ---\n";

    regina::Triangulation<4> tri;
    auto *p = tri.newPentachoron();
    auto *q = tri.newPentachoron();
    p->join(4, q, regina::Perm<5>());
    q->join(0, q, regina::Perm<5>(1, 0, 2, 3, 4));

    Skeleton<4, 2> skeleton(tri);

    int target = -1;
    for (size_t i = 0; i < skeleton.numFaces(); ++i) {
        const auto *face = skeleton.getNodes()[i].face;
        bool vertexCollision = face->vertex(0) == face->vertex(1) ||
                               face->vertex(0) == face->vertex(2) ||
                               face->vertex(1) == face->vertex(2);
        bool edgeCollision = face->edge(0) == face->edge(1) ||
                             face->edge(0) == face->edge(2) ||
                             face->edge(1) == face->edge(2);
        if (vertexCollision && !edgeCollision) {
            target = static_cast<int>(i);
            break;
        }
    }
    EXPECT_EQ(target != -1, true,
              "found a triangular face with colliding vertices but "
              "non-colliding edges");
    if (target == -1)
        return;

    std::cout << "  target ambient face index: " << target << "\n";

    KnottedSurface embedding(skeleton);
    bool added = embedding.addFace(target);
    std::cout << "  addFace(" << target << ") returned "
              << (added ? "true" : "false") << "\n";

    if (!added) {
        std::cout << green << "  PASS: " << resetColor
                  << "addFace() correctly rejected the degenerate face\n";
        ++passed;
        return;
    }

    Correspondence correspondence;
    correspondence[embedding.triangulation().simplex(0)] =
        skeleton.getNodes()[target].face;
    std::string why;
    bool ok = isGenuinelyEmbedded(correspondence, &why);
    EXPECT_EQ(ok, false,
              "addFace() accepted a face with internally-colliding vertices, "
              "which isGenuinelyEmbedded correctly flags as not embedded "
              "(reason: " + (ok ? std::string("n/a") : why) + ")");
}

// Direct unit test of hasUnexplainedSelfCollision(), independent of the DFS:
// face 6 (no self-gluing entry) must be flagged, and face 7 (a legitimate
// self-gluing, facet 0 <-> 1) must NOT be -- pinning down both the fix and
// the no-over-rejection guarantee without going through addFace() or the
// search graph at all.
//
// Disabled along with hasUnexplainedSelfCollision() itself, which is
// commented out in embeddedsubmanifold.cpp (Phase 2 of addFace()). Kept
// here rather than deleted so it can be restored with the routine.
//
// void test_has_unexplained_self_collision() {
//     std::cout << "\n--- hasUnexplainedSelfCollision() unit check ---\n";
//
//     regina::Triangulation<4> tri;
//     auto *p = tri.newPentachoron();
//     auto *q = tri.newPentachoron();
//     p->join(4, q, regina::Perm<5>());
//     q->join(0, q, regina::Perm<5>(1, 0, 2, 3, 4));
//
//     Skeleton<4, 2> skeleton(tri);
//     const auto &nodes = skeleton.getNodes();
//
//     EXPECT_EQ((EmbeddedSubmanifold<4, 2>::hasUnexplainedSelfCollision(
//                   nodes[6].face, nodes[6].gluings)),
//               true, "face 6 (no self-gluing entry) IS flagged as an "
//                     "unexplained collision");
//     EXPECT_EQ((EmbeddedSubmanifold<4, 2>::hasUnexplainedSelfCollision(
//                   nodes[7].face, nodes[7].gluings)),
//               false, "face 7 (legitimate self-gluing facet 0<->1) is NOT "
//                      "flagged");
// }

// Documents that the incompleteness this test was originally written for is
// GONE. It used to record that face 6 was permanently excluded from the
// search graph by hasUnexplainedSelfCollision() (see
// ADDFACE_VERTEX_COLLISION_BUG.md), making the genuinely-embedded pair
// {6, 7} unreachable. That vertex-level filter has since been disabled along
// with Phase 2 of addFace(), so face 6 is reachable again and the pair is no
// longer lost. Both halves are asserted below, so this test now pins the
// restored behaviour rather than the old gap.
void test_buildgraph_known_incompleteness() {
    std::cout << "\n--- Known incompleteness: {6,7} unreachable despite "
                 "being genuinely embedded ---\n";

    regina::Triangulation<4> tri;
    auto *p = tri.newPentachoron();
    auto *q = tri.newPentachoron();
    p->join(4, q, regina::Perm<5>());
    q->join(0, q, regina::Perm<5>(1, 0, 2, 3, 4));

    Skeleton<4, 2> skeleton(tri);
    Graph graph = buildTestGraph(skeleton);

    bool face6Reachable = false;
    for (int skelIdx : graph.graphToSkel)
        if (skelIdx == 6)
            face6Reachable = true;
    EXPECT_EQ(face6Reachable, true,
              "face 6 is reachable again now that the vertex-level "
              "self-collision filter is disabled");

    // But {6, 7} really is genuinely embedded, outside the graph filter.
    KnottedSurface embedding(skeleton);
    bool added6 = embedding.addFace(6);
    bool added7 = embedding.addFace(7);
    Correspondence correspondence;
    correspondence[embedding.triangulation().simplex(0)] =
        skeleton.getNodes()[6].face;
    correspondence[embedding.triangulation().simplex(1)] =
        skeleton.getNodes()[7].face;
    bool ok = isGenuinelyEmbedded(correspondence);
    EXPECT_EQ(added6 && added7 && ok, true,
              "{6, 7} together is genuinely embedded, and is now reachable "
              "by the search again");
}

// The false->true direction of isEmbedded()'s tracking: face 6's own local
// vertices 0 and 1 collide at the same ambient point with no self-gluing
// entry explaining it (see test_single_face_internal_vertex_collision and
// test_has_unexplained_self_collision), so isEmbedded() must be false right
// after addFace(6) alone. test_buildgraph_known_incompleteness already
// established that {6, 7} together IS genuinely embedded (face 7's own
// legitimate self-gluing ends up identifying face 6's two colliding
// preimages in subtri_ once both are joined), so isEmbedded() must flip back
// to true once face 7 is added on top.
void test_isembedded_false_to_true_transition() {
    std::cout
        << "\n--- isEmbedded(): false->true transition (faces 6, 7) ---\n";

    regina::Triangulation<4> tri;
    auto *p = tri.newPentachoron();
    auto *q = tri.newPentachoron();
    p->join(4, q, regina::Perm<5>());
    q->join(0, q, regina::Perm<5>(1, 0, 2, 3, 4));

    Skeleton<4, 2> skeleton(tri);
    KnottedSurface embedding(skeleton);

    bool added6 = embedding.addFace(6);
    EXPECT_EQ(added6, true,
              "addFace(6) is accepted (Condition 1 doesn't catch this "
              "collision)");
    EXPECT_EQ(embedding.isEmbedded(), false,
              "isEmbedded() is false right after adding the self-colliding "
              "face 6 alone");

    bool added7 = embedding.addFace(7);
    EXPECT_EQ(added7, true, "addFace(7) is accepted");
    EXPECT_EQ(embedding.isEmbedded(), true,
              "isEmbedded() flips back to true once face 7's self-gluing "
              "identifies face 6's two colliding preimages in subtri_");
}

// LIFO stress test: addFace()/removeFace() must exactly restore
// isEmbedded()'s internal tracking (not just its boolean value) across
// repeated add/remove/re-add cycles, since removeFace()'s rollback replays
// the matching addFace()'s undo log rather than recomputing from scratch.
void test_isembedded_lifo_add_remove_add() {
    std::cout << "\n--- isEmbedded(): LIFO add/remove/re-add stress test "
                 "(faces 6, 7) ---\n";

    regina::Triangulation<4> tri;
    auto *p = tri.newPentachoron();
    auto *q = tri.newPentachoron();
    p->join(4, q, regina::Perm<5>());
    q->join(0, q, regina::Perm<5>(1, 0, 2, 3, 4));

    Skeleton<4, 2> skeleton(tri);
    KnottedSurface embedding(skeleton);

    EXPECT_EQ(embedding.isEmbedded(), true,
              "freshly constructed embedding is trivially embedded");

    embedding.addFace(6);
    EXPECT_EQ(embedding.isEmbedded(), false, "false after adding face 6 alone");

    embedding.addFace(7);
    EXPECT_EQ(embedding.isEmbedded(), true, "true again after adding face 7");
    size_t sizeAfterBoth = embedding.triangulation().size();

    embedding.removeFace(7);
    EXPECT_EQ(embedding.isEmbedded(), false,
              "removing face 7 restores the singularity it had resolved");

    // Re-add face 7: must reproduce the identical state, not just the same
    // boolean -- if removeFace()'s rollback under- or over-corrected the
    // DSU/registry state, a second addFace(7) could behave differently the
    // second time around (e.g. spuriously merge/not-merge).
    embedding.addFace(7);
    EXPECT_EQ(embedding.isEmbedded(), true, "re-adding face 7 reproduces true");
    EXPECT_EQ(embedding.triangulation().size(), sizeAfterBoth,
              "re-adding face 7 reproduces the same subtri_ size");

    embedding.removeFace(7);
    embedding.removeFace(6);
    EXPECT_EQ(embedding.isEmbedded(), true,
              "back to a pristine empty state after removing both faces");
    EXPECT_EQ(embedding.triangulation().size(), size_t{0},
              "subtri_ is empty again");
}

// ─────────────────────────────────────────────────────────────────────────────
// Section D: CobordismBuilder/knotbuilder-derived 4-manifolds. Depth-capped,
// but generously -- correctness coverage takes priority over speed here.
// ─────────────────────────────────────────────────────────────────────────────

// Calibrated empirically during development: the disc-derived 4-manifold
// below has ~64 candidate 2-faces with average graph degree ~12 (dense --
// each pentachoron contributes 10 triangles, heavily cross-linked by the
// cobordism's own layered structure). Connected-subgraph counts grow ~9x per
// extra depth level (445 at depth 2, 2.8M at depth 6), so a depth cap needs
// to stay modest here even though it can be generous for the much sparser
// hand-rolled pentachora cases above.
constexpr int kDerivedMaxDepth = 5;

void auditCobordismSurface(const std::string &label,
                           const regina::Triangulation<2> &surface,
                           int layers2, int layers3, bool cap) {
    std::cout << "\n--- " << label << " ---\n";

    if (!CobordismBuilder<2>::isOrdered(surface)) {
        std::cout << "  (skipped: base surface is not ordered, "
                     "CobordismBuilder<2> requires this)\n";
        return;
    }

    CobordismBuilder<2> cob2(surface);
    auto &tri3 = cob2.thicken(layers2);

    CobordismBuilder<3> cob3(tri3);
    if (layers3 > 0)
        cob3.thicken(layers3);
    if (cap)
        cob3.cone();

    auto result = auditAllEmbeddings(cob3.getCobordism(), kDerivedMaxDepth,
                                     /*checkBoundaryHomology=*/false,
                                     /*violationCap=*/25,
                                     /*checkOrientability=*/true);
    reportAudit(label, result);
}

void test_cobordism_disc() {
    auditCobordismSurface("CobordismBuilder: disc, thicken(1)x2",
                          regina::Example<2>::disc(), 1, 1, false);
}

void test_cobordism_mobius() {
    auditCobordismSurface("CobordismBuilder: mobius, thicken(1)x2",
                          regina::Example<2>::mobius(), 1, 1, false);
}

void test_cobordism_annulus() {
    auditCobordismSurface("CobordismBuilder: annulus, thicken(1)x2",
                          regina::Example<2>::annulus(), 1, 1, false);
}

// Example<2>::orientable(genus,punctures)/nonOrientable(...) are NOT ordered
// as constructed (verified directly -- CobordismBuilder<2>::isOrdered()
// returns false for every genus/puncture combination tried), and
// CobordismBuilder<2> can only auto-order dim-3 input (see
// CobordismBuilder::CobordismBuilder()'s `if constexpr (dim == 3)` branch);
// for dim 2 an unordered surface is a hard error. torus() and kb() are
// ordered as constructed and add genuinely different (closed, rather than
// punctured) topology to the coverage here.
void test_cobordism_torus() {
    auditCobordismSurface("CobordismBuilder: torus, thicken(1)x2",
                          regina::Example<2>::torus(), 1, 1, false);
}

void test_cobordism_kb() {
    auditCobordismSurface("CobordismBuilder: Klein bottle, thicken(1)x2",
                          regina::Example<2>::kb(), 1, 1, false);
}

// Calibrated empirically: knotbuilder's block-based construction produces
// noticeably larger, denser triangulations than the cobordism-only examples
// above (e.g. the trefoil's cone() has 237 candidate 2-faces at average
// graph degree ~18, vs. disc's 64 faces at degree ~12) -- connected-subgraph
// counts here grow ~15x per extra depth level rather than ~9x (2.3K at
// depth 2, 7.5M at depth 5), so this needs a noticeably smaller cap than
// kDerivedMaxDepth to stay fast.
constexpr int kKnotMaxDepth = 3;

void auditKnotSurface(const std::string &label, const std::string &pdCode,
                      bool reduce) {
    std::cout << "\n--- " << label << " ---\n";

    auto pd = knotbuilder::parsePDCode(pdCode);
    auto result0 = knotbuilder::buildLink(pd);

    knotbuilder::TriangulationWithLink result =
        reduce ? knotbuilder::reduceVertices(result0.tri, result0.edges)
              : result0;

    CobordismBuilder<3> cob(result.tri);
    cob.cone();

    auto audit = auditAllEmbeddings(cob.getCobordism(), kKnotMaxDepth,
                                    /*checkBoundaryHomology=*/true,
                                    /*violationCap=*/25,
                                    /*checkOrientability=*/false,
                                    /*checkOrientedBoundary=*/true);
    reportAudit(label, audit);
}

// PD codes reused verbatim from knotbuilder_test.cpp, where they're already
// validated (checked to produce a valid, closed triangulation of S^3 with
// the expected number of link components) before being used in any
// CobordismBuilder pipeline test there.
const char *kTrefoilPD = "1 4 2 5 3 6 4 1 5 2 6 3";  // 3_1, 1 component
const char *kHopfLinkPD = "1 4 2 3 3 2 4 1";         // 2-component link

void test_knotbuilder_trefoil_cone() {
    auditKnotSurface("knotbuilder: trefoil, cone()", kTrefoilPD, false);
}

void test_knotbuilder_trefoil_cone_reduced() {
    auditKnotSurface("knotbuilder: trefoil, reduceVertices() + cone()",
                     kTrefoilPD, true);
}

void test_knotbuilder_hopflink_cone() {
    auditKnotSurface("knotbuilder: Hopf link, cone()", kHopfLinkPD, false);
}

void test_knotbuilder_hopflink_cone_reduced() {
    auditKnotSurface("knotbuilder: Hopf link, reduceVertices() + cone()",
                     kHopfLinkPD, true);
}

// ─────────────────────────────────────────────────────────────────────────────
// Section D2: direct unit tests of EdgeComplement::linkingNumberWith().
//
// Nothing else in this suite reaches that method: the cone audits above cap
// at kKnotMaxDepth == 3 triangles, far too few to close a petal at the apex,
// so addFace()'s transverse-self-intersection branch is never entered. That
// is why linkingNumberWith() silently returned 0 for every link below --
// drilling curve A leaves an ideal (torus-cusp) vertex, and the homology it
// used to compute treated that vertex as an ordinary point, coning off the
// cusp instead of truncating it and killing the very meridian class the
// linking number is read from.
//
// These tests pin exact values rather than just "nonzero", so a future
// change that gets the sign convention or the truncation subtly wrong is
// caught rather than passing on a lucky !=0.
// ─────────────────────────────────────────────────────────────────────────────

// PD codes generated by Regina (regina::ExampleLink::torus(2,2n) / whitehead(),
// via Link::pdData()), then checked to build a valid closed triangulation with
// exactly 2 components.
const char *kT24PD = "1 5 2 8 3 7 4 6 5 3 6 2 7 1 8 4"; // (2,4) torus, lk = 2
const char *kT26PD = "1 7 2 12 3 9 4 8 5 11 6 10 7 3 8 2 9 5 10 4 11 1 12 6";
                                                        // (2,6) torus, lk = 3
const char *kWhiteheadPD = "1 7 2 8 3 10 4 9 5 2 6 3 8 5 9 4 10 6 7 1";
                                                        // lk = 0, NOT split

// `bothDirections` additionally checks lk(B,A), which exercises a genuinely
// different code path (each side drills its own curve and tracks the other
// through the drilling). It roughly doubles the cost, so it is enabled only
// for the cheapest link here -- HomologicalData on the unsimplified drilled
// complement is the dominant expense and grows steeply with triangulation
// size (the (2,6) case below is ~17s on its own).
void expectLinkingNumber(const char *label, const char *pdCode, long expected,
                         bool bothDirections = false) {
    auto result = knotbuilder::buildLink(knotbuilder::parsePDCode(pdCode));
    Link link(result.tri, result.edges);

    EXPECT_EQ(link.countComponents(), 2,
              std::string(label) + ": builds a 2-component link");
    if (link.countComponents() != 2)
        return;

    EXPECT_EQ(link.comps_[0].linkingNumberWith(link.comps_[1]), expected,
              std::string(label) + ": lk(A,B)");
    if (bothDirections)
        EXPECT_EQ(link.comps_[1].linkingNumberWith(link.comps_[0]), expected,
                  std::string(label) + ": lk(B,A) agrees");
}

void test_linking_number_hopf() {
    std::cout << "\n--- linkingNumberWith(): Hopf link, lk = 1 ---\n";
    expectLinkingNumber("Hopf", kHopfLinkPD, 1, /*bothDirections=*/true);
}

void test_linking_number_torus24() {
    std::cout << "\n--- linkingNumberWith(): (2,4) torus link, lk = 2 ---\n";
    expectLinkingNumber("(2,4) torus", kT24PD, 2);
}

void test_linking_number_torus26() {
    std::cout << "\n--- linkingNumberWith(): (2,6) torus link, lk = 3 ---\n";
    expectLinkingNumber("(2,6) torus", kT26PD, 3);
}

// The Whitehead link is the case that distinguishes "computes the linking
// number" from "detects whether the pair is separable": its two components
// are genuinely, inseparably linked, yet lk = 0. A correct linkingNumberWith()
// must report 0 here -- under-reporting is what keeps addFace()'s rejection
// sound (it never prunes a configuration it hasn't proved is transverse).
void test_linking_number_whitehead() {
    std::cout << "\n--- linkingNumberWith(): Whitehead link, lk = 0 "
                 "(linked, but algebraically zero) ---\n";
    expectLinkingNumber("Whitehead", kWhiteheadPD, 0);
}

// ─────────────────────────────────────────────────────────────────────────────
// Section E: end-to-end regression test for the originally reported bug --
// coning a triangulated S^3 traced by a non-slice knot to a point used to
// produce an "embedded disk" bounded by that knot: topologically a disk
// (cone on any knot is one), but not locally flat at the cone point, since
// the disk's intersection with the apex's link is the knotted curve itself,
// not an unknot. addFace() must now reject the last cone-triangle that
// would close the apex's petal, rather than silently accepting the disk.
// ─────────────────────────────────────────────────────────────────────────────

// Reconstructs the "coning surface" bounded by `edges` in coned -- one
// triangle {apex, i, j} per edge {i, j}, in whichever original tetrahedron
// each edge fronts -- and returns their ambient Skeleton<4,2> face indices,
// in the same order as `edges`. Relies on CobordismBuilder<3>::cone()'s own
// construction (see cobordismbuilder.cpp): one new pentachoron per original
// tetrahedron, in matching index order, with local vertices 0..3 preserved
// verbatim from the original tetrahedron and local vertex 4 always the
// single shared apex.
std::vector<int> coneTriangleIndices(const Skeleton<4, 2> &skeleton,
                                     const regina::Triangulation<4> &coned,
                                     const std::vector<const regina::Edge<3> *>
                                         &edges) {
    std::unordered_map<const regina::Face<4, 2> *, int> faceToSkelIdx;
    for (size_t i = 0; i < skeleton.numFaces(); ++i)
        faceToSkelIdx[skeleton.getNodes()[i].face] = static_cast<int>(i);

    std::vector<int> result;
    result.reserve(edges.size());
    for (const regina::Edge<3> *e : edges) {
        auto emb = e->front();
        auto *pent = coned.simplex(emb.tetrahedron()->index());
        regina::Perm<4> ordering =
            regina::FaceNumbering<3, 1>::ordering(emb.edge());
        int lv0 = ordering[0];
        int lv1 = ordering[1];
        int t = regina::FaceNumbering<4, 2>::triangleNumber[lv0][lv1][4];
        const regina::Face<4, 2> *triangle = pent->triangle(t);
        result.push_back(faceToSkelIdx.at(triangle));
    }
    return result;
}

void test_cone_on_trefoil_rejected() {
    std::cout << "\n--- Coning the trefoil to a point: apex closure is "
                 "rejected (non-locally-flat) ---\n";

    auto pd = knotbuilder::parsePDCode(kTrefoilPD);
    auto result = knotbuilder::buildLink(pd);

    CobordismBuilder<3> cob(result.tri);
    auto &coned = cob.cone();

    Skeleton<4, 2> skeleton(coned);
    KnottedSurface embedding(skeleton);

    std::vector<int> coneFaces =
        coneTriangleIndices(skeleton, coned, result.edges);
    EXPECT_EQ(coneFaces.empty(), false,
              "the trefoil's traced edges give at least one cone-triangle");

    bool sawRejection = false;
    for (int f : coneFaces) {
        if (!embedding.addFace(f)) {
            sawRejection = true;
            break;
        }
    }
    EXPECT_EQ(sawRejection, true,
              "addFace() rejects one of the cone-triangles that would "
              "close the apex's petal into the (knotted) trefoil curve");
}

// Same construction, but with the unknot: coning a genuine unknot to a
// point IS locally flat everywhere (the cone point's link intersection is
// an honest unknot), so every cone-triangle must be accepted.
void test_cone_on_unknot_accepted() {
    std::cout << "\n--- Coning the unknot to a point: apex closure is "
                 "accepted (locally flat) ---\n";

    // A single-crossing unknot diagram (a Reidemeister-1 kink), reusing
    // knotbuilder's own PD-code conventions.
    auto pd = knotbuilder::parsePDCode("1 2 2 1");
    auto result = knotbuilder::buildLink(pd);

    CobordismBuilder<3> cob(result.tri);
    auto &coned = cob.cone();

    Skeleton<4, 2> skeleton(coned);
    KnottedSurface embedding(skeleton);

    std::vector<int> coneFaces =
        coneTriangleIndices(skeleton, coned, result.edges);

    bool allAccepted = true;
    for (int f : coneFaces) {
        if (!embedding.addFace(f)) {
            allAccepted = false;
            break;
        }
    }
    EXPECT_EQ(allAccepted, true,
              "every cone-triangle over the unknot is accepted -- the apex "
              "is genuinely locally flat");
}

// ─────────────────────────────────────────────────────────────────────────────
// Section E2: unlinked self-intersections (pl_enumeration_draft §4.5,
// KnottedSurface::isResolvable()/isAcceptable()).
//
// Coning a link L in S^3 to the interior apex of B^4 = cone(S^3) gives one
// disc per component of L, all meeting at the apex: one closed petal per
// component, and the trace T_v(S) there is L itself. So which L is coned decides
// every case: an unlink is resolvable, a Hopf link is pruned by
// P_transverse, and a Whitehead link passes every prune yet is not an
// unlink -- the case only the certificate can refuse.
// ─────────────────────────────────────────────────────────────────────────────

// 0-indexed or Regina (1-indexed) PD codes, each through parsePDCode().
// kUnlink3PD is Link(3) with two R2 moves; its group was checked to be
// Free(3) in Regina, independently of anything here.
const char *kUnlink2PD = "0 3 1 2 1 3 0 2";
const char *kUnlink3PD = "3 1 4 2 6 1 3 2 7 5 8 6 8 5 7 4";

// Adds every cone-triangle over `pdCode`'s link to `embedding`, reporting
// whether addFace() accepted all of them.
bool addConeOver(KnottedSurface &embedding, const Skeleton<4, 2> &skeleton,
                 const regina::Triangulation<4> &coned,
                 const std::vector<const regina::Edge<3> *> &edges) {
    for (int f : coneTriangleIndices(skeleton, coned, edges))
        if (!embedding.addFace(f))
            return false;
    return true;
}

void expectConeOnUnlinkResolvable(const char *label, const char *pdCode,
                                  size_t components) {
    std::cout << "\n--- Coning the " << label
              << ": a resolvable self-intersection at the apex ---\n";
    auto result = knotbuilder::buildLink(knotbuilder::parsePDCode(pdCode));
    CobordismBuilder<3> cob(result.tri);
    auto &coned = cob.cone();
    Skeleton<4, 2> skeleton(coned);

    PetalCache cache;
    SelfIntersectionCensus census;
    KnottedSurface off(skeleton);
    KnottedSurface on({.resolveUnlinked = true, .census = &census}, skeleton,
                      cache);

    std::string l(label);
    EXPECT_EQ(addConeOver(off, skeleton, coned, result.edges), true,
              l + ": every cone-triangle is accepted by addFace() -- the "
                  "petals are unknots with pairwise linking number 0");
    EXPECT_EQ(addConeOver(on, skeleton, coned, result.edges), true,
              l + ": ...with the options set, too");

    EXPECT_EQ(off.satisfies(BoundaryCondition::proper), true,
              l + ": the cone is properly embedded (boundary = the link)");
    EXPECT_EQ(off.isEmbedded(), false,
              l + ": it is NOT 2-embedded: " + std::to_string(components) +
                  " petals meet at the apex");
    EXPECT_EQ(off.singularVertexCount(), size_t{1},
              l + ": ...at exactly one vertex, the apex");
    EXPECT_EQ(off.isResolvable(), true,
              l + ": the apex's trace T_v(S) is certified the unlink");
    EXPECT_EQ(off.isAcceptable(), false,
              l + ": without --resolve-unlinked it is still rejected");
    EXPECT_EQ(on.isAcceptable(), true,
              l + ": with it, the surface is accepted");

    auto tubed = KnottedSurface::tubedSurfaceType(on.triangulation());
    EXPECT_EQ(on.triangulation().countComponents(), components,
              l + ": the abstract surface Delta[S] is the " +
                  std::to_string(components) +
                  " discs apart -- what the resolution embeds");
    EXPECT_EQ(tubed.genus, 0, l + ": ...of total genus 0");
    EXPECT_EQ(tubed.punctures, static_cast<int>(components),
              l + ": ...with one boundary circle per component");

    on.tallySelfIntersection();
    EXPECT_EQ(census.singular.load(), 1LL, l + ": census counts it singular");
    EXPECT_EQ(census.interiorUnlinked.load(), 1LL,
              l + ": ...as interior and certified unlinked");
    const long long checkedBefore = census.audited.load();
    EXPECT_EQ(on.isSmoothAtBoundary(), true,
              l + ": its boundary-vertex petals (cones on two points) are "
                  "unknotted, so it is smooth at the boundary");
    EXPECT_EQ(census.audited.load(), checkedBefore + 1,
              l + ": the boundary check is recorded in the census");
    EXPECT_EQ(census.auditKnotted.load(), 0LL,
              l + ": ...and never rejected it (isAcceptable() ran it too)");
}

void test_cone_on_unlink2_resolvable() {
    expectConeOnUnlinkResolvable("2-component unlink", kUnlink2PD, 2);
}

void test_cone_on_unlink3_resolvable() {
    expectConeOnUnlinkResolvable("3-component unlink", kUnlink3PD, 3);
}

void test_cone_on_whitehead_not_resolvable() {
    std::cout << "\n--- Coning the Whitehead link: passes every prune, but "
                 "is not resolvable ---\n";
    auto result = knotbuilder::buildLink(knotbuilder::parsePDCode(kWhiteheadPD));
    CobordismBuilder<3> cob(result.tri);
    auto &coned = cob.cone();
    Skeleton<4, 2> skeleton(coned);

    PetalCache cache;
    SelfIntersectionCensus census;
    KnottedSurface on({.resolveUnlinked = true, .census = &census}, skeleton,
                      cache);
    EXPECT_EQ(addConeOver(on, skeleton, coned, result.edges), true,
              "addFace() accepts every cone-triangle: unknotted components, "
              "linking number 0, so neither P_smooth nor P_transverse fires");
    EXPECT_EQ(on.isEmbedded(), false, "two petals meet at the apex");
    EXPECT_EQ(on.isResolvable(), false,
              "the Whitehead link is not certified an unlink");
    EXPECT_EQ(on.isAcceptable(), false,
              "so even with --resolve-unlinked the surface is rejected");
    on.tallySelfIntersection();
    EXPECT_EQ(census.interiorUncertified.load(), 1LL,
              "census: interior, not certified");
}

void test_cone_on_hopf_still_pruned() {
    std::cout << "\n--- Coning the Hopf link: still pruned by "
                 "P_transverse ---\n";
    auto result = knotbuilder::buildLink(knotbuilder::parsePDCode(kHopfLinkPD));
    CobordismBuilder<3> cob(result.tri);
    auto &coned = cob.cone();
    Skeleton<4, 2> skeleton(coned);

    PetalCache cache;
    KnottedSurface on({.resolveUnlinked = true}, skeleton, cache);
    EXPECT_EQ(addConeOver(on, skeleton, coned, result.edges), false,
              "addFace() rejects the triangle closing the second petal "
              "(linking number 1), whatever the options");
}

void test_boundary_vertex_self_intersection_not_resolvable() {
    std::cout << "\n--- Two petals at a BOUNDARY vertex are never "
                 "resolvable ---\n";
    auto result = knotbuilder::buildLink(knotbuilder::parsePDCode(kUnlink2PD));
    CobordismBuilder<3> cob(result.tri);
    auto &coned = cob.cone();
    Skeleton<4, 2> skeleton(coned);

    // Two boundary triangles sharing a vertex but no edge: each is its own
    // petal there, so they meet at a vertex on the boundary -- where the
    // boundary curves themselves would touch, which no perturbation rel
    // boundary can undo.
    std::optional<std::pair<int, int>> pair;
    for (size_t a = 0; a < skeleton.numFaces() && !pair; ++a) {
        const auto *fa = skeleton.getNodes()[a].face;
        if (!fa->isBoundary())
            continue;
        for (size_t b = a + 1; b < skeleton.numFaces() && !pair; ++b) {
            const auto *fb = skeleton.getNodes()[b].face;
            if (!fb->isBoundary())
                continue;
            int sharedVertices = 0;
            bool sharedEdge = false;
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    if (fa->vertex(i) == fb->vertex(j))
                        ++sharedVertices;
                    if (fa->edge(i) == fb->edge(j))
                        sharedEdge = true;
                }
            }
            if (sharedVertices == 1 && !sharedEdge)
                pair.emplace(static_cast<int>(a), static_cast<int>(b));
        }
    }
    EXPECT_EQ(pair.has_value(), true,
              "fixture: two boundary triangles meeting in exactly one vertex");
    if (!pair)
        return;

    PetalCache cache;
    KnottedSurface on({.resolveUnlinked = true}, skeleton, cache);
    EXPECT_EQ(on.addFace(pair->first) && on.addFace(pair->second), true,
              "both triangles are added");
    EXPECT_EQ(on.isEmbedded(), false, "they meet at their shared vertex");
    EXPECT_EQ(on.isResolvable(), false,
              "a self-intersection at a boundary vertex is not resolvable");
    EXPECT_EQ(on.isAcceptable(), false,
              "so it is rejected even with --resolve-unlinked");

    // The census's multi-open split, by where the vertex lies relative to
    // the search side: -1 (unknown) records neither half; its own boundary
    // component makes it search-side; any other makes it far-side, where
    // this lone vertex with exactly two open petals is the "simple" case.
    long bc = -1;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            if (skeleton.getNodes()[pair->first].face->vertex(i) ==
                skeleton.getNodes()[pair->second].face->vertex(j))
                bc = static_cast<long>(skeleton.getNodes()[pair->first]
                                           .face->vertex(i)
                                           ->boundaryComponent()
                                           ->index());
    for (long side : {-1L, bc, bc + 1}) {
        SelfIntersectionCensus census;
        census.searchSideBoundary = side;
        PetalCache sideCache;
        KnottedSurface tallied({.census = &census}, skeleton, sideCache);
        tallied.addFace(pair->first);
        tallied.addFace(pair->second);
        tallied.tallySelfIntersection();
        const std::string l = "search side " + std::to_string(side) + ": ";
        const bool far = side == bc + 1;
        EXPECT_EQ(census.multiOpen.load(), 1LL, l + "multi-open");
        EXPECT_EQ(census.multiOpenSearchSide.load(), side == bc ? 1LL : 0LL,
                  l + "search-side split");
        EXPECT_EQ(census.multiOpenFar.load(), far ? 1LL : 0LL,
                  l + "far-side split");
        EXPECT_EQ(census.multiOpenFarClean.load(), far ? 1LL : 0LL,
                  l + "no other singular vertex, so clean");
        EXPECT_EQ(census.multiOpenFarSimple.load(), far ? 1LL : 0LL,
                  l + "one vertex, two open petals: simple");
        EXPECT_EQ(census.farConfigs.size(), size_t(far ? 1 : 0),
                  l + "one far-side configuration");
    }
}

// Smoothness at the boundary (pl_enumeration_draft §4.3): the post hoc check
// isAcceptable() applies at boundary vertices, which addFace() never checks.
//
// Ambient: the cone on a 3-ball B, whose apex lies on the boundary of the
// resulting 4-ball. B is knotbuilder's S^3 for a knot K minus one tetrahedron
// containing an edge e of K, so K \ e is a proper arc in B, knotted exactly
// when K is. The cone on that arc is proper (the arc lies in B, part of the
// boundary, and its end cone-edges in the cone on dB), embedded, and passes
// every addFace() check, since all its vertices are boundary vertices. At the
// apex it has one open petal whose trace is the arc. So the boundary check is
// the only thing that can reject it -- and must, for the trefoil.
void test_boundary_filter_rejects_knotted_boundary_petal() {
    std::cout << "\n--- A knotted petal at a boundary vertex is rejected by "
                 "the boundary check ---\n";
    struct Case {
        const char *label;
        const char *pd;
        bool smooth;
    };
    for (const Case &c : {Case{"R1-kinked unknot", "1 2 2 1", true},
                          Case{"trefoil", kTrefoilPD, false}}) {
        std::string l(c.label);
        auto result = knotbuilder::buildLink(knotbuilder::parsePDCode(c.pd));

        // Built and modified in place, never moved: face pointers into it
        // are taken below, and each arc edge is recorded positionally first
        // because removeSimplex() rebuilds the skeleton.
        std::optional<regina::Triangulation<3>> ball;
        std::vector<const regina::Edge<3> *> arc;
        for (const auto *cut : result.edges) {
            for (const auto &emb : *cut) {
                const auto *tet = emb.tetrahedron();
                bool distinct = true;
                for (int i = 0; i < 4 && distinct; ++i)
                    for (int j = i + 1; j < 4 && distinct; ++j)
                        distinct = tet->vertex(i) != tet->vertex(j);
                if (!distinct)
                    continue;
                ball.emplace(result.tri);
                std::vector<std::pair<regina::Tetrahedron<3> *, int>> desc;
                bool ok = true;
                for (const auto *k : result.edges) {
                    if (k == cut)
                        continue;
                    bool found = false;
                    for (const auto &e2 : *k)
                        if (e2.tetrahedron()->index() != tet->index()) {
                            desc.emplace_back(ball->tetrahedron(
                                                  e2.tetrahedron()->index()),
                                              e2.edge());
                            found = true;
                            break;
                        }
                    ok = ok && found;
                }
                if (!ok) {
                    ball.reset();
                    continue;
                }
                ball->removeSimplex(ball->tetrahedron(tet->index()));
                if (!ball->isBall()) {
                    ball.reset();
                    continue;
                }
                for (auto [t, i] : desc)
                    arc.push_back(t->edge(i));
                break;
            }
            if (ball)
                break;
        }
        EXPECT_EQ(ball.has_value(), true,
                  l + ": fixture, a ball containing the knot minus one edge");
        if (!ball)
            continue;

        CobordismBuilder<3> cob(*ball);
        auto &coned = cob.cone();
        Skeleton<4, 2> skeleton(coned);
        PetalCache cache;
        SelfIntersectionCensus census;
        KnottedSurface s({.census = &census}, skeleton, cache);
        EXPECT_EQ(addConeOver(s, skeleton, coned, arc), true,
                  l + ": addFace() accepts every cone-triangle (all its "
                      "vertices are boundary vertices, which it never checks)");
        EXPECT_EQ(s.satisfies(BoundaryCondition::proper), true,
                  l + ": the cone on the arc is proper");
        EXPECT_EQ(s.isEmbedded(), true, l + ": and embedded");
        EXPECT_EQ(s.isSmoothAtBoundary(), c.smooth,
                  l + (c.smooth ? ": its apex petal is unknotted"
                                : ": its apex petal (the knotted arc) is "
                                  "knotted"));
        EXPECT_EQ(s.isAcceptable(), c.smooth,
                  l + (c.smooth ? ": so it is accepted"
                                : ": so it is rejected, embedded or not"));
        EXPECT_EQ(census.auditKnotted.load(), c.smooth ? 0LL : 2LL,
                  l + ": the census records each rejection");
    }
}

// The case the boundary guard in vertexUnlinked_() exists for, which the test
// above cannot reach (its petals are open, so the closedness check already
// refuses them): CLOSED petals whose traces form an unlink, at a BOUNDARY
// vertex.
//
// Ambient: the cone on a 3-ball B -- knotbuilder's S^3 for the 2-component
// unlink, punctured (Triangulation<3>::puncture(): a triangle is thickened
// into a pillow and a ball removed from inside it, so the new boundary sphere
// has only new vertices and the link stays in the interior). The cone on a
// ball is a 4-ball whose apex lies on its boundary (in the cone on dB). The
// cone on the link then has one CLOSED petal per component at the apex (its
// spokes run into the interior, since the link misses dB), and the petal link
// is the unlink. Everything the certificate checks passes; only the apex being
// on the boundary should stop it.
void test_boundary_vertex_closed_unlinked_petals_not_resolvable() {
    std::cout << "\n--- Closed, unlinked petals at a BOUNDARY vertex are "
                 "still not resolvable ---\n";
    auto result = knotbuilder::buildLink(knotbuilder::parsePDCode(kUnlink2PD));

    // Built and modified in place, never moved: face pointers into it are
    // taken below. puncture() rebuilds the skeleton but keeps every existing
    // tetrahedron, so each link edge is recorded as (tetrahedron, local edge)
    // first and re-resolved afterwards.
    std::optional<regina::Triangulation<3>> ball;
    ball.emplace(result.tri);
    std::vector<std::pair<regina::Tetrahedron<3> *, int>> desc;
    for (const auto *e : result.edges) {
        const auto &emb = e->front();
        desc.emplace_back(ball->tetrahedron(emb.tetrahedron()->index()),
                          emb.edge());
    }
    ball->puncture();
    std::vector<const regina::Edge<3> *> linkEdges;
    for (auto [t, i] : desc)
        linkEdges.push_back(t->edge(i));
    bool linkInterior = true;
    for (const auto *e : linkEdges)
        linkInterior = linkInterior && !e->vertex(0)->isBoundary() &&
                       !e->vertex(1)->isBoundary();
    EXPECT_EQ(ball->isBall(), true, "fixture: the punctured S^3 is a ball");
    EXPECT_EQ(linkInterior, true,
              "fixture: the link misses the new boundary sphere");
    if (!ball->isBall() || !linkInterior)
        ball.reset();
    if (!ball)
        return;

    CobordismBuilder<3> cob(*ball);
    auto &coned = cob.cone();
    Skeleton<4, 2> skeleton(coned);

    PetalCache cache;
    SelfIntersectionCensus census;
    KnottedSurface on({.resolveUnlinked = true, .census = &census}, skeleton,
                      cache);
    EXPECT_EQ(addConeOver(on, skeleton, coned, linkEdges), true,
              "every cone-triangle is added");
    EXPECT_EQ(on.satisfies(BoundaryCondition::proper), true,
              "the cone on the link is proper (its boundary, the link, "
              "lies in B, part of the boundary)");
    EXPECT_EQ(on.singularVertexCount(), size_t{1},
              "it meets itself at exactly one vertex, the apex");
    EXPECT_EQ(on.isResolvable(), false,
              "closed, unlinked petals at a boundary vertex are refused");
    EXPECT_EQ(on.isAcceptable(), false,
              "so even --resolve-unlinked rejects it");
    on.tallySelfIntersection();
    EXPECT_EQ(census.boundaryUnlinked.load(), 1LL,
              "census: the capped boundary test DOES certify it -- so the "
              "refusal above is the boundary guard's doing, not the "
              "certificate's");
}

// ─────────────────────────────────────────────────────────────────────────────
// Section F: hereditariness stress test for the new checks. enumerate_cis.h's
// filtered DFS *requires* the predicate be genuinely hereditary (for every
// connected U* satisfying it, every connected subset of U* must too) --
// a violation would silently drop results rather than just misclassify one
// state, so this gets its own dedicated check: for every face accepted while
// building up a random walk, removing it (if the remainder stays connected)
// and re-adding it via the same overridden addFace() must also succeed.
// ─────────────────────────────────────────────────────────────────────────────

void test_hereditariness_stress() {
    std::cout << "\n--- Hereditariness stress test: no accepted face's "
                 "removal+re-add ever newly fails ---\n";

    // Deliberately a small triangulation (not the trefoil's dense
    // reduceVertices()+cone() graph, which is fast to audit exhaustively
    // but too expensive to also replay through O(steps^2) addFace()/
    // removeFace() calls here) -- a single-crossing unknot's cone still
    // has a genuinely closing/reopening petal at the apex, which is what
    // this stress test needs to exercise, without the cost of repeatedly
    // recognizing a much larger knot complement.
    auto pd = knotbuilder::parsePDCode("1 2 2 1");
    auto result = knotbuilder::buildLink(pd);

    CobordismBuilder<3> cob(result.tri);
    auto &coned = cob.cone();

    Skeleton<4, 2> skeleton(coned);
    Graph graph = buildTestGraph(skeleton);
    KnottedSurface embedding(skeleton);

    std::mt19937 rng(12345);
    std::vector<int> path; // skeleton face indices, in addition order

    auto isStillConnected = [&](const std::vector<int> &faces) {
        if (faces.size() <= 1)
            return true;
        std::unordered_set<int> faceSet(faces.begin(), faces.end());
        std::unordered_set<int> reached{faces[0]};
        std::vector<int> stack{faces[0]};
        while (!stack.empty()) {
            int f = stack.back();
            stack.pop_back();
            for (const auto &g : skeleton.getNodes()[f].gluings) {
                int u = static_cast<int>(g.dstIndex);
                if (faceSet.contains(u) && reached.insert(u).second)
                    stack.push_back(u);
            }
        }
        return reached.size() == faces.size();
    };

    int steps = 0, removalChecks = 0, attempts = 0;
    while (steps < 20 && attempts < 2000) {
        ++attempts;
        // Try a random not-yet-present face adjacent to the current path
        // (or any face, if path is empty); addFace()'s own Phase 1 rejects
        // anything not actually addable, so a failed attempt just retries.
        int candidate;
        if (path.empty()) {
            candidate =
                static_cast<int>(rng() % skeleton.numFaces());
        } else {
            int anchor = path[rng() % path.size()];
            const auto &gluings = skeleton.getNodes()[anchor].gluings;
            if (gluings.empty())
                break;
            candidate =
                static_cast<int>(gluings[rng() % gluings.size()].dstIndex);
        }
        if (std::ranges::find(path, candidate) != path.end())
            continue;
        if (!embedding.addFace(candidate))
            continue;
        path.push_back(candidate);
        ++steps;

        // Pick a random already-added face and check heredity: if removing
        // it leaves the rest connected, re-adding it (after removing it)
        // must succeed, since it was already known addable in a superset.
        int idx = static_cast<int>(rng() % path.size());
        int victim = path[idx];
        std::vector<int> withoutVictim = path;
        withoutVictim.erase(withoutVictim.begin() + idx);
        if (!isStillConnected(withoutVictim))
            continue;

        // Roll back to just before victim was added (LIFO), removing
        // everything added after it too, then replay everything except
        // victim, then attempt to add victim back on top.
        std::vector<int> after(path.begin() + idx + 1, path.end());
        for (auto it = after.rbegin(); it != after.rend(); ++it)
            embedding.removeFace(*it);
        embedding.removeFace(victim);

        ++removalChecks;
        bool reAdded = embedding.addFace(victim);
        EXPECT_EQ(reAdded, true,
                  "removing an already-accepted face (leaving the rest "
                  "connected) and re-adding it never newly fails");

        if (!reAdded) {
            // Restore path to a consistent (victim absent) state so the
            // loop can keep going without redoing all the bookkeeping.
            path.erase(path.begin() + idx);
            for (int f : after) {
                embedding.addFace(f);
                path.push_back(f);
            }
            continue;
        }

        for (int f : after)
            embedding.addFace(f);
    }

    EXPECT_EQ(removalChecks > 0, true,
              "at least one removal+re-add heredity check actually ran");
    std::cout << "  (" << steps << " random additions, " << removalChecks
              << " heredity checks)\n";
}

// ─────────────────────────────────────────────────────────────────────────────

void run(const std::string &name, void (*fn)()) {
    std::cout << bold << "\n=== " << name << " ===" << resetColor << "\n";
    std::cout.flush();
    auto start = std::chrono::steady_clock::now();
    try {
        fn();
    } catch (const std::exception &e) {
        std::cout << red << "  EXCEPTION: " << e.what() << resetColor << "\n";
        ++failed_count;
    }
    double secs = std::chrono::duration<double>(
                      std::chrono::steady_clock::now() - start)
                      .count();
    std::cout << "  (" << secs << "s)\n";
    std::cout.flush();
}

int main() {
    run("single_pentachoron_no_gluing", test_single_pentachoron_no_gluing);
    run("single_pentachoron_self_gluing", test_single_pentachoron_self_gluing);
    run("two_pentachora_one_shared_facet", test_two_pentachora_one_shared_facet);
    run("two_pentachora_shared_facet_and_self_gluing",
        test_two_pentachora_shared_facet_and_self_gluing);
    run("two_pentachora_two_shared_facets", test_two_pentachora_two_shared_facets);
    run("foursphere_doubled_simplex", test_foursphere_doubled_simplex);
    run("single_face_internal_vertex_collision",
        test_single_face_internal_vertex_collision);
    // Disabled with hasUnexplainedSelfCollision() itself; see above.
    // run("has_unexplained_self_collision", test_has_unexplained_self_collision);
    run("buildgraph_known_incompleteness", test_buildgraph_known_incompleteness);
    run("isembedded_false_to_true_transition",
        test_isembedded_false_to_true_transition);
    run("isembedded_lifo_add_remove_add", test_isembedded_lifo_add_remove_add);

    run("cobordism_disc", test_cobordism_disc);
    run("cobordism_mobius", test_cobordism_mobius);
    run("cobordism_annulus", test_cobordism_annulus);
    run("cobordism_torus", test_cobordism_torus);
    run("cobordism_kb", test_cobordism_kb);

    run("knotbuilder_trefoil_cone", test_knotbuilder_trefoil_cone);
    run("knotbuilder_trefoil_cone_reduced", test_knotbuilder_trefoil_cone_reduced);
    run("knotbuilder_hopflink_cone", test_knotbuilder_hopflink_cone);
    run("knotbuilder_hopflink_cone_reduced", test_knotbuilder_hopflink_cone_reduced);
    run("linking_number_hopf", test_linking_number_hopf);
    run("linking_number_torus24", test_linking_number_torus24);
    run("linking_number_torus26", test_linking_number_torus26);
    run("linking_number_whitehead", test_linking_number_whitehead);

    run("cone_on_trefoil_rejected", test_cone_on_trefoil_rejected);
    run("cone_on_unknot_accepted", test_cone_on_unknot_accepted);
    run("cone_on_unlink2_resolvable", test_cone_on_unlink2_resolvable);
    run("cone_on_unlink3_resolvable", test_cone_on_unlink3_resolvable);
    run("cone_on_whitehead_not_resolvable",
        test_cone_on_whitehead_not_resolvable);
    run("cone_on_hopf_still_pruned", test_cone_on_hopf_still_pruned);
    run("boundary_vertex_self_intersection_not_resolvable",
        test_boundary_vertex_self_intersection_not_resolvable);
    run("boundary_vertex_closed_unlinked_petals_not_resolvable",
        test_boundary_vertex_closed_unlinked_petals_not_resolvable);
    run("boundary_filter_rejects_knotted_boundary_petal",
        test_boundary_filter_rejects_knotted_boundary_petal);
    run("hereditariness_stress", test_hereditariness_stress);

    std::cout << bold << "\n=== Summary: " << passed << " passed, "
              << failed_count << " failed ===" << resetColor << "\n";
    return failed_count > 0 ? 1 : 0;
}

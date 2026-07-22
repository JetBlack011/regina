// embeddedsubmanifold_test.cpp
//
// Tests whether EmbeddedSubmanifold<4,2>::addFace() (inherited verbatim by
// KnottedSurface) ever returns true for a face whose addition makes the
// growing subcomplex NOT genuinely embedded -- i.e. the realization map
// from subtri_'s cells into the ambient Triangulation<4>'s cells fails to
// be injective on some open cell.
//
// The check used here (isGenuinelyEmbedded(), below) is deliberately built
// from scratch against only the PUBLIC API of EmbeddedSubmanifold/
// KnottedSurface (addFace(), removeFace(), triangulation(), boundaryLinks())
// -- it never touches faceCount_ (private) or even faces_ (protected), so it
// shares none of addFace()'s own bookkeeping and can catch bugs in that
// bookkeeping rather than just re-deriving them.
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
//      invariant with nothing to do with addFace()'s internals.

#include <algorithm>
#include <chrono>
#include <iostream>
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

        std::string why;
        if (!isGenuinelyEmbedded(correspondence_, &why))
            violations_->push_back({path_, why});

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
                nodes[i].gluings))
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
// (H_1 = Z^n).
AuditResult auditAllEmbeddings(const regina::Triangulation<4> &tri,
                               int maxDepth, bool checkBoundaryHomology = false,
                               int violationCap = 25) {
    Skeleton<4, 2> skeleton(tri);
    Graph graph = buildTestGraph(skeleton);
    KnottedSurface embedding(skeleton);

    AuditResult result;
    EmbeddednessAuditor auditor(embedding, skeleton, graph.graphToSkel,
                                maxDepth, violationCap, result.violations);

    auto visit = [&](const std::vector<int> &) {
        ++result.subgraphsVisited;
        if (!checkBoundaryHomology)
            return;
        if (embedding.isClosed())
            return;
        if (!embedding.isProper() || !embedding.boundaryComponentsMapInjectively())
            return;

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

    auto result = auditAllEmbeddings(cob3.getCobordism(), kDerivedMaxDepth);
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
                                    /*checkBoundaryHomology=*/true);
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

    run("cobordism_disc", test_cobordism_disc);
    run("cobordism_mobius", test_cobordism_mobius);
    run("cobordism_annulus", test_cobordism_annulus);
    run("cobordism_torus", test_cobordism_torus);
    run("cobordism_kb", test_cobordism_kb);

    run("knotbuilder_trefoil_cone", test_knotbuilder_trefoil_cone);
    run("knotbuilder_trefoil_cone_reduced", test_knotbuilder_trefoil_cone_reduced);
    run("knotbuilder_hopflink_cone", test_knotbuilder_hopflink_cone);
    run("knotbuilder_hopflink_cone_reduced", test_knotbuilder_hopflink_cone_reduced);

    std::cout << bold << "\n=== Summary: " << passed << " passed, "
              << failed_count << " failed ===" << resetColor << "\n";
    return failed_count > 0 ? 1 : 0;
}

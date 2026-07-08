// surfacefinder_invariants_test.cpp
//
// John was skeptical of the correctness of SurfaceFinder's DFS after the
// hasUnresolvableConflict() rework (see knottedsurfaces.h/surfacefinder.h
// history). This file does two things that no other test file does:
//
//   1. Checks KnottedSurface's internal bookkeeping (emb_/inv_/indices_,
//      numBoundaryFacets_, improperEdges_, the union-find/imageRootCount_/
//      selfIntersections_ self-intersection tracking, candidatesByFacet_)
//      against independently-recomputed ground truth -- not just at the
//      end of a search, but after *every single mutation* of the real,
//      unmodified production DFS (SurfaceFinder::extend_/findSurfaces),
//      via the debug hook wired up in surfacefinder.h. See
//      KnottedSurface::checkInvariants()'s own comment for why its ground
//      truth deliberately avoids re-deriving the same algorithm a second
//      time (e.g. self-intersection is cross-checked against Regina's own,
//      separately-implemented vertex skeleton on surface_, not a second
//      hand-rolled union-find).
//
//   2. Cross-checks final search results against an independent
//      brute-force ground truth (enumerate every triangle subset, check
//      validity via addTriangle(..., checkSelfIntersection=false) +
//      hasSelfIntersection() + isConnected() + the condition-specific
//      predicate) on triangulations *small enough* for that to be
//      tractable -- computed BEFORE running the real search each time
//      (not after), so as not to bias reading the brute-force result
//      through whatever the DFS already found.
//
// Triangulations are chosen to stress exactly the mechanisms the recent
// investigation touched: self-folded triangles (Klein bottle, Mobius band),
// triangle pairs sharing more than one edge (torus, annulus -- the original
// AdjList bug's failure mode), and genuine multi-triangle branch points at
// a single vertex (Weeks manifold, and "baa"'s already-known 12 Mobius
// bands). Dimension 3 is used wherever it's easier to hand-build/reuse than
// dimension 4, per John's suggestion.

#include <iostream>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <triangulation/dim2.h>
#include <triangulation/dim3.h>
#include <triangulation/dim4.h>
#include <triangulation/example2.h>
#include <triangulation/example3.h>
#include <unistd.h>
#include <vector>

#include "cobordismbuilder.h"
#include "surfacefinder.h"

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

namespace {

const char *condName(SurfaceCondition c) {
    switch (c) {
    case SurfaceCondition::all:
        return "--all";
    case SurfaceCondition::boundary:
        return "--boundary";
    case SurfaceCondition::closed:
        return "--closed";
    }
    return "?";
}

// Does `surf` (fully assembled, deferred checking) satisfy `cond`? Mirrors
// SurfaceFinder::checkWin_()'s own predicate exactly -- this one spot is
// unavoidably "the same formula", since there's only one definition of
// what each condition means; the actual independence in this file comes
// from *never trusting the DFS's own accounting* of which triangle sets
// reach this predicate; see bruteForceGroundTruth's enumeration.
template <int dim>
bool satisfiesCondition(const KnottedSurface<dim> &surf,
                        SurfaceCondition cond) {
    switch (cond) {
    case SurfaceCondition::all:
        return true;
    case SurfaceCondition::boundary:
        return surf.isProper();
    case SurfaceCondition::closed:
        return surf.isClosed();
    }
    return false;
}

// Independent brute-force ground truth: every subset of tri's triangles,
// checked via KnottedSurface::addTriangle(..., checkSelfIntersection=false)
// (deferred whole-set check, order-independent -- see the investigation
// this supports) + hasSelfIntersection() + isConnected() +
// satisfiesCondition(). Computed with its own throwaway SurfaceFinder just
// to source real adjacency data (adjacencyOf) -- buildGluingEdges_ itself
// is checked elsewhere (candidatesByFacet_ inside checkInvariants()), so
// reusing it here for convenience isn't circular for what THIS function is
// actually verifying (the DFS's traversal/pruning decisions).
//
// Deliberately called *before* the real SurfaceFinder search runs at every
// call site below, so the brute-force result is never read with the DFS's
// own answer already in mind.
template <int dim>
std::optional<std::set<std::set<int>>>
bruteForceGroundTruth(const regina::Triangulation<dim> &tri,
                      SurfaceCondition cond, int maxTriangles = 20) {
    int n = tri.countTriangles();
    if (n > maxTriangles || n == 0)
        return std::nullopt;

    SurfaceFinder<dim> finder(tri, cond);
    std::vector<regina::Triangle<dim> *> triangles;
    for (auto *t : tri.triangles())
        triangles.push_back(t);

    std::set<std::set<int>> valid;
    unsigned long total = 1UL << n;
    for (unsigned long mask = 1; mask < total; ++mask) {
        std::vector<int> subset;
        for (int i = 0; i < n; ++i)
            if (mask & (1UL << i))
                subset.push_back(i);

        KnottedSurface<dim> surf(&tri);
        bool ok = true;
        for (int i : subset) {
            if (!surf.addTriangle(triangles[i],
                                  finder.adjacencyOf(triangles[i]), false)) {
                ok = false;
                break;
            }
        }
        if (!ok || surf.hasSelfIntersection() || !surf.surface().isConnected())
            continue;
        if (!satisfiesCondition(surf, cond))
            continue;

        valid.insert(std::set<int>(subset.begin(), subset.end()));
    }
    return valid;
}

// Runs the real SurfaceFinder search with checkInvariants() wired in after
// every single mutation, then (if ground truth was computed) diffs the
// final result exactly against it. `groundTruth` must already have been
// computed by the caller *before* calling this, so the brute-force
// enumeration is never influenced by having already seen the DFS's answer.
template <int dim>
void checkSearch(const std::string &name, const regina::Triangulation<dim> &tri,
                 SurfaceCondition cond,
                 const std::optional<std::set<std::set<int>>> &groundTruth) {
    std::string label = name + " " + condName(cond);
    std::cout << "\n--- " << label << " ---\n";

    SurfaceFinder<dim> finder(tri, cond);

    long hookCalls = 0, violations = 0;
    finder.setDebugHook([&](const char *stage, const KnottedSurface<dim> &s) {
        ++hookCalls;
        std::ostringstream diag;
        if (!s.checkInvariants(diag)) {
            ++violations;
            if (violations <= 5) {
                std::cout << red << "  INVARIANT VIOLATION after stage '"
                          << stage << "' (hook call " << hookCalls << "):\n"
                          << diag.str() << resetColor;
            }
        }
    });

    auto &surfaces = finder.findSurfaces();

    EXPECT_EQ(violations, 0L,
              label + ": zero checkInvariants() violations across " +
                  std::to_string(hookCalls) +
                  " mutations of the live production DFS");

    if (!groundTruth) {
        std::cout << "  (skipped brute-force cross-check: too large)\n";
        return;
    }

    std::set<std::set<int>> found;
    for (const auto &s : surfaces) {
        std::set<int> idx;
        for (auto *t : s.surface().triangles())
            idx.insert(s.image(t)->index());
        found.insert(idx);
    }

    int missing = 0, bogus = 0;
    for (const auto &s : *groundTruth)
        if (!found.count(s))
            ++missing;
    for (const auto &s : found)
        if (!groundTruth->count(s))
            ++bogus;

    std::cout << "  brute-force ground truth: " << groundTruth->size()
              << " surfaces; DFS found: " << found.size() << " (" << missing
              << " missing, " << bogus << " bogus)\n";

    EXPECT_EQ(found.size(), groundTruth->size(),
              label + ": DFS found the same number of surfaces as brute-force "
                      "ground truth");
    EXPECT_EQ(missing, 0, label + ": zero surfaces missing vs. ground truth");
    EXPECT_EQ(bogus, 0, label + ": zero bogus surfaces vs. ground truth");
}

// Computes brute-force ground truth for all 3 conditions *before* running
// any real search (see checkSearch's comment), then runs+checks each.
template <int dim>
void checkTriangulation(const std::string &name,
                        const regina::Triangulation<dim> &tri,
                        int maxBruteForceTriangles = 20) {
    std::cout << "\n=== " << name << " (" << tri.countTriangles()
              << " triangles) ===\n";

    for (SurfaceCondition cond :
         {SurfaceCondition::all, SurfaceCondition::boundary,
          SurfaceCondition::closed}) {
        auto groundTruth =
            bruteForceGroundTruth(tri, cond, maxBruteForceTriangles);
        checkSearch(name, tri, cond, groundTruth);
    }
}

regina::Triangulation<3> thickened(const regina::Triangulation<2> &surf) {
    CobordismBuilder<2> cob(surf);
    return cob.thicken();
}

} // namespace

int main() {
    // isoSig "baa": already independently confirmed (separate session) to
    // have exactly 282 --all surfaces including 12 Mobius bands from
    // genuine branch points that the pre-fix DFS silently dropped. Small
    // enough (10 triangles) for full brute-force + full invariant coverage.
    checkTriangulation("baa (dim 4)", regina::Triangulation<4>("baa"));

    // Weeks manifold: closed, 1 vertex total -- every one of its 18
    // triangles' all 3 corners map to that single ambient vertex, an
    // extreme "everything touches one vertex" stress case, and it
    // independently has a genuine degree-3 branch point (3 pairwise
    // non-edge-adjacent triangles sharing the vertex). Small enough (18
    // triangles) for full brute force.
    checkTriangulation("Weeks manifold (dim 3)", regina::Example<3>::weeks());

    // Self-folded triangles (a triangle with 2 of its own edges identified
    // to the same ambient edge) -- Klein bottle and Mobius band, thickened
    // one dimension up. Originally the first known trigger of the AdjList
    // multi-edge bug.
    checkTriangulation("thickened Klein bottle (dim 3)",
                       thickened(regina::Example<2>::kb()));
    checkTriangulation("thickened Mobius band (dim 3)",
                       thickened(regina::Example<2>::mobius()));

    // Triangle pairs sharing *more than one* edge with each other -- the
    // torus and annulus's exact original AdjList bug repro.
    checkTriangulation("thickened torus (dim 3)",
                       thickened(regina::Example<2>::torus()));
    checkTriangulation("thickened annulus (dim 3)",
                       thickened(regina::Example<2>::annulus()));

    // Smaller thickened cases, cheap to include.
    checkTriangulation("thickened disc (dim 3)",
                       thickened(regina::Example<2>::disc()));

    // thickened sphere (28 triangles) deliberately omitted here: an
    // *unconstrained* findSurfaces() over it (every triangle as a
    // potential seed, no fixed base) runs into the already-documented
    // DFS-call blowup on richer topologies (see the frontier-dedup
    // tradeoff writeup) -- confirmed to grow into multi-gigabyte memory
    // use within seconds when tried here, unrelated to anything this file
    // is investigating. Not a new finding, just the wrong tool for this
    // file's job (checking bookkeeping correctness, not re-litigating a
    // known, separate performance issue).

    std::cout << "\n"
              << bold << (failed_count > 0 ? red : green) << "=== " << passed
              << " passed, " << failed_count << " failed ===" << resetColor
              << "\n";
    return failed_count > 0 ? 1 : 0;
}

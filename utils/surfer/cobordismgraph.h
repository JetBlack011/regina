//
//  cobordismgraph.h
//
//  Created by John Teague on 08/02/2026.
//

#ifndef COBORDISMGRAPH_H

#define COBORDISMGRAPH_H

#include <functional>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

#include "surfacesearch.h"

/*! \file utils/surfer/cobordismgraph.h
 *  \brief The name/genus resolution graph verifyslicegenus.cpp builds up
 *  across its --input rows: given a set of witnessed cobordisms between
 *  named knots/links (each an edge at a specific genus), works out which
 *  rows' slice genus that pins down exactly.
 *
 *  Extracted out of verifyslicegenus.cpp specifically so this logic (previously
 *  untested except via full CLI runs) is unit-testable on its own -- see
 *  tests/cobordismgraph_test.cpp.
 */

namespace cobordismgraph {

/** One --input row: a name to verify/bound the slice genus of, plus its
 * literature bounds and the PD code to search from. */
struct InputRow {
  std::string name;
  std::string pdNotation;
  int lo = 0, hi = 0; // literature genus bounds; lo == hi except for a
                      // handful of [lo;hi] range rows
  int crossings = 0;
};

/** One row of --output: a name's current resolution status and (if
 * resolved) how it was witnessed. */
struct OutputRow {
  std::string knot;
  int resolvedGenus = 0;
  std::string status;       // resolved | range | unresolved | skipped
  std::string witnessKind;  // direct | cobordism | incidental-direct |
                            // incidental-cobordism | propagated | none
  std::string witnessPairSig;
  std::string viaKnot;
  int viaEdgeGenus = 0;
  std::string dependsOn;
  int literatureLo = 0;
  int literatureHi = 0;
};

struct GraphEdge {
  std::string other;
  int genus;
  std::string pairSig;
};

/** name -> every witnessed cobordism edge touching it, both directions. */
using Graph = std::unordered_map<std::string, std::vector<GraphEdge>>;

void addEdge(Graph &graph, const std::string &a, const std::string &b,
            int genus, const std::string &pairSig);

/**
 * Whether `graph` already has an (a, b) edge at exactly `genus` -- see
 * recordWitness()'s use of this for why exact-genus (not "any genus <=
 * this one") is the right dedup key: resolvable() checks an exact
 * equality (target == g + h or target == h - g), not an inequality, so
 * two edges to the same neighbor with different genus are NOT
 * interchangeable, but two with the SAME genus (however many redundant
 * surfaces happened to witness it) are.
 */
bool hasEdgeWithGenus(const Graph &graph, const std::string &a,
                      const std::string &b, int genus);

struct ResolveInfo {
  std::string viaKnot;
  int viaEdgeGenus;
  std::string pairSig; // the resolving edge's own pairSig, if any
};

/**
 * Returns the edge (if any) that pins `name`'s genus to exactly `target`,
 * given `knownGenus`'s currently-established values: an edge (name, other,
 * g) with other's genus known as h pins target only when target == g + h
 * (constructive direction) or target == h - g (contrapositive direction) --
 * merely falling inside [h-g, h+g] is consistent, not conclusive.
 */
std::optional<ResolveInfo>
resolvable(const std::string &name, int target, const Graph &graph,
          const std::unordered_map<std::string, int> &knownGenus);

/** Walks from `viaKnot` back through outputRows' own via_knot chain to
 * whatever ultimately bottomed out (a direct witness, or the Unknot,
 * which has no outputRows entry of its own), joining the path with ';'. */
std::string buildDependsOn(
    const std::string &viaKnot,
    const std::unordered_map<std::string, OutputRow> &outputRows);

/**
 * Sweeps for new resolutions unlockable purely from the current
 * knownGenus/graph state, repeating until a full pass makes no further
 * progress. Two sources of candidates, in one sweep:
 *
 *  - `pending`: --input rows not yet resolved -- resolving one here means
 *    the main loop's own `if (knownGenus.contains(row.name)) continue;`
 *    skips it without ever running a search for it.
 *  - every other name already in `outputRows` that isn't yet resolved --
 *    since --output is a single unified table shared across every run
 *    ever pointed at it, this is typically a knot/link from a *different*
 *    --input than this run's own. Its target is its own already-recorded
 *    literatureLo/Hi, from whenever it was first processed.
 *
 * Either way, each newly-resolved name is recorded into `outputRows` with
 * witness_kind "propagated" (no search run for it this session).
 *
 * Returns every name newly resolved by this call, in resolution order --
 * deliberately not printed here (this function stays pure/side-effect-free,
 * same as recordWitness(), for its own unit-testability); the caller
 * decides whether/how to announce them.
 */
std::vector<std::string>
propagateGraph(const std::vector<InputRow> &pending, const Graph &graph,
              std::unordered_map<std::string, int> &knownGenus,
              std::unordered_map<std::string, OutputRow> &outputRows);

/** The result of recordWitness(): whether it resolved `subject` outright,
 * and whether it detected a mathematical impossibility (a witness
 * implying a genus below `subject`'s own literature lower bound) that
 * should halt the whole program -- see verifyslicegenus.cpp's own
 * flagFatalBug()/haltIfFatalBugDetected(), which the caller is
 * responsible for invoking; this struct only reports the fact and
 * message, it has no side effects of its own on program control flow. */
struct WitnessOutcome {
  bool resolved = false;
  bool fatalBug = false;
  std::string fatalBugMessage;
};

/**
 * Attempts to resolve `subject`'s genus (using its own literature bounds
 * subjectLo/subjectHi and candidate `target`, normally subjectHi) given a
 * newly-found witness of genus `witnessGenus` -- either a direct witness
 * (viaName empty: subject's own boundary alone) or a cobordism to
 * `viaName`. Handles the fatal-bug check, the resolves-now check, edge
 * recording (skipped for a direct witness), and improvement-upper-bound
 * logging to stdout.
 *
 * `capturePairSig` is only ever invoked lazily, at most once, and only
 * when a pairSig is actually about to be used (writing a resolved
 * OutputRow, or adding a genuinely new graph edge) -- never
 * unconditionally, since capturing a pairSig is expensive (see
 * SurfaceFoundInfo::capturePairSig's own doc comment in surfacesearch.h)
 * and the overwhelming majority of witnesses seen in practice are
 * duplicates of an edge already recorded at the same genus.
 *
 * Pure with respect to program control flow: never stops a search or
 * halts the program itself (unlike the search-side branches this
 * generalizes, which additionally call
 * e.requestStop()/e.skipRemainingBoundaryProcessing() on the caller's
 * SurfaceSearch when `subject` is the row currently being searched) --
 * the caller decides what to do with a resolved/fatalBug WitnessOutcome.
 */
WitnessOutcome recordWitness(
    const std::string &subject, int subjectLo, int subjectHi, int target,
    const std::string &viaName, int witnessGenus,
    const std::function<std::string()> &capturePairSig,
    const char *witnessKind, Graph &graph,
    std::unordered_map<std::string, int> &knownGenus,
    std::unordered_map<std::string, OutputRow> &outputRows);

/**
 * One non-search-side ambient boundary component's identity, as classified by
 * splitBoundary(): `name` is always populated (describeBoundary_()
 * guarantees either a single curve name or, for a multi-curve component,
 * a linkName -- never neither), but `safe` says whether that name is
 * trustworthy for a genus deduction. A single curve (one component) is
 * always safe -- no orientation ambiguity is possible with only one
 * component. A multi-curve linkName is only safe when
 * identify::isOrientationSafeName() says so ("Unknot" or an
 * "<n>-component unlink" -- split components don't interact, so
 * orientation doesn't matter); a genuinely linked multi-component name is
 * NOT safe, since different orientations of the same link can have very
 * different true slice genus while sharing one complement (see
 * isOrientationSafeName()'s own doc comment) -- using it for a deduction
 * would silently conflate them.
 */
struct BoundarySide {
  std::string name;
  bool safe;
};

/**
 * Splits a SurfaceBoundaryInfo::boundaryComponents grouping into "the
 * search side" (this row's own side) and every other ambient boundary
 * component the surface touches, each classified via BoundarySide above.
 *
 * The ambient component == searchSideBC is only accepted as the search
 * side if its identified name (the single curve name, or -- for more than
 * one curve -- the linkName) also equals `rowOwnName`. Geometric position
 * (component == searchSideBC) alone is NOT sufficient: the DFS is free to
 * add faces beyond the seeded collar that also touch that same ambient
 * boundary component, so a surface can land the right *count* of curves
 * there (matching this row's own component count) while actually tracing
 * out a completely different, unrelated link -- silently misattributing
 * that other link's genus to this row. `rowOwnName` (identify() run once
 * on this row's own link, before searching) is the only way to catch
 * that: if the names don't match, component == searchSideBC is treated as
 * just another "other" side, exactly like any non-searchSideBC component.
 */
struct BoundarySplit {
  size_t searchCurveCount = 0; // 0 if the search side has no boundary here
  std::vector<BoundarySide> otherSides;
};

BoundarySplit splitBoundary(
    const std::vector<BoundaryComponentNames> &boundaryComponents,
    size_t searchSideBC, const std::string &rowOwnName);

/**
 * A row's own PD-tagged diagram edges (knotbuilder::TriangulationWithLink's
 * `edges`/`reversed`), translated into directed pairs of *vertex indices*
 * within some triangulation combinatorially isomorphic to the row's own
 * knotbuilder triangulation -- normally the ambient search-side boundary
 * component, rebuilt via `BoundaryComponent<4>::build()` (see
 * buildRowOrientation()).
 *
 * Keyed by vertex *index* rather than raw `Vertex<3>*`, since a found
 * surface's own KnottedSurface instance builds its own independent copy of
 * that boundary triangulation (`build()` is called once per KnottedSurface,
 * not shared) -- confirmed empirically that separate `build()` calls on the
 * same ambient boundary component are index-stable (same tetrahedron/vertex
 * numbering every time), so an index survives across objects even though a
 * pointer would not.
 */
struct RowOrientation {
  std::unordered_map<size_t, size_t> headOf; // tail vertex index -> head vertex index
};

/**
 * Builds `rowEdges`/`rowReversed`'s RowOrientation against `searchSideTri`.
 *
 * `searchSideTri` need not be, and normally is not, the same C++ object
 * `rowEdges` themselves live in (typically `rowEdges` live in knotbuilder's
 * own output triangulation, while `searchSideTri` is
 * `ambientTri.boundaryComponent(searchSideBC)->build()`) -- they only need
 * to be *combinatorially isomorphic*, which is established once here via
 * `Triangulation<3>::isIsomorphicTo()` (confirmed empirically this session
 * to correctly identify the row's own tagged edges' corresponding vertices
 * in the rebuilt boundary triangulation -- raw edge/vertex *index*
 * correspondence between the two does NOT hold in general, only up to this
 * isomorphism). This computation is the expensive part of this feature
 * (an isomorphism search), so it is deliberately factored out to run once
 * per row rather than once per found surface -- matchesRowOrientation()
 * below is then cheap.
 *
 * \throws regina::InvalidArgument if `rowEdges` is empty, or if its own
 * triangulation is not isomorphic to `searchSideTri` (should not happen
 * when `searchSideTri` is genuinely built from the ambient boundary
 * component the row's diagram was seeded into).
 */
RowOrientation buildRowOrientation(
    const std::vector<const regina::Edge<3> *> &rowEdges,
    const std::vector<bool> &rowReversed,
    const regina::Triangulation<3> &searchSideTri);

/**
 * Whether `curves` (one found surface's own induced boundary curves on the
 * row's search-side ambient boundary component, from
 * KnottedSurface::orientedBoundaryLinks(), for one arbitrary choice of the
 * surface's two orientations) matches `row` -- either everywhere or
 * nowhere, across every curve (component). A mixed pattern (some curves
 * match, some don't) means the surface's own relative orientation between
 * components doesn't match this row's PD convention, and is rejected here
 * exactly as a fully-mismatched pattern would be reversed by simply
 * picking the surface's other orientation -- see this feature's design
 * notes for the full derivation of why "all or nothing" is the correct
 * acceptance criterion.
 *
 * Also rejects (returns false) if any curve edge isn't one of `row`'s own
 * tagged edges at all -- should not arise when the search's own protected-
 * boundary-component mechanism (see EmbeddingSearch) has kept the search
 * side's edge set fixed to exactly the row's own diagram, but checked
 * defensively rather than assumed.
 */
bool matchesRowOrientation(const RowOrientation &row,
                           const std::vector<OrientedCurve> &curves);

} // namespace cobordismgraph

#endif // COBORDISMGRAPH_H

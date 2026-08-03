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
 * halts the program itself (unlike the "mine" branches this generalizes,
 * which additionally call e.requestStop()/e.skipRemainingBoundaryProcessing()
 * on the caller's SurfaceSearch when `subject` is the row currently being
 * searched) -- the caller decides what to do with a resolved/fatalBug
 * WitnessOutcome.
 */
WitnessOutcome recordWitness(
    const std::string &subject, int subjectLo, int subjectHi, int target,
    const std::string &viaName, int witnessGenus,
    const std::function<std::string()> &capturePairSig,
    const char *witnessKind, Graph &graph,
    std::unordered_map<std::string, int> &knownGenus,
    std::unordered_map<std::string, OutputRow> &outputRows);

/**
 * One non-"mine" ambient boundary component's identity, as classified by
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
 * Splits a SurfaceBoundaryInfo::boundaryComponents grouping into "this
 * row's own side" (the ambient component == knotSideBC, identified
 * geometrically -- never by name, so it needs no safety check of its own)
 * and every other ambient boundary component the surface touches, each
 * classified via BoundarySide above.
 */
struct BoundarySplit {
  size_t mineCurveCount = 0; // 0 if this row's own side has no boundary here
  std::vector<BoundarySide> otherSides;
};

BoundarySplit splitBoundary(
    const std::vector<BoundaryComponentNames> &boundaryComponents,
    size_t knotSideBC);

} // namespace cobordismgraph

#endif // COBORDISMGRAPH_H

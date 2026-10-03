//
//  preconditions.h
//
//  What a found surface must satisfy to witness anything about its row.
//

#ifndef SURFER_COBOUND_PRECONDITIONS_H
#define SURFER_COBOUND_PRECONDITIONS_H

#include <atomic>
#include <map>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include <triangulation/dim3.h>

#include "cobound/cobordisms/cobordism.h"
#include "cobound/outgoing/fromdatabase.h"
#include "cobound/outgoing/outgoinglink.h"
#include "cobound/search/incoming.h"
#include "surfer/enumeration/surfacesearch.h"

/*! \file utils/surfer/cobound/search/preconditions.h
 *  \brief Per found surface: its boundary split into the incoming side and
 *  the others, by geometry, and its incoming curves' orientation against the
 *  row's.
 */

namespace cobordismgraph {

/* Boundary classification */

/**
 * One non-search-side ambient boundary component's identity, as classified
 * by splitBoundary().
 */
struct BoundarySide {
    std::string name;
    int components = 1;
};

/**
 * Splits a SurfaceBoundaryInfo::boundaryComponents grouping into "the
 * search side" (this row's own side) and every other ambient boundary
 * component the surface touches.
 */
struct BoundarySplit {
    size_t searchCurveCount = 0; // 0 if the search side has no boundary here
    bool unnamedSide = false;
    /**< Set when a non-search-side component carried no name at all, which
         describeBoundary_() never produces; the caller treats it as a bug. */
    std::vector<BoundarySide> otherSides;
};

/**
 * The search side is identified by geometry, never by name.
 *
 * Component `searchSideBC` is the search side. Every search is seeded: the
 * seed is L x {0} and no other triangle with an edge in that boundary
 * component is ever searchable, so its curves are L by construction
 * (asserted once per row). Identified names are deliberately not consulted:
 * they are not canonical (a census hit's "#N" varies from one
 * identification to the next), and comparing them once silently discarded
 * whole rows.
 */
BoundarySplit
splitBoundary(const std::vector<BoundaryComponentNames> &boundaryComponents,
              size_t searchSideBC);

/** How a surface's search-side boundary compares with the row's orientation;
 * see classifyRowOrientation(). */
enum class OrientationVerdict {
    match,        /**< Some choice of orientation on each surface component
                       induces the row's own orientation. */
    mismatch,     /**< Some surface component's curves induce an orientation
                       pattern no flip of that component can fix: the surface
                       witnesses a different oriented variant of the link. */
    foreignEdge,  /**< A search-side edge is not one of the row's own. */
    incoherentCurve, /**< A single curve's edges disagree in direction, or a
                          curve's surface component is unknown. */
};

/** judgeRowOrientation()'s answer: the verdict, and the flips it implies. */
struct RowOrientationJudgement {
    OrientationVerdict verdict = OrientationVerdict::mismatch;
    /** Per surface component met by a curve, +1 when its curves run as the
     *  row's link does, -1 when against: every component's, when the
     *  verdict is match. */
    std::map<size_t, int> flips;
    /** No non-empty curve at all (the verdict is then mismatch). */
    bool noCurves = false;

    /** The flips when every component is consistently oriented: with
     *  verdict match, or (an empty map) with no curves at all; else
     *  nullopt. farside::incomingFlips() is this. */
    std::optional<std::map<size_t, int>> consistentFlips() const {
        if (verdict == OrientationVerdict::match || noCurves)
            return flips;
        return std::nullopt;
    }
};

/**
 * Compares one found surface's induced boundary orientation on the
 * search side with `row`: the one walk behind classifyRowOrientation() and
 * farside::incomingFlips().
 *
 * `curves` come from KnottedSurface::orientedBoundaryLinks(), which orients
 * each CONNECTED COMPONENT of the surface independently and arbitrarily;
 * `surfaceComponentOf` (KnottedSurface::boundaryEdgeSurfaceComponent())
 * says which component each edge belongs to. So the curves are grouped by
 * surface component, and each group must match `row` all at once or be
 * reversed all at once; different groups are independent. That is exactly
 * the freedom the surface has: its components can be oriented separately
 * and still tubed into one oriented surface, because a split two-component
 * unlink bounds an oriented annulus for either relative orientation.
 *
 * `foreignEdge` and `incoherentCurve` cannot happen for a correctly built
 * row in a seeded search; the caller treats them as bugs.
 */
RowOrientationJudgement judgeRowOrientation(
    const RowOrientation &row, const std::vector<OrientedCurve> &curves,
    const std::map<const regina::Edge<3> *, size_t> &surfaceComponentOf);

/** judgeRowOrientation()'s verdict. */
OrientationVerdict classifyRowOrientation(
    const RowOrientation &row, const std::vector<OrientedCurve> &curves,
    const std::map<const regina::Edge<3> *, size_t> &surfaceComponentOf);

} // namespace cobordismgraph

namespace farside {
class DiagramNamer;
}

namespace rowsearch {

/** Why a surface does or does not witness anything about its row. */
enum class Gate {
    accepted,
    nonOrientable,       /**< Impossible: orientableOnly prunes these. */
    unnamedSide,         /**< Impossible: a far side with no name at all. */
    searchSideBroken,    /**< Impossible when seeded. */
    orientation,         /**< Witnesses another oriented variant of the row. */
    orientationBroken,   /**< Impossible: an incoherent or foreign curve. */
    multiFarSide,        /**< Impossible: S^3 x I has two boundary components. */
};

/** The name --rejection-sample-log records a rejection under. */
const char *gateReason(Gate gate);

/** A found surface, judged against its row by gateSurface(). */
struct GatedSurface {
    Gate gate = Gate::accepted;
    cobordismgraph::BoundarySplit split;
    /** Captured only once the search side has passed (so from the
     *  orientation gate on). Per ambient boundary component, its oriented
     *  curves, each surface component oriented independently. */
    std::vector<std::pair<size_t, std::vector<OrientedCurve>>> orientedLinks;
    /** Which surface component each boundary edge lies on. */
    std::map<const regina::Edge<3> *, size_t> surfaceOf;
    std::vector<OrientedCurve> searchSideCurves;
    /** Each surface component's flip against the row (its incoming curves'
     *  judgeRowOrientation()): complete for an accepted surface. */
    std::map<size_t, int> flips;

    bool accepted() const { return gate == Gate::accepted; }
};

/**
 * Judges `info` against `row`, in this order: orientable; the boundary split
 * (search side by geometry, never by name); the search side holding the
 * row's component count; the row's own orientation, per surface component
 * (cobordismgraph::classifyRowOrientation()); at most one far side.
 */
GatedSurface gateSurface(const SurfaceBoundaryInfo &info, const RowBuild &row);

/**
 * The name a witness records for an accepted surface's one far side: the
 * namer's name, normalized (identify() decorates a translated census hit as
 * "4_1 (m004 : #1)" while the tables call it "4_1"). With exact names on
 * (`namer->exactNamesOn()`), a multi-component far side takes its exact
 * oriented name when there is one: two surfaces whose far sides are
 * different orientation variants of one link must be two witnesses.
 *
 * \pre `g` is accepted with exactly one far side.
 */
std::string farSideName(const GatedSurface &g, const farside::DiagramNamer *namer);

/**
 * Every surface the drain describes lands in exactly one of these, and at
 * row end they must add up to what the search accepted. A surface can never
 * vanish between the search and the witness record without being counted.
 */
struct RowAccounting {
    std::atomic<long long> described{0};
    std::atomic<long long> recorded{0};
    std::atomic<long long> duplicate{0};
    std::atomic<long long> orientation{0}; // another oriented variant
    // Impossible for a correct build; any nonzero count halts the run.
    std::atomic<long long> nonOrientable{0};
    std::atomic<long long> searchSideBroken{0};
    std::atomic<long long> orientationBroken{0};
    std::atomic<long long> multiFarSide{0};
    std::atomic<long long> unnamedSide{0};

    /** Counts one surface rejected by `gate` (not Gate::accepted). */
    void reject(Gate gate);

    long long impossible() const {
        return nonOrientable + searchSideBroken + orientationBroken +
               multiFarSide + unnamedSide;
    }
    long long bucketed() const {
        return recorded + duplicate + orientation + impossible();
    }
    /** Surfaces were described, yet none reached the witness record. That
     *  can be genuine (every one witnesses another oriented variant), but it
     *  is also what a broken gate looks like, so it licenses no negative. */
    bool nothingExamined() const {
        return described > 0 && recorded + duplicate == 0;
    }

    /**
     * Why this row's surfaces do not add up, or empty if they do. Every
     * surface the search accepted must have been described by the drain
     * (unless the drain was deliberately cut short), and every described
     * surface must sit in exactly one bucket, none of them impossible.
     */
    std::string failure(long long accepted, long long rebuildFailed,
                        bool drainSkipped) const;

    /** The body of the `accounting:` line tools/orchestrate/dispatch.py
     *  parses (RE_ACCOUNTING): "accepted N, described N, ..., ok". */
    std::string summary(long long accepted, bool drainSkipped) const;
};

/**
 * What a search keeps one cobordism per: the cobordism's identity
 * (cobordismgraph::witnessIdentity(), the one identity every dedupe uses),
 * then -- the grouping a goal run's graph needs -- which row components and
 * how many outgoing curves each surface component carries, as a canonical
 * string (surface components are unlabelled, so the entries are sorted).
 */
std::string keptKey(const cobordismgraph::Witness &w, const farside::OutgoingLink &link,
                    const farside::WitnessRedrawer &row);
} // namespace rowsearch

#endif // SURFER_COBOUND_PRECONDITIONS_H

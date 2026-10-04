//
//  preconditions.h
//
//  What a found surface must satisfy to be a cobordism from the incoming link.
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
 *  incoming link's.
 */

namespace search {

/* Boundary classification */

/**
 * One ambient boundary component other than the incoming one, as classified
 * by splitBoundary().
 */
struct BoundarySide {
    std::string name;
    int components = 1;
};

/**
 * Splits a SurfaceBoundaryInfo::boundaryComponents grouping into "the
 * incoming side" (the incoming link's) and every other ambient boundary
 * component the surface touches.
 */
struct BoundarySplit {
    size_t searchCurveCount = 0; // 0 if the incoming side has no boundary here
    bool unnamedSide = false;
    /**< Set when a component other than the incoming one carried no name at all, which
         describeBoundary_() never produces; the caller treats it as a bug. */
    std::vector<BoundarySide> otherSides;
};

/**
 * The incoming side is told by geometry, never by name.
 *
 * Component `incomingBC` is the incoming side. Every search is seeded: the
 * seed is L x {0} and no other triangle with an edge in that boundary
 * component is ever searchable, so its curves are L by construction
 * (asserted once per search). Names are deliberately not consulted: they
 * are not canonical (a census hit's "#N" varies from one naming to the
 * next), and comparing them once silently discarded whole searches.
 */
BoundarySplit
splitBoundary(const std::vector<BoundaryComponentNames> &boundaryComponents,
              size_t incomingBC);

/** How a surface's incoming boundary compares with the incoming link's orientation;
 * see classifyIncomingOrientation(). */
enum class OrientationVerdict {
    match,        /**< Some choice of orientation on each surface component
                       induces the incoming link's own orientation. */
    mismatch,     /**< Some surface component's curves induce an orientation
                       pattern no flip of that component can fix: the surface
                       is a cobordism from a different oriented variant of the link. */
    foreignEdge,  /**< An incoming edge is not one of the incoming link's. */
    incoherentCurve, /**< A single curve's edges disagree in direction, or a
                          curve's surface component is unknown. */
};

/** judgeIncomingOrientation()'s answer: the verdict, and the flips it implies. */
struct IncomingOrientationJudgement {
    OrientationVerdict verdict = OrientationVerdict::mismatch;
    /** Per surface component met by a curve, +1 when its curves run as the
     *  incoming link does, -1 when against: every component's, when the
     *  verdict is match. */
    std::map<size_t, int> flips;
    /** No non-empty curve at all (the verdict is then mismatch). */
    bool noCurves = false;

    /** The flips when every component is consistently oriented: with
     *  verdict match, or (an empty map) with no curves at all; else
     *  nullopt. outgoing::incomingFlips() is this. */
    std::optional<std::map<size_t, int>> consistentFlips() const {
        if (verdict == OrientationVerdict::match || noCurves)
            return flips;
        return std::nullopt;
    }
};

/**
 * Compares one found surface's induced boundary orientation on the
 * incoming side with `incoming`: the one walk behind classifyIncomingOrientation() and
 * outgoing::incomingFlips().
 *
 * `curves` come from KnottedSurface::orientedBoundaryLinks(), which orients
 * each CONNECTED COMPONENT of the surface independently and arbitrarily;
 * `surfaceComponentOf` (KnottedSurface::boundaryEdgeSurfaceComponent())
 * says which component each edge belongs to. So the curves are grouped by
 * surface component, and each group must match `incoming` all at once or be
 * reversed all at once; different groups are independent. That is exactly
 * the freedom the surface has: its components can be oriented separately
 * and still tubed into one oriented surface, because a split two-component
 * unlink bounds an oriented annulus for either relative orientation.
 *
 * `foreignEdge` and `incoherentCurve` cannot happen for a correctly built
 * incoming link in a seeded search; the caller treats them as bugs.
 */
IncomingOrientationJudgement judgeIncomingOrientation(
    const IncomingOrientation &incoming, const std::vector<OrientedCurve> &curves,
    const std::map<const regina::Edge<3> *, size_t> &surfaceComponentOf);

/** judgeIncomingOrientation()'s verdict. */
OrientationVerdict classifyIncomingOrientation(
    const IncomingOrientation &incoming, const std::vector<OrientedCurve> &curves,
    const std::map<const regina::Edge<3> *, size_t> &surfaceComponentOf);

} // namespace search

namespace outgoing {
class OutgoingNamer;
}

namespace search {

/** Why a surface is, or is not, a cobordism from its incoming link. */
enum class Gate {
    accepted,
    nonOrientable,       /**< Impossible: orientableOnly prunes these. */
    unnamedSide,         /**< Impossible: an outgoing link with no name at all. */
    incomingBroken,    /**< Impossible when seeded. */
    orientation,         /**< From another oriented variant of the incoming link. */
    orientationBroken,   /**< Impossible: an incoherent or foreign curve. */
    multiOutgoing,        /**< Impossible: S^3 x I has two boundary components. */
};

/** The name --rejection-sample-log records a rejection under. */
const char *gateReason(Gate gate);

/** A found surface, judged against its incoming link by gateSurface(). */
struct GatedSurface {
    Gate gate = Gate::accepted;
    search::BoundarySplit split;
    /** Captured only once the incoming side has passed (so from the
     *  orientation gate on). Per ambient boundary component, its oriented
     *  curves, each surface component oriented independently. */
    std::vector<std::pair<size_t, std::vector<OrientedCurve>>> orientedLinks;
    /** Which surface component each boundary edge lies on. */
    std::map<const regina::Edge<3> *, size_t> surfaceOf;
    std::vector<OrientedCurve> incomingCurves;
    /** Each surface component's flip against the incoming link (its incoming curves'
     *  judgeIncomingOrientation()): complete for an accepted surface. */
    std::map<size_t, int> flips;

    bool accepted() const { return gate == Gate::accepted; }
};

/**
 * Judges `info` against `thickened`, in this order: orientable; the boundary split
 * (incoming side by geometry, never by name); the incoming side holding the
 * incoming link's component count; its own orientation, per surface component
 * (search::classifyIncomingOrientation()); at most one outgoing link.
 */
GatedSurface gateSurface(const SurfaceBoundaryInfo &info, const IncomingThickening &thickened);

/**
 * The name a cobordism records for an accepted surface's one outgoing link: the
 * namer's name, normalized (census::nameComplement() decorates a translated census hit as
 * "4_1 (m004 : #1)" while the tables call it "4_1"). With names on
 * (`namer->orientedNamesOn()`), a multi-component outgoing link takes its
 * oriented name when there is one: two surfaces whose outgoing links are
 * different orientation variants of one link must be two cobordisms.
 *
 * \pre `g` is accepted with exactly one outgoing link.
 */
std::string nameOutgoing(const GatedSurface &g, const outgoing::OutgoingNamer *namer);

/**
 * Every surface the drain describes lands in exactly one of these, and at
 * the search's end they must add up to what the search accepted. A surface can never
 * vanish between the search and the cobordism record without being counted.
 */
struct SearchAccounting {
    std::atomic<long long> described{0};
    std::atomic<long long> recorded{0};
    std::atomic<long long> duplicate{0};
    std::atomic<long long> orientation{0}; // another oriented variant
    // Impossible for a correct build; any nonzero count halts the run.
    std::atomic<long long> nonOrientable{0};
    std::atomic<long long> incomingBroken{0};
    std::atomic<long long> orientationBroken{0};
    std::atomic<long long> multiOutgoing{0};
    std::atomic<long long> unnamedSide{0};

    /** Counts one surface rejected by `gate` (not Gate::accepted). */
    void reject(Gate gate);

    long long impossible() const {
        return nonOrientable + incomingBroken + orientationBroken +
               multiOutgoing + unnamedSide;
    }
    long long bucketed() const {
        return recorded + duplicate + orientation + impossible();
    }
    /** Surfaces were described, yet none reached the cobordism record. That
     *  can be genuine (every one is a cobordism from another oriented variant), but it
     *  is also what a broken gate looks like, so it licenses no negative. */
    bool nothingExamined() const {
        return described > 0 && recorded + duplicate == 0;
    }

    /**
     * Why this search's surfaces do not add up, or empty if they do. Every
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
 * (cobordisms::cobordismIdentity(), the one identity every dedupe uses),
 * then -- the grouping a goal run's graph needs -- which incoming components and
 * how many outgoing curves each surface component carries, as a canonical
 * string (surface components are unlabelled, so the entries are sorted).
 */
std::string keptKey(const cobordisms::Cobordism &w, const outgoing::OutgoingLink &link,
                    const outgoing::OutgoingReader &reader);
} // namespace search

#endif // SURFER_COBOUND_PRECONDITIONS_H

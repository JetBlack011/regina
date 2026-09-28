#ifndef ROWSEARCH_H
#define ROWSEARCH_H

// The row pipeline, shared by verifyslicegenus and cascadesearch.
//
// A "row" is one link, searched for surfaces running from it across
// S^3 x I. These are the pieces every driver needs and none may duplicate,
// since each is a place a surface can silently go missing:
//
//   - buildRow(): knotbuilder's T for the PD code, thickened into S^3 x I,
//     the collar L x [0, collarLayers] as the seed, and the row map taking
//     L onto the seed's own edges (with its setup-time checks);
//   - gateSurface(): the checks a found surface must pass to witness
//     anything about this row (orientable, search side intact, the row's own
//     oriented variant, one far side);
//   - RowAccounting: every surface the drain describes lands in exactly one
//     bucket, and the buckets must add up to what the search accepted;
//   - RowWatchdog: ends a row's search at its surface target or deadline;
//   - conditionFor(): which BoundaryCondition a row is searched under.
//
// Drivers keep their own printing, witness records and conclusions.

#include <atomic>
#include <chrono>
#include <functional>
#include <map>
#include <optional>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include <triangulation/dim3.h>
#include <triangulation/dim4.h>

#include "cobordismbuilder.h"
#include "cobordismgraph.h"
#include "embeddedsubmanifold.h"
#include "knotbuilder/knotbuilder.h"
#include "surfacesearch.h"

namespace farside {
class DiagramNamer;
}

namespace rowsearch {

/** Which BoundaryCondition to search a row under; see conditionFor(). */
enum class BoundaryConditionMode { automatic, connected, proper };

/**
 * The BoundaryCondition for a row of `componentCount` components.
 *
 * `connected` (one curve per ambient boundary component) prunes far harder,
 * but it caps the curve count on the far side too, so a knot row searched
 * under it can never find a knot-to-link cobordism. It is honoured only for
 * a knot: a multi-component link can never meet it on its own search side.
 * `automatic` is `connected` for a knot and `proper` otherwise.
 */
BoundaryCondition conditionFor(BoundaryConditionMode mode, int componentCount);

/**
 * A row's ambient S^3 x I and its seed. Filled in place by buildRow() and
 * never moved: a DiagramNamer and a SurfaceSearch built from it hold
 * pointers into `link.tri` and `cob`.
 */
struct RowBuild {
    RowBuild() = default;
    RowBuild(const RowBuild &) = delete;
    RowBuild &operator=(const RowBuild &) = delete;

    knotbuilder::PDCode pdcode;
    knotbuilder::TriangulationWithLink link; /**< T, and L's edges in it. */
    std::optional<CobordismBuilder<3>> cob;
    regina::Triangulation<4> tri;  /**< The search's ambient. */
    std::vector<int> seedFaces;    /**< The collar; empty without one. */
    size_t searchSideBC = 0;       /**< The incoming boundary, T x {0}. */
    int componentCount = 1;        /**< Of the row's link. */

    std::optional<cobordismgraph::RowOrientation> orientation;
    /**< The row map: L's edges and PD orientation in search-side terms. */
    std::vector<size_t> searchEdges;
    /**< The row's own link on the search side, as sorted edge indices of that
         boundary component's built triangulation. Seeded, the seed's own
         edges there (L x {0}); unseeded, the image of L under the row map,
         which splitBoundary() filters on. */
};

/**
 * Builds `row` for PD code `pdNotation`: buildAmbient(), then orientRow().
 *
 * \throws regina::InvalidArgument as either does.
 */
void buildRow(const std::string &pdNotation, int thickenLayers,
              int collarLayers, bool useCone, RowBuild &row);

/**
 * The ambient alone: T, `thickenLayers` thickenings with a collar through
 * the first `collarLayers`, and an optional cone. Leaves `orientation` and
 * `searchEdges` empty.
 *
 * \throws regina::InvalidArgument for an unparseable or unbuildable PD code.
 */
void buildAmbient(const std::string &pdNotation, int thickenLayers,
                  int collarLayers, bool useCone, RowBuild &row);

/**
 * The row map for an ambient built by buildAmbient(). Checks, once, that the
 * search side holds exactly L's edges in L's number of components.
 *
 * \throws regina::InvalidArgument for a row map that cannot be built or
 * fails those checks.
 */
void orientRow(RowBuild &row);

/** Why a surface does or does not witness anything about its row. */
enum class Gate {
    accepted,
    nonOrientable,       /**< Impossible: orientableOnly prunes these. */
    unnamedSide,         /**< Impossible: a far side with no name at all. */
    searchSideElsewhere, /**< Unseeded only: another link on the search side. */
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
std::string farSideName(const GatedSurface &g, const RowBuild &row,
                        const farside::DiagramNamer *namer);

/**
 * Every surface the drain describes lands in exactly one of these, and at
 * row end they must add up to what the search accepted. A surface can never
 * vanish between the search and the witness record without being counted.
 */
struct RowAccounting {
    std::atomic<long long> described{0};
    std::atomic<long long> recorded{0};
    std::atomic<long long> duplicate{0};
    std::atomic<long long> orientation{0};         // another oriented variant
    std::atomic<long long> searchSideElsewhere{0}; // unseeded only
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
        return recorded + duplicate + orientation + searchSideElsewhere +
               impossible();
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

/** When RowWatchdog ends a row's search. Unset limits never fire. */
struct WatchdogLimits {
    std::optional<long long> surfaceTarget;
    std::optional<double> rowSeconds;
    std::optional<double> sweepSeconds;
    std::chrono::steady_clock::time_point sweepStart{};
    std::optional<double> quiescenceSeconds;
    /** Milliseconds since the row last taught us something new; read only
     *  with `quiescenceSeconds`. */
    std::function<long long()> idleMillis;

    bool any() const {
        return surfaceTarget || rowSeconds || sweepSeconds || quiescenceSeconds;
    }
};

/**
 * Polls every 200 ms, while the search and its drain run, and calls
 * `endRow(why)` at most once: "surface-target", "timeout" or "quiescent".
 * The surface target is checked before the clocks, so a row reaching it in
 * the same tick as a deadline records "surface-target": the two mean
 * different things to anyone later reading the negative.
 *
 * No thread is started when no limit is set.
 */
class RowWatchdog {
public:
    RowWatchdog(WatchdogLimits limits, std::function<void(const char *)> endRow);
    ~RowWatchdog();
    RowWatchdog(const RowWatchdog &) = delete;
    RowWatchdog &operator=(const RowWatchdog &) = delete;

    /** The live count of surfaces satisfying the boundary condition, from
     *  SearchCallbacks::onProgress (the only place it is handed out). */
    void publishSatisfying(long long count) {
        satisfying_.store(count, std::memory_order_relaxed);
    }
    /** Stops polling and joins; idempotent. Call once the search returns. */
    void stop();

private:
    WatchdogLimits limits_;
    std::function<void(const char *)> endRow_;
    std::atomic<long long> satisfying_{0};
    std::atomic<bool> done_{false};
    std::thread thread_;
};

} // namespace rowsearch

#endif // ROWSEARCH_H

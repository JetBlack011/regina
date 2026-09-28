#include "rowsearch.h"

#include "collar.h"
#include "farsidecurves.h"
#include "farsidenaming.h"
#include "linkcomplement.h"

namespace rowsearch {

BoundaryCondition conditionFor(BoundaryConditionMode mode, int componentCount) {
    switch (mode) {
    case BoundaryConditionMode::proper:
        return BoundaryCondition::proper;
    case BoundaryConditionMode::connected:
    case BoundaryConditionMode::automatic:
    default:
        return componentCount == 1 ? BoundaryCondition::connected
                                   : BoundaryCondition::proper;
    }
}

void buildAmbient(const std::string &pdNotation, int thickenLayers,
                  int collarLayers, bool useCone, RowBuild &row) {
    row.pdcode = knotbuilder::parsePDCode(pdNotation);
    row.link = knotbuilder::buildLink(row.pdcode);

    auto &[t2, edges2, reversed2] = row.link;
    Link linkGrouping(t2, edges2);
    row.componentCount = linkGrouping.countComponents();

    std::vector<int> edgeIndices;
    edgeIndices.reserve(edges2.size());
    for (const regina::Edge<3> *e : edges2)
        edgeIndices.push_back(static_cast<int>(e->index()));

    // Indices are preserved across CobordismBuilder's internal copy of t2
    // (see CobordismBuilder::baseTriangulation()), so edgeIndices still name
    // L's edges in the cobordism's base. The collar must be extended on
    // every layer it is meant to cover: CollarBuilder::addLayer() captures
    // only the most recently built layer's prisms.
    row.cob.emplace(t2);
    CobordismBuilder<3> &cob = *row.cob;
    CollarBuilder collarBuilder(edgeIndices);
    for (int i = 0; i < thickenLayers; ++i) {
        cob.thicken();
        if (i < collarLayers)
            collarBuilder.addLayer(cob);
    }
    if (useCone)
        cob.cone();

    row.searchSideBC = cob.baseBoundaryComponent()->index();
    row.tri = cob.getCobordism();

    if (collarLayers > 0)
        for (regina::Triangle<4> *t : collarBuilder.resolve())
            row.seedFaces.push_back(static_cast<int>(t->index()));
}

void orientRow(RowBuild &row) {
    const auto &edges2 = row.link.edges;
    const auto &reversed2 = row.link.reversed;
    // The seed's own edges on the search side: exactly L x {0}.
    if (!row.seedFaces.empty())
        row.searchEdges = farside::boundaryEdgesOf(row.tri, row.seedFaces,
                                                   row.searchSideBC);
    row.orientation = cobordismgraph::buildRowOrientation(
        edges2, reversed2, row.tri.boundaryComponent(row.searchSideBC)->build(),
        row.seedFaces.empty() ? nullptr : &row.searchEdges);
    if (row.seedFaces.empty())
        row.searchEdges = row.orientation->edges;

    // Setup-time checks on the row's own link, in place of any per-surface
    // ones: the search side is fixed from here on.
    if (row.searchEdges.size() != edges2.size())
        throw regina::InvalidArgument(
            "the search side holds " + std::to_string(row.searchEdges.size()) +
            " link edges, the diagram " + std::to_string(edges2.size()));
    if (row.orientation->components !=
        static_cast<size_t>(row.componentCount))
        throw regina::InvalidArgument(
            "the search-side link has " +
            std::to_string(row.orientation->components) +
            " components, the diagram " + std::to_string(row.componentCount));
}

void buildRow(const std::string &pdNotation, int thickenLayers,
              int collarLayers, bool useCone, RowBuild &row) {
    buildAmbient(pdNotation, thickenLayers, collarLayers, useCone, row);
    orientRow(row);
}

const char *gateReason(Gate gate) {
    switch (gate) {
    case Gate::accepted: return "accepted";
    case Gate::nonOrientable: return "non-orientable";
    case Gate::unnamedSide: return "unnamed-side";
    case Gate::searchSideElsewhere: return "search-side-elsewhere";
    case Gate::searchSideBroken: return "search-side-broken";
    case Gate::orientation: return "orientation";
    case Gate::orientationBroken: return "orientation-broken";
    case Gate::multiFarSide: return "multi-far-side";
    }
    return "unknown";
}

GatedSurface gateSurface(const SurfaceBoundaryInfo &info, const RowBuild &row) {
    GatedSurface g;
    if (!info.orientable) {
        // orientableOnly=true prunes these during the search.
        g.gate = Gate::nonOrientable;
        return g;
    }

    // Seeded, the search side is L by construction (asserted once at row
    // setup), so splitBoundary() just takes component searchSideBC.
    // Unseeded, it filters on L's own edges, setwise.
    g.split = cobordismgraph::splitBoundary(
        info.boundaryComponents, row.searchSideBC,
        row.seedFaces.empty() ? &row.searchEdges : nullptr);
    if (g.split.unnamedSide) {
        g.gate = Gate::unnamedSide;
        return g;
    }
    if (g.split.searchCurveCount != static_cast<size_t>(row.componentCount)) {
        // Unseeded: a different link on the search side, so whatever this
        // surface witnesses, it isn't about this row.
        g.gate = row.seedFaces.empty() ? Gate::searchSideElsewhere
                                       : Gate::searchSideBroken;
        return g;
    }

    // Orientation: a surface component whose curves induce a pattern no flip
    // of that component can fix witnesses a DIFFERENT oriented variant of
    // this link -- the L6a3{0}/L6a3{1} misattribution. Judged per surface
    // component, since each can be oriented independently
    // (classifyRowOrientation()).
    g.orientedLinks = info.captureOrientedBoundaryLinks();
    g.surfaceOf = info.captureBoundaryEdgeSurfaceComponent();
    bool foundSearchSide = false;
    for (auto &[c, curves] : g.orientedLinks) {
        if (c == row.searchSideBC) {
            g.searchSideCurves = curves;
            foundSearchSide = true;
            break;
        }
    }
    const cobordismgraph::OrientationVerdict verdict =
        foundSearchSide ? cobordismgraph::classifyRowOrientation(
                              *row.orientation, g.searchSideCurves, g.surfaceOf)
                        : cobordismgraph::OrientationVerdict::incoherentCurve;
    if (verdict == cobordismgraph::OrientationVerdict::mismatch) {
        g.gate = Gate::orientation;
        return g;
    }
    if (verdict != cobordismgraph::OrientationVerdict::match) {
        g.gate = Gate::orientationBroken;
        return g;
    }

    if (g.split.otherSides.size() > 1)
        g.gate = Gate::multiFarSide;
    return g;
}

std::string farSideName(const GatedSurface &g, const RowBuild &row,
                        const farside::DiagramNamer *namer) {
    const cobordismgraph::BoundarySide &far = g.split.otherSides.front();
    std::string name = cobordismgraph::normalizeIdentifiedName(far.name);
    if (far.components > 1 && namer && namer->exactNamesOn()) {
        if (auto flips = farside::incomingFlips(
                *row.orientation, g.searchSideCurves, g.surfaceOf)) {
            for (const auto &[bc, curves] : g.orientedLinks)
                if (namer->handles(bc))
                    if (auto n = namer->orientedName(curves, g.surfaceOf, *flips))
                        name = *n;
        }
    }
    return name;
}

void RowAccounting::reject(Gate gate) {
    switch (gate) {
    case Gate::accepted: break;
    case Gate::nonOrientable:
        nonOrientable.fetch_add(1, std::memory_order_relaxed);
        break;
    case Gate::unnamedSide:
        unnamedSide.fetch_add(1, std::memory_order_relaxed);
        break;
    case Gate::searchSideElsewhere:
        searchSideElsewhere.fetch_add(1, std::memory_order_relaxed);
        break;
    case Gate::searchSideBroken:
        searchSideBroken.fetch_add(1, std::memory_order_relaxed);
        break;
    case Gate::orientation:
        orientation.fetch_add(1, std::memory_order_relaxed);
        break;
    case Gate::orientationBroken:
        orientationBroken.fetch_add(1, std::memory_order_relaxed);
        break;
    case Gate::multiFarSide:
        multiFarSide.fetch_add(1, std::memory_order_relaxed);
        break;
    }
}

std::string RowAccounting::failure(long long accepted, long long rebuildFailed,
                                   bool drainSkipped) const {
    const long long nDescribed = described.load();
    if (bucketed() != nDescribed)
        return std::to_string(nDescribed) + " surfaces described but " +
               std::to_string(bucketed()) + " accounted for";
    if (rebuildFailed > 0)
        return std::to_string(rebuildFailed) +
               " accepted surfaces failed to rebuild in the drain";
    if (!drainSkipped && nDescribed != accepted)
        return std::to_string(accepted) + " surfaces accepted but " +
               std::to_string(nDescribed) + " described";
    if (impossible() > 0)
        return std::to_string(impossible()) +
               " surfaces hit a state that cannot occur (non-orientable " +
               std::to_string(nonOrientable.load()) + ", search side " +
               std::to_string(searchSideBroken.load()) + ", orientation " +
               std::to_string(orientationBroken.load()) + ", multi far side " +
               std::to_string(multiFarSide.load()) + ", unnamed side " +
               std::to_string(unnamedSide.load()) + ")";
    return {};
}

std::string RowAccounting::summary(long long accepted, bool drainSkipped) const {
    return "accepted " + std::to_string(accepted) + ", described " +
           std::to_string(described.load()) + ", recorded " +
           std::to_string(recorded.load()) + ", duplicate " +
           std::to_string(duplicate.load()) + ", other-orientation " +
           std::to_string(orientation.load()) + ", search-side-elsewhere " +
           std::to_string(searchSideElsewhere.load()) + ", impossible " +
           std::to_string(impossible()) + ", drain " +
           (drainSkipped ? "skipped" : "complete") + ", " +
           (nothingExamined() ? "WARNING" : "ok");
}

RowWatchdog::RowWatchdog(WatchdogLimits limits,
                         std::function<void(const char *)> endRow)
    : limits_(std::move(limits)), endRow_(std::move(endRow)) {
    if (!limits_.any())
        return;
    thread_ = std::thread([this] {
        const auto rowDeadline =
            std::chrono::steady_clock::now() +
            std::chrono::duration<double>(limits_.rowSeconds.value_or(0));
        while (!done_.load(std::memory_order_relaxed)) {
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            if (done_.load(std::memory_order_relaxed))
                break;
            // Checked before the clocks: see the class comment.
            if (limits_.surfaceTarget &&
                satisfying_.load(std::memory_order_relaxed) >=
                    *limits_.surfaceTarget) {
                endRow_("surface-target");
                break;
            }
            const auto now = std::chrono::steady_clock::now();
            if (limits_.rowSeconds && now >= rowDeadline) {
                endRow_("timeout");
                break;
            }
            if (limits_.sweepSeconds &&
                now - limits_.sweepStart >
                    std::chrono::duration<double>(*limits_.sweepSeconds)) {
                endRow_("timeout");
                break;
            }
            // Quiescence: this row has stopped teaching us anything new.
            // Meaningful during the drain too, which is where witnesses are
            // actually identified.
            if (limits_.quiescenceSeconds &&
                limits_.idleMillis() >
                    static_cast<long long>(*limits_.quiescenceSeconds * 1000)) {
                endRow_("quiescent");
                break;
            }
        }
    });
}

void RowWatchdog::stop() {
    done_.store(true, std::memory_order_relaxed);
    if (thread_.joinable())
        thread_.join();
}

RowWatchdog::~RowWatchdog() { stop(); }

} // namespace rowsearch

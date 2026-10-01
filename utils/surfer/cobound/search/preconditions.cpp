//
//  preconditions.cpp
//

#include "cobound/search/preconditions.h"

#include <algorithm>
#include <optional>
#include <unordered_map>

#include "cobound/outgoing/outgoingnamer.h"
#include "cobound/outgoing/outgoinglink.h"
#include "linknaming/names.h"

namespace cobordismgraph {

/* Boundary classification */

BoundarySplit
splitBoundary(const std::vector<BoundaryComponentNames> &boundaryComponents,
              size_t searchSideBC,
              const std::vector<size_t> *requiredSearchEdges) {
    BoundarySplit result;
    for (const auto &info : boundaryComponents) {
        if (info.component == searchSideBC) {
            if (requiredSearchEdges && info.edgeIndices != *requiredSearchEdges)
                result.searchSideRejected = true;
            else
                result.searchCurveCount = info.curveNames.size();
            continue;
        }

        std::optional<std::string> name =
            info.curveNames.size() == 1
                ? std::optional<std::string>(info.curveNames.front())
                : info.linkName;
        if (!name) {
            // describeBoundary_() never leaves a multi-curve component
            // without a linkName. Reported rather than skipped: skipping
            // would silently turn a cobordism into a "direct" witness.
            result.unnamedSide = true;
            continue;
        }

        result.otherSides.push_back(
            {*name, static_cast<int>(info.curveNames.size())});
    }
    return result;
}

OrientationVerdict classifyRowOrientation(
    const RowOrientation &row, const std::vector<OrientedCurve> &curves,
    const std::map<const regina::Edge<3> *, size_t> &surfaceComponentOf) {
    std::unordered_map<size_t, bool> componentMatch;
    for (const OrientedCurve &curve : curves) {
        if (curve.empty())
            continue;

        std::optional<bool> curveMatch;
        for (const OrientedEdge &oe : curve) {
            auto it = row.tailOf.find(oe.edge->index());
            if (it == row.tailOf.end())
                return OrientationVerdict::foreignEdge;
            const regina::Vertex<3> *tail =
                oe.reversed ? oe.edge->vertex(1) : oe.edge->vertex(0);
            bool edgeMatches = tail->index() == it->second;
            if (!curveMatch)
                curveMatch = edgeMatches;
            else if (*curveMatch != edgeMatches)
                return OrientationVerdict::incoherentCurve;
        }

        auto comp = surfaceComponentOf.find(curve.front().edge);
        if (comp == surfaceComponentOf.end())
            return OrientationVerdict::incoherentCurve;
        auto [slot, inserted] = componentMatch.emplace(comp->second, *curveMatch);
        if (!inserted && slot->second != *curveMatch)
            return OrientationVerdict::mismatch;
    }
    return componentMatch.empty() ? OrientationVerdict::mismatch
                                  : OrientationVerdict::match;
}

} // namespace cobordismgraph

namespace rowsearch {

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

namespace {

// Which row components and how many far-side curves each surface component
// carries, as a canonical string: surface components are unlabelled, so the
// per-component entries are sorted.
std::string groupingOf(const farside::OutgoingLink &link,
                       const farside::WitnessRedrawer &row) {
  std::map<size_t, std::pair<std::vector<size_t>, int>> bySurface;
  for (size_t i = 0; i < link.incomingFirstEdge.size(); ++i)
    bySurface[link.incomingSurfaceComponent[i]].first.push_back(
        row.rowComponentOf(link.incomingFirstEdge[i]));
  for (size_t sc : link.surfaceComponent) ++bySurface[sc].second;
  std::vector<std::string> parts;
  for (auto &[sc, entry] : bySurface) {
    std::sort(entry.first.begin(), entry.first.end());
    std::string s;
    for (size_t rc : entry.first) s += std::to_string(rc) + '.';
    parts.push_back(s + ':' + std::to_string(entry.second));
  }
  std::sort(parts.begin(), parts.end());
  std::string out;
  for (const std::string &p : parts) out += p + '|';
  return out;
}

} // namespace

std::string keptKey(const cobordismgraph::Witness &w, const farside::OutgoingLink &link,
                    const farside::WitnessRedrawer &row) {
    return cobordismgraph::witnessIdentity(w) + '\x1f' + groupingOf(link, row);
}

} // namespace rowsearch

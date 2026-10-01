//
//  preconditions.cpp
//

#include "cobound/search/preconditions.h"

#include <optional>
#include <unordered_map>

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

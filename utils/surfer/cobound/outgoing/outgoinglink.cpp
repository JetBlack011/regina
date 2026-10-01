//
//  outgoinglink.cpp
//

#include "cobound/outgoing/outgoinglink.h"

#include "cobound/search/preconditions.h"

#include <algorithm>

namespace farside {

OutgoingCurve outgoingCurve(const OrientedCurve &curve) {
    OutgoingCurve out;
    out.reserve(curve.size());
    for (const OrientedEdge &oe : curve)
        out.push_back({oe.edge, oe.reversed});
    return out;
}

std::optional<std::map<size_t, int>> incomingFlips(
    const cobordismgraph::RowOrientation &row,
    const std::vector<OrientedCurve> &incomingCurves,
    const std::map<const regina::Edge<3> *, size_t> &surfaceComponentOf) {
    return cobordismgraph::judgeRowOrientation(row, incomingCurves, surfaceComponentOf)
        .consistentFlips();
}

std::optional<OutgoingLink> orientedOutgoingLink(
    const KnottedSurface &surface, const OutgoingMap &map,
    const cobordismgraph::RowOrientation &row, size_t incomingBC) {
    return orientedOutgoingLink(surface.orientedBoundaryLinks(),
                                surface.boundaryEdgeSurfaceComponent(), map,
                                row, incomingBC);
}

std::optional<OutgoingLink> orientedOutgoingLink(
    const std::vector<std::pair<size_t, std::vector<OrientedCurve>>> &oriented,
    const std::map<const regina::Edge<3> *, size_t> &surfaceOf,
    const OutgoingMap &map, const cobordismgraph::RowOrientation &row,
    size_t incomingBC, std::string *why) {
    std::optional<std::map<size_t, int>> flips;
    for (const auto &[bc, curves] : oriented)
        if (bc == incomingBC) flips = incomingFlips(row, curves, surfaceOf);
    if (!flips) {
        if (why) *why = "incoming orientation is inconsistent";
        return std::nullopt;
    }
    return orientedOutgoingLink(oriented, surfaceOf, map, *flips, incomingBC, why);
}

std::optional<OutgoingLink> orientedOutgoingLink(
    const std::vector<std::pair<size_t, std::vector<OrientedCurve>>> &oriented,
    const std::map<const regina::Edge<3> *, size_t> &surfaceOf,
    const OutgoingMap &map, const std::map<size_t, int> &flips, size_t incomingBC,
    std::string *why) {
    OutgoingLink out;
    for (const auto &[bc, curves] : oriented) {
        if (bc != incomingBC) continue;
        for (const OrientedCurve &curve : curves) {
            if (curve.empty()) continue;
            out.incomingFirstEdge.push_back(curve.front().edge->index());
            out.incomingSurfaceComponent.push_back(surfaceOf.at(curve.front().edge));
        }
    }
    for (const auto &[bc, curves] : oriented) {
        if (bc != map.boundaryComponent()) continue;
        for (const OrientedCurve &curve : curves) {
            if (curve.empty()) continue;
            size_t comp = surfaceOf.at(curve.front().edge);
            auto f = flips.find(comp);
            if (f == flips.end()) { // cannot happen: every component meets the row
                if (why) *why = "a surface component misses the row";
                return std::nullopt;
            }
            knotbuilder::EdgeCycle cyc = map.carry(outgoingCurve(curve));
            if (f->second < 0) {
                std::reverse(cyc.begin(), cyc.end());
                for (auto &de : cyc) de.reversed = !de.reversed;
            }
            out.curves.push_back(std::move(cyc));
            out.surfaceComponent.push_back(comp);
        }
    }
    return out;
}

} // namespace farside

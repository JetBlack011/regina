//
//  outgoinglink.cpp
//

#include "cobound/outgoing/outgoinglink.h"

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
    std::map<size_t, int> flips;
    for (const OrientedCurve &curve : incomingCurves) {
        if (curve.empty()) continue;
        std::optional<bool> match;
        for (const OrientedEdge &oe : curve) {
            auto it = row.tailOf.find(oe.edge->index());
            if (it == row.tailOf.end()) return std::nullopt;
            const regina::Vertex<3> *tail =
                oe.reversed ? oe.edge->vertex(1) : oe.edge->vertex(0);
            bool m = tail->index() == it->second;
            if (match && *match != m) return std::nullopt;
            match = m;
        }
        auto comp = surfaceComponentOf.find(curve.front().edge);
        if (comp == surfaceComponentOf.end()) return std::nullopt;
        int flip = *match ? 1 : -1;
        auto [slot, fresh] = flips.emplace(comp->second, flip);
        if (!fresh && slot->second != flip) return std::nullopt;
    }
    return flips;
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
    size_t incomingBC) {
    std::optional<std::map<size_t, int>> flips;
    for (const auto &[bc, curves] : oriented)
        if (bc == incomingBC) flips = incomingFlips(row, curves, surfaceOf);
    if (!flips) return std::nullopt;

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
            auto f = flips->find(comp);
            if (f == flips->end()) return std::nullopt; // cannot happen: every
                                                        // component meets the row
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

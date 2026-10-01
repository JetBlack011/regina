//
//  farsidecurves.cpp
//

#include "cobound/outgoing/farsidecurves.h"

#include <algorithm>

#include "diagramtriangulation/thickening/prism.h"

namespace farside {

OutgoingMap::OutgoingMap(const regina::Triangulation<3> &knotT,
                         const CobordismBuilder<3> &cob) {
    const regina::Triangulation<3> &base = cob.baseTriangulation();
    const regina::Triangulation<4> &w = cob.getCobordism();
    if (w.countBoundaryComponents() != 2)
        throw regina::InvalidArgument(
            "OutgoingMap: the thickening must have two boundary components");
    const size_t incoming = cob.baseBoundaryComponent()->index();
    bc_ = incoming == 0 ? 1 : 0;
    const regina::BoundaryComponent<4> *bc = w.boundaryComponent(bc_);

    // CobordismBuilder orders its own copy of knotT, relabelling vertices
    // within tetrahedra but keeping tetrahedron indices: recover each
    // tetrahedron's relabelling as the isomorphism base -> knotT that fixes
    // every tetrahedron.
    std::optional<regina::Isomorphism<3>> relabel;
    base.findAllIsomorphisms(knotT, [&](const regina::Isomorphism<3> &iso) {
        for (size_t i = 0; i < base.size(); ++i)
            if (iso.simpImage(i) != i) return false;
        relabel = iso;
        return true;
    });
    if (!relabel)
        throw regina::InvalidArgument(
            "OutgoingMap: the thickening's base is not knotT relabelled");

    std::unordered_map<const regina::Edge<4> *, size_t> localEdge;
    for (size_t k = 0; k < bc->countEdges(); ++k) localEdge.emplace(bc->edge(k), k);
    std::unordered_map<const regina::Vertex<4> *, size_t> localVertex;
    for (size_t k = 0; k < bc->countVertices(); ++k)
        localVertex.emplace(bc->vertex(k), k);

    for (size_t i = 0; i < base.size(); ++i) {
        const regina::Simplex<4> *top = cob.currentTopSimplex(base.tetrahedron(i), 3);
        const regina::Perm<4> p = relabel->facetPerm(i);
        const regina::Tetrahedron<3> *t = knotT.tetrahedron(i);
        for (int a = 0; a < 4; ++a) {
            const regina::Vertex<4> *wv =
                top->vertex(SimplicialPrism<4>::localVertex(a, true));
            auto lv = localVertex.find(wv);
            if (lv == localVertex.end())
                throw regina::InvalidArgument(
                    "OutgoingMap: a top vertex is not on the outgoing boundary");
            vertexToT_[lv->second] = t->vertex(p[a])->index();
            for (int b = a + 1; b < 4; ++b) {
                const regina::Edge<4> *we =
                    top->edge(SimplicialPrism<4>::localVertex(a, true),
                              SimplicialPrism<4>::localVertex(b, true));
                auto le = localEdge.find(we);
                if (le == localEdge.end())
                    throw regina::InvalidArgument(
                        "OutgoingMap: a top edge is not on the outgoing boundary");
                edgeToT_[le->second] = t->edge(p[a], p[b])->index();
            }
        }
    }
    if (edgeToT_.size() != bc->countEdges() || vertexToT_.size() != bc->countVertices())
        throw regina::InvalidArgument(
            "OutgoingMap: the top prisms do not cover the outgoing boundary");
    tTail_.resize(knotT.countEdges());
    for (size_t e = 0; e < knotT.countEdges(); ++e)
        tTail_[e] = knotT.edge(e)->vertex(0)->index();
}

knotbuilder::EdgeCycle OutgoingMap::carry(const OrientedCurve &curve) const {
    knotbuilder::EdgeCycle out;
    out.reserve(curve.size());
    for (const OrientedEdge &oe : curve) {
        size_t e = edgeToT_.at(oe.edge->index());
        size_t tail =
            vertexToT_.at((oe.reversed ? oe.edge->vertex(1) : oe.edge->vertex(0))->index());
        out.push_back({e, tTail_[e] != tail});
    }
    return out;
}

knotbuilder::EdgeCycle OutgoingMap::carryCycle(
    const std::vector<const regina::Edge<3> *> &edges) const {
    // Chain the edges by shared vertices: each vertex of a closed curve
    // meets exactly two of them.
    std::unordered_map<size_t, std::vector<size_t>> at; // vertex -> edge positions
    for (size_t i = 0; i < edges.size(); ++i) {
        at[edges[i]->vertex(0)->index()].push_back(i);
        at[edges[i]->vertex(1)->index()].push_back(i);
    }
    for (const auto &[v, es] : at)
        if (es.size() != 2)
            throw regina::InvalidArgument("carryCycle: not a single closed curve");
    OrientedCurve curve;
    curve.reserve(edges.size());
    std::vector<bool> used(edges.size(), false);
    size_t i = 0;
    size_t tail = edges[0]->vertex(0)->index();
    for (size_t k = 0; k < edges.size(); ++k) {
        used[i] = true;
        bool reversed = edges[i]->vertex(0)->index() != tail;
        curve.push_back({edges[i], reversed});
        size_t head = edges[i]->vertex(reversed ? 0 : 1)->index();
        const auto &next = at.at(head);
        size_t j = next[0] == i ? next[1] : next[0];
        if (k + 1 < edges.size() && used[j])
            throw regina::InvalidArgument("carryCycle: not a single closed curve");
        i = j;
        tail = head;
    }
    return carry(curve);
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
            knotbuilder::EdgeCycle cyc = map.carry(curve);
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

std::vector<size_t> boundaryEdgesOf(const regina::Triangulation<4> &tri,
                                    const std::vector<int> &faces,
                                    size_t bcIndex) {
    const regina::BoundaryComponent<4> *bc = tri.boundaryComponent(bcIndex);
    std::unordered_map<const regina::Edge<4> *, size_t> local;
    for (size_t k = 0; k < bc->countEdges(); ++k) local.emplace(bc->edge(k), k);
    std::vector<size_t> edges;
    for (int f : faces) {
        const regina::Triangle<4> *t = tri.triangle(f);
        for (int i = 0; i < 3; ++i)
            if (auto it = local.find(t->edge(i)); it != local.end())
                edges.push_back(it->second);
    }
    std::ranges::sort(edges);
    edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
    return edges;
}

} // namespace farside

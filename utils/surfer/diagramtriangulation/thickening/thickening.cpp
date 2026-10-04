//
//  cobordismbuilder.cpp
//
//  Created by John Teague on 04/12/2025.
//

#include "diagramtriangulation/thickening/thickening.h"

#include <algorithm>
#include <cassert>
#include <optional>
#include <unordered_map>
#include <unordered_set>

#include "diagramtriangulation/thickening/collar.h"
#include "linknaming/complement/edgecycles.h"

template <int dim>
CobordismBuilder<dim>::CobordismBuilder(const regina::Triangulation<dim> &tri)
    : tri_(tri) {
    if (isOrdered(tri_))
        return;

    if constexpr (dim == 3) {
        if (!tri_.order() || !isOrdered(tri_))
            throw regina::InvalidArgument(
                "CobordismBuilder::CobordismBuilder(): triangulation "
                "could not be ordered.");
    } else {
        throw regina::InvalidArgument(
            "CobordismBuilder::CobordismBuilder(): triangulation is not "
            "ordered, and automatic ordering is only implemented for "
            "dim == 3.");
    }
}

template <int dim>
bool CobordismBuilder<dim>::isOrdered(const regina::Triangulation<dim> &tri) {
    for (const auto &s : tri.simplices()) {
        for (int f = 0; f <= dim; ++f) {
            if (s->adjacentSimplex(f) == nullptr)
                continue;
            regina::Perm<dim + 1> g = s->adjacentGluing(f);
            std::vector<int> a;
            for (int i = 0; i < dim + 1; ++i) {
                if (i == f)
                    continue;
                a.push_back(g[i]);
            }

            if (!std::ranges::is_sorted(a)) {
                return false;
            }
        }
    }

    return true;
}

template <int dim>
template <int d>
regina::Triangulation<d> &CobordismBuilder<dim>::glueBoundaries(
    regina::Triangulation<d> &tri, int bdryIndex1, int bdryIndex2,
    const regina::Isomorphism<d - 1> &iso) {
    using Gluing = std::tuple<regina::Simplex<d> *, int,
                              regina::Simplex<d> *, regina::Perm<d + 1>>;
    // Defer joins to avoid modifying the triangulation while iterating
    // boundary components.
    std::vector<Gluing> gluings;
    const regina::BoundaryComponent<d> *bdry1 =
        tri.boundaryComponent(bdryIndex1);
    const regina::BoundaryComponent<d> *bdry2 =
        tri.boundaryComponent(bdryIndex2);

    for (int i = 0; i < bdry1->size(); ++i) {
        // Since these are boundary faces, they belong to exactly 1 simplex
        const auto &emb1 = bdry1->facet(i)->front();
        const auto &emb2 = bdry2->facet(iso.simpImage(i))->front();
        regina::Simplex<d> *s1 = emb1.simplex();
        regina::Simplex<d> *s2 = emb2.simplex();

        regina::Perm<d + 1> i1 = emb1.vertices();
        regina::Perm<d + 1> i2 = emb2.vertices();
        regina::Perm<d> p = iso.facetPerm(i);
        std::array<int, d + 1> isoPerm;

        for (int j = 0; j < d; ++j) {
            isoPerm[j] = p[j];
        }
        isoPerm[d] = d;

        gluings.emplace_back(s1, emb1.face(), s2,
                             i2 * isoPerm * i1.inverse());
    }

    for (const auto &[src, srcFacet, dst, gluingPerm] : gluings) {
        src->join(srcFacet, dst, gluingPerm);
    }

    return tri;
}

template <int dim>
template <int d>
regina::Triangulation<d> CobordismBuilder<dim>::glueTriangulations(
    const regina::Triangulation<d> &tri1, int bdryIndex1,
    const regina::Triangulation<d> &tri2, int bdryIndex2,
    const regina::Isomorphism<d - 1> &iso) {
    regina::Triangulation<d> tri;
    tri.insertTriangulation(tri1);
    tri.insertTriangulation(tri2);

    return glueBoundaries<d>(tri, bdryIndex1,
                             tri1.countBoundaryComponents() + bdryIndex2, iso);
}

template <int dim>
regina::Triangulation<dim + 1> &CobordismBuilder<dim>::thicken_() {
    // Make a new prism for each simplex in the triangulation. This is
    // its own layer: it gets fully glued together internally below,
    // independently of any previous layer.
    PrismMap newPrisms;
    newPrisms.reserve(tri_.size());
    for (const auto *s : tri_.simplices()) {
        newPrisms.emplace(s, cob_);
    }

    // Now glue the prisms together along their walls according to the
    // gluing of the original triangulation
    std::set<std::pair<const regina::Simplex<dim> *, int>> visited;

    for (size_t i = 0; i < tri_.size(); ++i) {
        for (int facet = 0; facet < dim + 1; ++facet) {
            const regina::Simplex<dim> *s = tri_.simplex(i);
            const regina::Simplex<dim> *adj = s->adjacentSimplex(facet);

            if (adj == nullptr || visited.contains({s, facet}))
                continue;

            int adjFacet = s->adjacentFacet(facet);

            newPrisms.at(s).glue(facet, newPrisms.at(adj), adjFacet);

            visited.insert({s, facet});
            visited.insert({adj, adjFacet});
        }
    }

    // If there is a previous layer, stitch this layer's bottom onto its
    // top, simplex by simplex: since both layers thicken the same base
    // triangulation, the seam uses the same base simplex on both sides
    // and needs no relabelling of base vertices.
    if (hasPreviousLayer_) {
        for (const auto *s : tri_.simplices()) {
            topPrisms_.at(s).stitchTop(newPrisms.at(s));
        }
    } else {
        // First layer only: see baseBoundaryComponent()'s doc comment.
        // simplex(0) of any base simplex's prism holds every one of that
        // base simplex's own "bottom" vertices plus the "top" copy of
        // base vertex 0 -- captured here since this is the only point at
        // which the *first* layer's own prisms are directly at hand
        // (topPrisms_ is overwritten wholesale by every later call).
        const regina::Simplex<dim> *anyBase = *tri_.simplices().begin();
        baseBoundaryFacetSimplex_ = newPrisms.at(anyBase).simplex(0);
    }

    topPrisms_ = std::move(newPrisms);
    hasPreviousLayer_ = true;

    // Validity is guaranteed by construction (a layer of prisms over a
    // valid, ordered triangulation, glued face-for-face onto the previous
    // layer), and checking it here is the dominant cost of thicken() (a
    // full vertex-link recognition pass over the whole accumulated
    // cobordism, on every single layer). Debug-only.
    assert(cob_.isValid() &&
           "CobordismBuilder::thicken(): resulting triangulation is not "
           "valid.");

    return cob_;
}

template <int dim>
regina::BoundaryComponent<dim + 1> *
CobordismBuilder<dim>::baseBoundaryComponent() const {
    if (!baseBoundaryFacetSimplex_)
        throw regina::InvalidArgument(
            "CobordismBuilder::baseBoundaryComponent(): no thicken() call "
            "has been made");
    return baseBoundaryFacetSimplex_->template face<dim>(dim + 1)
        ->boundaryComponent();
}

template class CobordismBuilder<2>;
template class CobordismBuilder<3>;

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

knotbuilder::EdgeCycle OutgoingMap::carry(const OutgoingCurve &curve) const {
    knotbuilder::EdgeCycle out;
    out.reserve(curve.size());
    for (const OutgoingEdge &oe : curve) {
        size_t e = edgeToT_.at(oe.edge->index());
        size_t tail =
            vertexToT_.at((oe.reversed ? oe.edge->vertex(1) : oe.edge->vertex(0))->index());
        out.push_back({e, tTail_[e] != tail});
    }
    return out;
}

knotbuilder::EdgeCycle OutgoingMap::carryCycle(
    const std::vector<const regina::Edge<3> *> &edges) const {
    // One simple closed curve, run from the first edge's vertex(0).
    const auto steps = edgecycles::walkClosedCurve(edgecycles::endsOf(edges));
    if (!steps)
        throw regina::InvalidArgument("carryCycle: not a single closed curve");
    OutgoingCurve curve;
    curve.reserve(steps->size());
    for (const edgecycles::Step &s : *steps)
        curve.push_back({edges[s.pos], s.reversed});
    return carry(curve);
}

namespace {

// How many closed curves a link's edges form (edgecycles::countClosedCurves):
// the count Link::countComponents() (linkcomplement.h) gives, without
// building the edge-set Link. Edges that are not disjoint closed curves --
// a repeated edge, as that Link refuses, or a vertex meeting other than two
// -- are refused.
int countLinkComponents(const std::vector<const regina::Edge<3> *> &edges) {
    const std::optional<size_t> curves =
        edgecycles::countClosedCurves(edgecycles::endsOf(edges));
    if (!curves)
        throw regina::InvalidArgument(
            "buildAmbient(): the link's edges are not disjoint closed curves");
    return static_cast<int>(*curves);
}

} // namespace

void buildAmbient(const std::string &pdNotation, int thickenLayers,
                  int collarLayers, ThickenedLink &row) {
    row.pdcode = knotbuilder::parsePDCode(pdNotation);
    row.link = knotbuilder::buildLink(row.pdcode);

    auto &[t2, edges2, reversed2] = row.link;
    row.componentCount = countLinkComponents(edges2);

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

    row.incomingBC = cob.baseBoundaryComponent()->index();
    row.tri = cob.getCobordism();

    if (collarLayers > 0) {
        for (regina::Triangle<4> *t : collarBuilder.resolve())
            row.seedFaces.push_back(static_cast<int>(t->index()));
        // In index order, never the set's (address) order: a surface's
        // triangles are numbered from the seed's, so each component's
        // orientation, and so where each of its curves starts, follows it.
        std::ranges::sort(row.seedFaces);
    }
}

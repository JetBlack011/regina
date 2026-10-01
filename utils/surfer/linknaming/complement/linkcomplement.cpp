//
//  linkcomplement.cpp
//
//  Created by John Teague on 07/21/2026.
//

#include "linknaming/complement/linkcomplement.h"

#include "linknaming/complement/edgecycles.h"

#include <algorithm>
#include <map>
#include <set>

#include <triangulation/dim3/homologicaldata.h>

std::atomic<bool> simplifyComplements{true};

EdgeComplement::EdgeComplement(
    const regina::Triangulation<3> &tri,
    const std::vector<const regina::Edge<3> *> &edges)
    : tri_(&tri), edges_(edges) {}

std::map<size_t, std::set<int>> EdgeComplement::tetEdges_() const {
    // Keyed by tetrahedron INDEX, never by address: the drilling loops below
    // take the lowest-index tetrahedron's lowest local edge next, so they
    // pinch in an order that depends on the edges alone.
    std::map<size_t, std::set<int>> map;
    for (const regina::Edge<3> *edge : edges_)
        for (const regina::EdgeEmbedding<3> &emb : edge->embeddings())
            map[emb.tetrahedron()->index()].insert(emb.face());
    return map;
}

namespace {

// Pinches every edge `tetEdges` lists (tetrahedron index -> local edges) out
// of `complement`, lowest tetrahedron first. pinchEdge() only ever appends
// tetrahedra, never removes or renumbers one, so the indices stay valid.
void pinchAll(regina::Triangulation<3> &complement,
              std::map<size_t, std::set<int>> tetEdges) {
    while (!tetEdges.empty()) {
        auto &[tet, edges] = *tetEdges.begin();
        regina::Edge<3> *e = complement.tetrahedron(tet)->edge(*edges.begin());

        for (const regina::EdgeEmbedding<3> &emb : e->embeddings()) {
            auto it = tetEdges.find(emb.tetrahedron()->index());
            if (it == tetEdges.end())
                continue;
            it->second.erase(emb.face());
            if (it->second.empty())
                tetEdges.erase(it);
        }
        complement.pinchEdge(e);
    }
}

} // namespace

EdgeComplement &EdgeComplement::operator=(const EdgeComplement &other) {
    if (this != &other) {
        tri_ = other.tri_;
        edges_ = other.edges_;
    }
    return *this;
}

regina::Triangulation<3> EdgeComplement::buildComplement() const {
    regina::Triangulation<3> complement(*tri_);
    pinchAll(complement, tetEdges_());

    // complement.idealToFinite();
    if (simplifyComplements.load(std::memory_order_relaxed))
        complement.simplify();

    return complement;
}

std::pair<regina::Triangulation<3>, std::vector<const regina::Edge<3> *>>
EdgeComplement::drillTrackingEdges_(
    const std::vector<const regina::Edge<3> *> &trackEdges) const {
    regina::Triangulation<3> complement(*tri_);

    // trackEdges is tracked via (Tetrahedron*, localEdge) descriptors,
    // which stay valid across pinchEdge() (pinchEdge() only ever appends
    // new tetrahedra, never removes or renumbers existing ones).
    std::vector<std::pair<regina::Tetrahedron<3> *, int>> trackedDescs;
    trackedDescs.reserve(trackEdges.size());
    for (const regina::Edge<3> *e : trackEdges) {
        auto emb = e->front();
        trackedDescs.emplace_back(
            complement.tetrahedron(emb.tetrahedron()->index()), emb.edge());
    }

    // Drill this object's own edges_ -- as buildComplement() does, just
    // without simplify(): homology doesn't need simplification, and
    // simplifying would make tracking trackedDescs through it intractable.
    pinchAll(complement, tetEdges_());

    std::vector<const regina::Edge<3> *> tracked;
    tracked.reserve(trackedDescs.size());
    for (const auto &[tet, localEdge] : trackedDescs)
        tracked.push_back(tet->edge(localEdge));

    return {std::move(complement), std::move(tracked)};
}

long EdgeComplement::linkingNumberWith(const EdgeComplement &other) const {
    auto [complement, trackedOther] = drillTrackingEdges_(other.edges_);

    // Drilling curve A collapses it to a single ideal vertex -- a torus
    // cusp. H_1 must therefore be computed with that vertex *truncated*
    // (a genuine torus boundary), not treated as an ordinary 0-cell:
    // coning the cusp torus to a point kills every loop on it, including
    // A's meridian, which is precisely the class the linking number is
    // read off. A hand-rolled (vertices x edges, edges x triangles) pair
    // of boundary maps does exactly that wrong thing, and used to make
    // this routine return 0 for every link it was given.
    //
    // HomologicalData computes the correctly-truncated group and, unlike
    // Triangulation<3>::homology(), keeps it marked so a specific cycle
    // can be located inside it. Its standard cellular coordinates order
    // the 1-cells as "edges.begin() to edges.end(), followed by the ideal
    // edges of faces" (see its class documentation), so the triangulation's
    // own edges are a *prefix* of the coordinate space -- which is what
    // lets curve B's edge-indexed cycle vector below be used as-is, padded
    // with zeros over the trailing ideal-truncation cells.
    regina::HomologicalData hd(complement);
    const regina::MarkedAbelianGroup &h1 = hd.homology(1);

    // H_1 of a knot complement in S^3 is always Z, so this should now be
    // unreachable; kept as a defensive check rather than an assertion,
    // since reporting "no detected linking" stays sound (this routine
    // never falsely reports linking) if the precondition is ever violated.
    if (!h1.isZ())
        return 0;

    // Curve B must avoid the drilled cusp for the zero-padding described
    // above to be valid: an edge ending at the ideal vertex terminates on
    // one of the *new* truncation 0-cells instead of an ordinary vertex,
    // so such a cycle would no longer close up in these coordinates and
    // snfRep() would throw rather than return a wrong answer.
    // Geometrically this cannot happen -- an edge of B running into the
    // ideal vertex would mean B meets the drilled-out curve A,
    // contradicting this method's disjointness precondition -- so it is
    // checked, not assumed.
    for (const regina::Edge<3> *e : trackedOther)
        if (e->vertex(0)->isIdeal() || e->vertex(1)->isIdeal())
            throw regina::InvalidArgument(
                "EdgeComplement::linkingNumberWith(): other's curve meets "
                "the drilled ideal vertex, so the two curves are not "
                "disjoint");

    // Walk trackedOther's edges into an oriented cycle
    // (edgecycles::walkClosedCurve(), the walk Link::Link() also splits a
    // multi-component edge set with), then read off its class in H_1. Only
    // the magnitude is meaningful, so no canonical orientation needs to be
    // imposed: any consistent walk direction works. The walk must be
    // *oriented*: an unsigned sum of B's edges is not a cycle in these
    // coordinates and snfRep() would reject it.
    regina::Vector<regina::Integer> cycle(hd.countStandardCells(1));
    if (!trackedOther.empty()) {
        const auto steps =
            edgecycles::walkClosedCurve(edgecycles::endsOf(trackedOther));
        if (!steps)
            throw regina::InvalidArgument(
                "EdgeComplement::linkingNumberWith(): other's edges do "
                "not form a single closed curve");
        for (const edgecycles::Step &s : *steps)
            cycle[trackedOther[s.pos]->index()] += s.reversed ? -1 : 1;
    }

    return h1.snfRep(cycle)[0].abs().safeValue<long>();
}

std::vector<size_t> EdgeComplement::edgeIndices() const {
    std::vector<size_t> indices;
    indices.reserve(edges_.size());
    for (const regina::Edge<3> *e : edges_)
        indices.push_back(e->index());
    std::ranges::sort(indices);
    return indices;
}

bool operator<(const EdgeComplement &e1, const EdgeComplement &e2) {
    return e1.edges_.size() < e2.edges_.size();
}

std::ostream &operator<<(std::ostream &os, const EdgeComplement &e) {
    os << "{";
    for (int i = 0; i < e.edges_.size(); ++i) {
        os << e.edges_[i]->index();
        if (i != e.edges_.size() - 1) {
            os << ", ";
        }
    }
    return os << "}";
}

Link::Link(const regina::Triangulation<3> &tri,
          const std::vector<const regina::Edge<3> *> &edges)
    : EdgeComplement(tri, edges) {
    // Add the edge to the correct component (edgecycles::walkCurves()). The
    // walk takes its edges in index order, never in address order: each
    // component starts from its lowest-index edge, towards that edge's
    // vertex(1), and the components come out in the order of their lowest
    // edges. So the components and their edge sequences are a function of
    // the edge set alone.
    std::vector<const regina::Edge<3> *> sorted(edges.begin(), edges.end());
    std::ranges::sort(sorted, {}, [](const regina::Edge<3> *e) {
        return e->index();
    });
    if (std::ranges::adjacent_find(sorted) != sorted.end()) {
        throw regina::InvalidArgument(
            "Link::Link: Duplicate edges in link");
    }

    for (const auto &curve : edgecycles::walkCurves(edgecycles::endsOf(sorted))) {
        std::vector<const regina::Edge<3> *> compEdges;
        compEdges.reserve(curve.size());
        for (const edgecycles::Step &s : curve)
            compEdges.push_back(sorted[s.pos]);
        comps_.emplace_back(tri, compEdges);
    }
}

Link &Link::operator=(const Link &other) {
    if (this != &other) {
        EdgeComplement::operator=(other);
        comps_ = other.comps_;
    }
    return *this;
}

bool operator<(const Link &l1, const Link &l2) {
    if (l1.comps_.size() != l2.comps_.size())
        return l1.comps_.size() < l2.comps_.size();

    return static_cast<EdgeComplement>(l1) <
           static_cast<EdgeComplement>(l2);
}

std::ostream &operator<<(std::ostream &os, const Link &l) {
    os << "[";
    for (int i = 0; i < l.comps_.size(); ++i) {
        os << l.comps_[i];
        if (i != l.comps_.size() - 1) {
            os << ", ";
        }
    }
    return os << "]";
}

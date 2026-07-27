//
//  linkcomplement.cpp
//
//  Created by John Teague on 07/21/2026.
//

#include "linkcomplement.h"

#include <iostream>
#include <list>
#include <optional>

#include <triangulation/dim3/homologicaldata.h>

std::mutex censusLookupMutex;

namespace {

// Memoizes regina::Census::lookup() by isomorphism signature. The same
// boundary complement (e.g. a hyperbolic knot guaranteed to be a census
// hit) tends to recur across many found surfaces, and every real lookup
// reopens six on-disk census databases from scratch under
// censusLookupMutex -- so caching turns "one lookup per surface" into "one
// lookup per distinct complement", which is where nearly all of the
// mutex's contention actually comes from. Guarded by censusLookupMutex
// itself, which already has to serialize the underlying lookups.
std::unordered_map<std::string, std::optional<std::string>> censusNameCache;

// Returns the census name for `complement` (whose isoSig is `sig`), or
// nullopt if it has no census hit. Must be called with the fast
// recogniseHandlebody() path already ruled out by the caller.
std::optional<std::string> cachedCensusLookupName(
    const regina::Triangulation<3> &complement, const std::string &sig) {
    std::lock_guard<std::mutex> lock(censusLookupMutex);

    auto it = censusNameCache.find(sig);
    if (it != censusNameCache.end())
        return it->second;

    std::list<regina::CensusHit> hits = regina::Census::lookup(complement);
    std::optional<std::string> name =
        hits.empty() ? std::nullopt
                     : std::make_optional(hits.front().name());
    censusNameCache.emplace(sig, name);
    return name;
}

} // namespace

EdgeComplement::EdgeComplement(
    const regina::Triangulation<3> &tri,
    const std::vector<const regina::Edge<3> *> &edges)
    : tri_(&tri), edges_(edges) {
    // Add the edge to the correct tetrahedra
    for (const regina::Edge<3> *edge : edges) {
        for (const regina::EdgeEmbedding<3> &emb : edge->embeddings()) {
            tetEdges_[emb.tetrahedron()].insert(emb.face());
        }
    }
}

EdgeComplement &EdgeComplement::operator=(const EdgeComplement &other) {
    if (this != &other) {
        tri_ = other.tri_;
        edges_ = other.edges_;
        tetEdges_ = other.tetEdges_;
    }
    return *this;
}

regina::Triangulation<3> EdgeComplement::buildComplement() const {
    regina::Triangulation<3> complement(*tri_);
    std::unordered_map<regina::Tetrahedron<3> *, std::unordered_set<size_t>>
        complementTetEdges;

    for (const auto &[tet, edges] : tetEdges_) {
        complementTetEdges.emplace(complement.tetrahedron(tet->index()),
                                   edges);
    }

    while (!complementTetEdges.empty()) {
        auto &[tet, edges] = *complementTetEdges.begin();
        regina::Edge<3> *e = tet->edge(*edges.begin());

        for (const regina::EdgeEmbedding<3> &emb : e->embeddings()) {
            regina::Tetrahedron<3> *embTet = emb.tetrahedron();
            complementTetEdges[embTet].erase(emb.face());
            if (complementTetEdges[embTet].empty()) {
                complementTetEdges.erase(embTet);
            }
        }
        complement.pinchEdge(e);
    }

    // complement.idealToFinite();
    complement.simplify();

    return complement;
}

bool EdgeComplement::recognizeComplement() const {
    auto complement = buildComplement();
    if (complement.recogniseHandlebody() == 1) {
        std::cout << "      unknot, " << complement.isoSig() << "\n";
        return true;
    }

    std::string sig = complement.isoSig();
    if (auto name = cachedCensusLookupName(complement, sig)) {
        std::cout << "      recognized as " << *name << ", " << sig << "\n";
        return true;
    }
    return false;
}

std::string EdgeComplement::identify() const {
    auto complement = buildComplement();
    if (complement.recogniseHandlebody() == 1)
        return "Unknot";

    std::string sig = complement.isoSig();
    if (auto name = cachedCensusLookupName(complement, sig))
        return *name;
    return sig;
}

bool EdgeComplement::isUnknot() const {
    return buildComplement().recogniseHandlebody() == 1;
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

    // Drill this object's own edges_ -- identical loop to
    // buildComplement(), just without simplify(): homology doesn't need
    // simplification, and simplifying would make tracking trackedDescs
    // through it intractable.
    std::unordered_map<regina::Tetrahedron<3> *, std::unordered_set<size_t>>
        complementTetEdges;
    for (const auto &[tet, edges] : tetEdges_)
        complementTetEdges.emplace(complement.tetrahedron(tet->index()),
                                   edges);

    while (!complementTetEdges.empty()) {
        auto &[tet, edges] = *complementTetEdges.begin();
        regina::Edge<3> *e = tet->edge(*edges.begin());

        for (const regina::EdgeEmbedding<3> &emb : e->embeddings()) {
            regina::Tetrahedron<3> *embTet = emb.tetrahedron();
            complementTetEdges[embTet].erase(emb.face());
            if (complementTetEdges[embTet].empty()) {
                complementTetEdges.erase(embTet);
            }
        }
        complement.pinchEdge(e);
    }

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

    // Walk trackedOther's edges into an oriented cycle (the same
    // shared-vertex walk Link::Link() uses to split a multi-component edge
    // set into per-component knots), then read off its class in H_1. Only
    // the magnitude is meaningful, so no canonical orientation needs to be
    // imposed: any consistent walk direction works. The walk must be
    // *oriented*: an unsigned sum of B's edges is not a cycle in these
    // coordinates and snfRep() would reject it.
    regina::Vector<regina::Integer> cycle(hd.countStandardCells(1));
    if (!trackedOther.empty()) {
        std::vector<const regina::Edge<3> *> remaining(trackedOther.begin(),
                                                        trackedOther.end());
        const regina::Edge<3> *first = remaining.front();
        const regina::Vertex<3> *currVert = first->vertex(1);
        cycle[first->index()] += 1;
        remaining.erase(remaining.begin());

        while (!remaining.empty()) {
            bool found = false;
            for (size_t i = 0; i < remaining.size(); ++i) {
                const regina::Edge<3> *e = remaining[i];
                if (e->vertex(0) == currVert) {
                    cycle[e->index()] += 1;
                    currVert = e->vertex(1);
                } else if (e->vertex(1) == currVert) {
                    cycle[e->index()] -= 1;
                    currVert = e->vertex(0);
                } else {
                    continue;
                }
                remaining.erase(remaining.begin() +
                                static_cast<ptrdiff_t>(i));
                found = true;
                break;
            }
            if (!found)
                throw regina::InvalidArgument(
                    "EdgeComplement::linkingNumberWith(): other's edges do "
                    "not form a single closed curve");
        }
    }

    return h1.snfRep(cycle)[0].abs().safeValue<long>();
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
    // Add the edge to the correct component
    std::vector<std::vector<const regina::Edge<3> *>> edgesByComp;
    std::unordered_set<const regina::Edge<3> *> edgeSet(edges.begin(),
                                                        edges.end());
    if (edgeSet.size() != edges.size()) {
        throw regina::InvalidArgument(
            "Link::Link: Duplicate edges in link");
    }

    const regina::Vertex<3> *currVert;
    bool newComponent = true;
    while (!edgeSet.empty()) {
        if (newComponent) {
            const auto edge = *edgeSet.begin();
            edgesByComp.push_back({edge});
            currVert = edge->vertex(1);
            edgeSet.erase(edge);
        }

        newComponent = true;
        for (const regina::Edge<3> *edge : edgeSet) {
            if (edge->vertex(0) == currVert ||
                edge->vertex(1) == currVert) {
                edgesByComp.back().push_back(edge);
                currVert = edge->vertex(0) == currVert ? edge->vertex(1)
                                                       : edge->vertex(0);
                edgeSet.erase(edge);
                newComponent = false;
                break;
            }
        }
    }

    if (edgesByComp.size() == 1) {
        comps_.emplace_back(tri, edges);
        return;
    }

    for (const auto &compEdges : edgesByComp) {
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

void Link::recognizeComplement() const {
    if (EdgeComplement::recognizeComplement()) {
        return;
    }

    auto complement = buildComplement();
    ssize_t genus = complement.recogniseHandlebody();
    int numComponents = countComponents();

    if (genus != -1 && numComponents == 1) {
        std::cout << "[!] WARNING! Recognized as a genus " << genus
                  << " handlebody, " << complement.isoSig() << "\n";
        std::cout << "[!] This is almost definitely a bug, please "
                     "report it!\n";
    } else if (numComponents == 1) {
        std::cout << "      NOT unknot, " << complement.isoSig() << "\n";
    } else if (numComponents > 1) {
        for (int i = 0; i < numComponents; ++i) {
            regina::Triangulation<3> complement = buildComplement(i);
            ssize_t genus = complement.recogniseHandlebody();

            std::cout << "    Component " << i + 1 << ": ";
            if (genus == 1) {
                std::cout << "unknot, " << complement.isoSig() << "\n";
            } else if (genus != -1) {
                std::cout << "\n[!] WARNING! Recognized as a genus "
                          << genus << " handlebody, ";
                std::cout << "\n[!] This is almost definitely a bug, "
                             "please "
                             "report it!\n";
            } else {
                // Not any handlebody -- genuinely might be a census hit,
                // so (unlike the genus >= 0 cases above) this is the one
                // branch that actually needs the (cached) census lookup.
                std::string sig = complement.isoSig();
                if (auto name = cachedCensusLookupName(complement, sig)) {
                    std::cout << "      recognized as " << *name << ", "
                              << sig << "\n";
                } else {
                    std::cout << "NOT unknot, " << sig << "\n";
                }
            }
        }
    }
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

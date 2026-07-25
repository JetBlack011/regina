//
//  linkcomplement.cpp
//
//  Created by John Teague on 07/21/2026.
//

#include "linkcomplement.h"

#include <iostream>
#include <list>

std::mutex censusLookupMutex;

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

    std::list<regina::CensusHit> hits;
    {
        std::lock_guard<std::mutex> lock(censusLookupMutex);
        hits = regina::Census::lookup(complement);
    }
    if (!hits.empty()) {
        std::cout << "      recognized as " << hits.front().name() << ", "
                  << complement.isoSig() << "\n";
        return true;
    }
    return false;
}

std::string EdgeComplement::identify() const {
    auto complement = buildComplement();
    if (complement.recogniseHandlebody() == 1)
        return "Unknot";

    std::list<regina::CensusHit> hits;
    {
        std::lock_guard<std::mutex> lock(censusLookupMutex);
        hits = regina::Census::lookup(complement);
    }
    if (!hits.empty())
        return hits.front().name();
    return complement.isoSig();
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

    // M (vertices x edges, signed vertex incidence) and N (edges x
    // triangles, signed edge incidence) are boundary maps with
    // H_1 = ker(M)/im(N).
    //
    // Known limitation: drilling curve A leaves an ideal vertex, and
    // computing H_1 correctly in its presence requires treating that
    // vertex as truncated -- which this hand-rolled M/N does not do (and
    // Regina's own truncating routines would invalidate trackedOther's
    // tracked-edge descriptors). So h1.isZ() is used only as a confidence
    // check: if it fails, this reports "no detected linking" (0) rather
    // than risk a wrong nonzero value. That keeps the result sound (never
    // falsely reports linking) but not complete (a real link can be
    // missed).
    regina::MatrixInt m(complement.countVertices(), complement.countEdges());
    for (auto e : complement.edges()) {
        m.entry(e->vertex(0)->index(), e->index()) -= 1;
        m.entry(e->vertex(1)->index(), e->index()) += 1;
    }

    regina::MatrixInt n(complement.countEdges(), complement.countTriangles());
    for (auto t : complement.triangles()) {
        for (int j = 0; j < 3; ++j) {
            size_t e = t->edge(j)->index();
            if (t->edgeMapping(j).sign() > 0)
                ++n.entry(e, t->index());
            else
                --n.entry(e, t->index());
        }
    }

    regina::MarkedAbelianGroup h1(m, n);
    if (!h1.isZ())
        return 0;

    // Walk trackedOther's edges into an oriented cycle (the same
    // shared-vertex walk Link::Link() uses to split a multi-component edge
    // set into per-component knots), then read off its class in H_1. Only
    // the magnitude is meaningful, so no canonical orientation needs to be
    // imposed: any consistent walk direction works.
    regina::Vector<regina::Integer> cycle(complement.countEdges());
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
                // branch that actually needs the serialized lookup.
                std::list<regina::CensusHit> hits;
                {
                    std::lock_guard<std::mutex> lock(censusLookupMutex);
                    hits = regina::Census::lookup(complement);
                }
                if (!hits.empty()) {
                    std::cout << "      recognized as "
                              << hits.front().name() << ", "
                              << complement.isoSig() << "\n";
                } else {
                    std::cout << "NOT unknot, " << complement.isoSig()
                              << "\n";
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

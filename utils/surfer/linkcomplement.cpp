//
//  linkcomplement.cpp
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

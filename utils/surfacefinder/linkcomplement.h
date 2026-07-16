#ifndef LINKCOMPLEMENT_H

#define LINKCOMPLEMENT_H

#include <mutex>
#include <string>
#include <unordered_map>
#include <unordered_set>

#include <census/census.h>
#include <triangulation/dim2.h>
#include <triangulation/dim3.h>

// regina::Census::lookup()'s Tokyo Cabinet-backed database is reopened from
// scratch, with no internal locking, on every single call (see
// CensusDB::lookupKey() in census-impl.h) -- so unlike most of Regina's
// engine, it cannot safely be called from more than one thread at a time.
// This holds regardless of whether the census databases have already been
// loaded once: concurrent calls fail unpredictably with "Could not open
// Tokyo Cabinet database ... -> threading error" (observed directly).
// Every regina::Census::lookup() call below serializes through this mutex;
// buildComplement()'s pinch + simplify work is unaffected and stays
// parallel across callers.
inline std::mutex censusLookupMutex;

/** Knot/Link implementation, specialized for taking complements */
class EdgeComplement {
  private:
    const regina::Triangulation<3> *tri_;
    std::vector<const regina::Edge<3> *> edges_;
    std::unordered_map<regina::Tetrahedron<3> *, std::unordered_set<size_t>>
        tetEdges_;

  public:
    EdgeComplement(const regina::Triangulation<3> &tri,
                   const std::vector<const regina::Edge<3> *> &edges)
        : tri_(&tri), edges_(edges) {
        // Add the edge to the correct tetrahedra
        for (const regina::Edge<3> *edge : edges) {
            for (const regina::EdgeEmbedding<3> &emb : edge->embeddings()) {
                tetEdges_[emb.tetrahedron()].insert(emb.face());
            }
        }
    }

    EdgeComplement &operator=(const EdgeComplement &other) {
        if (this != &other) {
            tri_ = other.tri_;
            edges_ = other.edges_;
            tetEdges_ = other.tetEdges_;
        }
        return *this;
    }

    regina::Triangulation<3> buildComplement() const {
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

    bool recognizeComplement() const {
        auto complement = buildComplement();
        std::list<regina::CensusHit> hits;
        {
            std::lock_guard<std::mutex> lock(censusLookupMutex);
            hits = regina::Census::lookup(complement);
        }
        ssize_t genus = complement.recogniseHandlebody();

        if (!hits.empty()) {
            std::cout << "      recognized as " << hits.front().name() << ", "
                      << complement.isoSig() << "\n";
            return true;
        } else if (genus == 1) {
            std::cout << "      unknot, " << complement.isoSig() << "\n";
            return true;
        }
        return false;
    }

    // Non-printing counterpart to recognizeComplement(), for programmatic
    // use: the census name if Regina's census recognizes the complement,
    // "Unknot" if it's a genus-1 handlebody (the complement of an unknotted
    // component), or else the bare isoSig as a fallback identifier.
    std::string identify() const {
        auto complement = buildComplement();
        std::list<regina::CensusHit> hits;
        {
            std::lock_guard<std::mutex> lock(censusLookupMutex);
            hits = regina::Census::lookup(complement);
        }
        if (!hits.empty())
            return hits.front().name();
        if (complement.recogniseHandlebody() == 1)
            return "Unknot";
        return complement.isoSig();
    }

    friend bool operator<(const EdgeComplement &e1, const EdgeComplement &e2) {
        return e1.edges_.size() < e2.edges_.size();
    }

    friend std::ostream &operator<<(std::ostream &os, const EdgeComplement &e) {
        os << "{";
        for (int i = 0; i < e.edges_.size(); ++i) {
            os << e.edges_[i]->index();
            if (i != e.edges_.size() - 1) {
                os << ", ";
            }
        }
        return os << "}";
    }
};

class Knot : public EdgeComplement {
  public:
    Knot(const regina::Triangulation<3> &tri,
         const std::vector<const regina::Edge<3> *> &edges)
        : EdgeComplement(tri, edges) {}

    Knot &operator=(const Knot &other) {
        if (this != &other) {
            EdgeComplement::operator=(other);
        }
        return *this;
    }
};

class Link : public EdgeComplement {
  public:
    std::vector<Knot> comps_;

  public:
    Link(const regina::Triangulation<3> &tri,
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

    Link &operator=(const Link &other) {
        if (this != &other) {
            EdgeComplement::operator=(other);
            comps_ = other.comps_;
        }
        return *this;
    }

    regina::Triangulation<3> buildComplement(int component) const {
        return comps_[component].buildComplement();
    }

    regina::Triangulation<3> buildComplement() const {
        return EdgeComplement::buildComplement();
    }

    int countComponents() const { return comps_.size(); }

    void recognizeComplement() const {
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
                std::list<regina::CensusHit> hits;
                {
                    std::lock_guard<std::mutex> lock(censusLookupMutex);
                    hits = regina::Census::lookup(complement);
                }
                ssize_t genus = complement.recogniseHandlebody();

                std::cout << "    Component " << i + 1 << ": ";
                if (!hits.empty()) {
                    std::cout << "      recognized as " << hits.front().name()
                              << ", " << complement.isoSig() << "\n";
                } else if (genus == 1) {
                    std::cout << "unknot, " << complement.isoSig() << "\n";
                } else if (genus != -1) {
                    std::cout << "\n[!] WARNING! Recognized as a genus "
                              << genus << " handlebody, ";
                    std::cout << "\n[!] This is almost definitely a bug, "
                                 "please "
                                 "report it!\n";
                } else {
                    std::cout << "NOT unknot, " << complement.isoSig() << "\n";
                }
            }
        }
    }

    friend bool operator<(const Link &l1, const Link &l2) {
        if (l1.comps_.size() != l2.comps_.size())
            return l1.comps_.size() < l2.comps_.size();

        return static_cast<EdgeComplement>(l1) <
               static_cast<EdgeComplement>(l2);
    }

    friend std::ostream &operator<<(std::ostream &os, const Link &l) {
        os << "[";
        for (int i = 0; i < l.comps_.size(); ++i) {
            os << l.comps_[i];
            if (i != l.comps_.size() - 1) {
                os << ", ";
            }
        }
        return os << "]";
    }
};

#endif // LINKCOMPLEMENT_H

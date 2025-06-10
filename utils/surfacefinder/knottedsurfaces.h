//
//  knottedsurfaces.h
//
//  Created by John Teague on 06/19/2024.
//

#ifndef KNOTTED_SURFACES_H

#define KNOTTED_SURFACES_H

#include <gmpxx.h>
#include <link/link.h>
#include <triangulation/dim2.h>
#include <triangulation/dim3.h>
#include <triangulation/dim4.h>
#include <unistd.h>

#include <array>
#include <cstddef>
#include <cstring>
#include <iostream>
#include <ostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "census/census.h"
#include "triangulation/forward.h"
#include "triangulation/generic/boundarycomponent.h"
#include "triangulation/generic/face.h"
#include "triangulation/generic/faceembedding.h"
#include "triangulation/generic/triangulation.h"
#include "utilities/exception.h"

#include "gluing.h"

enum class SurfaceCondition { all, boundary, closed };

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
        auto hits = regina::Census::lookup(complement);
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
                auto hits = regina::Census::lookup(complement);
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

/** KnottedSurface Implementation */

template <int dim>
class KnottedSurface {
  private:
    const regina::Triangulation<dim> *tri_;
    /**< The dim-manifold triangulation in which this subdim-manifold is
     * embedded */
    regina::Triangulation<dim - 1> bdry_;
    regina::Triangulation<2> surface_;
    /**< The triangulation of the sub-triangulation induced by the embedding
     */

    std::unordered_map<regina::Triangle<2> *, const regina::Triangle<dim> *>
        emb_;
    /**< A mapping from the top-dimensional simplices of the
     * sub-triangulation into the subdim-simplices of the dim-triangulation
     */
    std::unordered_map<const regina::Triangle<dim> *, regina::Triangle<2> *>
        inv_;
    /**< A mapping from the subdim-simplices of the dim-triangulation to the
     * top-dimensional simplices of the sub-triangulation. Note that
     * inv_.at(emb_.at(face)) = face */

    std::array<int, 3> invariants_;

    std::set<int> indices_;

  public:
    std::unordered_set<const regina::Edge<dim> *> improperEdges_;
    /**< Notes whether an embedding is proper, in the sense that the image
     * of the boundary of the sub-triangulation is entirely contained within
     * the boundary of the dim-manifold triangulation */

    KnottedSurface(const regina::Triangulation<dim> *tri) : tri_(tri) {
        if (!tri_->isClosed()) {
            bdry_ = tri_->boundaryComponent(0)->build();
        }
    }

    KnottedSurface(const KnottedSurface &other)
        : tri_(other.tri_), surface_(other.surface_),
          improperEdges_(other.improperEdges_), invariants_(other.invariants_),
          indices_(other.indices_) {
        if (!tri_->isClosed()) {
            bdry_ = tri_->boundaryComponent(0)->build();
        }
        // Connect the wires...
        for (regina::Triangle<2> *t : other.surface_.triangles()) {
            const regina::Triangle<dim> *f = other.emb_.at(t);
            regina::Triangle<2> *newT = surface_.triangle(t->index());
            emb_[newT] = f;
            inv_[f] = newT;
        }
    }

    KnottedSurface &operator=(const KnottedSurface &other) {
        if (this != &other) {
            tri_ = other.tri_;
            surface_ = other.surface_;
            improperEdges_ = other.improperEdges_;
            invariants_ = other.invariants_;
            indices_ = other.indices_;

            // Connect the wires...
            for (regina::Triangle<2> *t : other.surface_.triangles()) {
                const regina::Triangle<dim> *f = other.emb_.at(t);
                regina::Triangle<2> *newT = surface_.triangle(t->index());
                emb_[newT] = f;
                inv_[f] = newT;
            }
        }

        return *this;
    }

    const regina::Triangulation<2> &surface() const { return surface_; }

    bool isClosed() const { return surface_.isClosed(); }

    bool isProper() const {
        return surface_.isClosed() || improperEdges_.empty();
    }

    bool hasSelfIntersection() const {
        std::unordered_set<const regina::Vertex<dim> *> images;
        for (const regina::Vertex<2> *v : surface_.vertices()) {
            const regina::Vertex<dim> *im = image(v);
            if (images.find(im) != images.end())
                return true;
            images.insert(im);
        }
        return false;
    }

    Link boundary() const {
        std::unordered_set<const regina::Edge<dim - 1> *> edges;

        // TODO: Make this O(1) somehow
        for (const regina::BoundaryComponent<dim> *triBoundaryComp :
             tri_->boundaryComponents()) {

            for (const regina::BoundaryComponent<2> *surfaceBoundaryComp :
                 surface_.boundaryComponents()) {
                for (const regina::Edge<2> *e : surfaceBoundaryComp->edges()) {
                    const regina::Edge<dim> *im = image(e);
                    for (int i = 0; i < triBoundaryComp->countEdges(); ++i) {
                        if (im == triBoundaryComp->edge(i)) {
                            edges.insert(bdry_.edge(i));
                            break;
                        }
                    }
                }
            }
        }

        std::vector<const regina::Edge<dim - 1> *> edgeList(edges.begin(),
                                                            edges.end());
        return {bdry_, edgeList};
    }

    template <int facedim>
    regina::Face<dim, facedim> *image(const regina::Face<2, facedim> *f) const {
        static_assert(0 <= facedim && facedim <= 2,
                      "Must have 0 <= facedim <= 2");

        // for (auto &[t, f] : emb_) {
        //     std::cout << "(" << t->index() << ", " << f->index() << ") ";
        // }
        // std::cout << "\n";

        const regina::FaceEmbedding<2, facedim> &fEmb = f->front();
        return emb_.at(fEmb.simplex())->template face<facedim>(fEmb.face());
    }

    const regina::Triangle<dim> *image(regina::Triangle<2> *t) const {
        return emb_.at(t);
    }

    regina::Triangle<2> *preimage(const regina::Triangle<dim> *f) const {
        auto search = inv_.find(f);
        if (search == inv_.end())
            return nullptr;
        return search->second;
    }

    /** Mutation methods */

    bool addTriangle(const regina::Triangle<dim> *f,
                     const typename GluingNode<dim>::AdjList &adj) {
        // TODO: Figure out a nice way for this method not to know what an
        // AdjList is
        regina::Triangle<2> *src = surface_.newSimplex();
        std::array<bool, 3> isBoundaryEdge = {true, true, true};

        emb_[src] = f;
        inv_[f] = src;

        for (const auto &[adjNode, g] : adj) {
            auto search = inv_.find(adjNode->f);
            if (search == inv_.end())
                continue;

            regina::Triangle<2> *dst = search->second;
            int dstFacet = g.gluing[g.srcFacet];

            // Saves us from catching an error
            if (dst->adjacentSimplex(dstFacet) != nullptr ||
                src->adjacentSimplex(g.srcFacet) != nullptr) {
                surface_.removeSimplex(src);
                emb_.erase(src);
                inv_.erase(f);
                return false;
            }

            isBoundaryEdge[g.srcFacet] = false;
            src->join(g.srcFacet, dst, g.gluing);
        }

        if (hasSelfIntersection()) {
            surface_.removeSimplex(src);
            emb_.erase(src);
            inv_.erase(f);
            return false;
        }

        for (int i = 0; i < 3; ++i) {
            if (!isBoundaryEdge[i]) {
                improperEdges_.erase(f->edge(i));
            } else if (isBoundaryEdge[i] && !f->edge(i)->isBoundary()) {
                improperEdges_.insert(f->edge(i));
            }
        }

        indices_.insert(f->index());
        updateInvariants_();

        return true;
    }

    void removeTriangle(const regina::Triangle<dim> *f) {
        regina::Triangle<2> *src = inv_.at(f);

        for (int i = 0; i < 3; ++i) {
            regina::Triangle<2> *dst = src->adjacentSimplex(i);

            if (dst != nullptr) {
                regina::Edge<dim> *edge = f->edge(i);
                if (!edge->isBoundary()) {
                    improperEdges_.insert(edge);
                }
            }
        }

        surface_.removeSimplex(src);
        emb_.erase(src);
        inv_.erase(f);
        indices_.erase(f->index());
        updateInvariants_();
    }

    std::string detail() const {
        // if (!detail_.empty()) return detail_;

        /* Generate surface details */
        bool isOrientable = invariants_[0];
        int punctures = invariants_[1];
        int genus = invariants_[2];
        std::ostringstream ans;

        if (!surface_.isConnected()) {
            throw regina::InvalidArgument("Surface must be connceted!");
        }

        if (isOrientable) {
            // Special names for surfaces with boundary:
            if (genus == 0 && punctures == 1)
                ans << "Disc";
            else if (genus == 0 && punctures == 2)
                ans << "Annulus";
            else {
                if (genus == 0)
                    ans << "Sphere";
                else if (genus == 1)
                    ans << "Torus";
                else
                    ans << "Orientable genus " << genus << " surface";

                if (punctures == 1)
                    ans << ", 1 puncture";
                else if (punctures > 1)
                    ans << ", " << punctures << " punctures";
            }
        } else {
            // Special names for surfaces with boundary:
            if (genus == 1 && punctures == 1)
                ans << "Möbius band";
            else {
                if (genus == 1)
                    ans << "Projective plane";
                else if (genus == 2)
                    ans << "Klein bottle";
                else
                    ans << "Non-orientable genus " << genus << " surface";

                if (punctures == 1)
                    ans << ", 1 puncture";
                else if (punctures > 1)
                    ans << ", " << punctures << " punctures";
            }
        }

        return ans.str();
    }

    friend bool operator==(const KnottedSurface &lhs,
                           const KnottedSurface &rhs) {
        return lhs.invariants_ == rhs.invariants_ &&
               lhs.indices_ == rhs.indices_;
    }

    friend bool operator<(const KnottedSurface &lhs,
                          const KnottedSurface &rhs) {
        return lhs.invariants_ < rhs.invariants_ ||
               (lhs.invariants_ == rhs.invariants_ &&
                lhs.indices_ < rhs.indices_);
    }

  private:
    void updateInvariants_() {
        bool isOrientable = surface_.isOrientable();
        int punctures = surface_.countBoundaryComponents();
        int genus = isOrientable ? (2 - surface_.eulerChar() - punctures) / 2
                                 : 2 - surface_.eulerChar() - punctures;

        invariants_[0] = isOrientable;
        invariants_[1] = punctures;
        invariants_[2] = genus;
    }
};

#endif

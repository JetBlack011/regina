//
//  surfacefinder.cpp
//
//  Created by John Teague on 06/19/2024.
//

#include <gmpxx.h>
#include <link/link.h>
#include <triangulation/dim2.h>
#include <triangulation/dim3.h>
#include <triangulation/dim4.h>
#include <unistd.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <cstring>
#include <iostream>
#include <map>
#include <ostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include "maths/perm.h"
#include "triangulation/example2.h"
#include "triangulation/forward.h"
#include "triangulation/generic/boundarycomponent.h"
#include "triangulation/generic/face.h"
#include "triangulation/generic/faceembedding.h"
#include "triangulation/generic/triangulation.h"
#include "utilities/exception.h"

/* Constants */
static const regina::Triangulation<2> GENUS_2_SURFACE =
    regina::Triangulation<2>::fromGluings(6, {{0, 2, 5, {2, 1, 0}},
                                              {0, 1, 1, {0, 2, 1}},
                                              {0, 0, 5, {1, 0, 2}},
                                              {1, 1, 2, {0, 2, 1}},
                                              {1, 0, 3, {0, 2, 1}},
                                              {2, 1, 3, {0, 2, 1}},
                                              {2, 0, 4, {0, 2, 1}},
                                              {3, 1, 4, {0, 2, 1}},
                                              {4, 1, 5, {0, 2, 1}}});

// template <int n>
// regina::Triangulation<3> knotComplement(const regina::Triangulation<3> t,
//                                         const Knot &k) {
//     return linkComplement(t, {k});
// }

enum SurfaceCondition { all, boundary, closed };

template <int n>
using Curve = std::vector<const regina::Edge<n> *>;

template <int n>
using DualCurve = std::vector<regina::Triangle<n> *>;

/** Gluing Implementation */

template <int dim, int subdim>
struct Gluing {
    regina::Face<dim, subdim> *src;
    int srcFacet;
    regina::Face<dim, subdim> *dst;
    regina::Perm<subdim + 1> gluing;

    Gluing() = default;

    Gluing(regina::Face<dim, subdim> *src, int srcFacet,
           regina::Face<dim, subdim> *dst, regina::Perm<subdim + 1> gluing)
        : src(src), srcFacet(srcFacet), dst(dst), gluing(gluing) {}

    friend std::ostream &operator<<(std::ostream &os, const Gluing &g) {
        return os << "(" << g.src->index() << ", " << g.srcFacet << ", "
                  << g.dst->index() << ", " << g.gluing << ")";
    }
};

/** Knot implementation */

class Knot {
   private:
    regina::Triangulation<3> tri_;
    std::unordered_map<const regina::Tetrahedron<3> *,
                       std::unordered_set<size_t>>
        tetEdges_;
    size_t numEdges_;

    bool isUnknot_ = false;

    /**
     * Finds a push off of the given knot onto the dual 1-skeleton of the
     * triangulation it's sitting in. Note that this may change the underlying
     * triangulation.
     */
    DualCurve<3> toDualCurve_(const Knot &knot);

   public:
    Knot(regina::Triangulation<3> &tri, const Curve<3> &edges)
        : tri_(tri), numEdges_(edges.size()) {
        for (const regina::Edge<3> *edge : edges) {
            addEdge_(tri_.edge(edge->index()));
        }

        if (edges.size() <= 4) isUnknot_ = true;
    }

    bool isUnknot() const { return isUnknot_; }

    /**
     * Simplifies a given knot so that each tetrahedron it runs through contains
     * at most one of the knot's edges.
     *
     * Algorithm:
     * 0) If we're the unknot (e.g. all edges run through a single tetrahedron)
     * or if every tetrahedron contains exactly one edge of the knot, then we're
     * done.
     * 1) Perform shrinkage. That is, suppose a tetrahedron contains more than
     * one edge. If those edges are connected end to end, we can isotope that
     * edge to go from its start point to its end point.
     * 2) Separate the remaining edges. After step 1), a tetrahedron containing
     * an edge of the knot contains either a single edge or two disconnected
     * edges. In the second case, perform a 1-4 move to separate the edges.
     */
    void simplify() {
        if (isUnknot_) return;
        shrink_();
        separate_();
    }

    friend std::ostream &operator<<(std::ostream &os, const Knot &k) {
        if (k.isUnknot_) {
            return os << "Unknot";
        }

        std::unordered_set<const regina::Edge<3> *> edges;
        for (const auto &[tet, tetEdges] : k.tetEdges_) {
            for (size_t e : tetEdges) {
                edges.insert(tet->edge(e));
            }
        }

        os << "{ ";
        for (const regina::Edge<3> *e : edges) {
            os << "(" << e->vertex(0)->index() << ", " << e->vertex(1)->index()
               << ") ";
        }
        return os << "}";

        //// wtf...
        // if (k.isUnknot_) {
        //     return os << "Unknot";
        // }
        // os << "{";
        // auto it = k.edges_.begin();
        // int next;
        // for (auto it = k.edges_.begin(); it != k.edges_.end() - 1; ++it) {
        //     int v0 = (*it)->vertex(0)->index();
        //     int v1 = (*it)->vertex(1)->index();
        //     int u0 = (*(it + 1))->vertex(0)->index();
        //     int u1 = (*(it + 1))->vertex(1)->index();
        //     if (v0 != u0 && v0 != u1) {
        //         os << v0 << ", ";
        //         next = v1 == u0 ? u0 : u1;
        //     } else {
        //         os << v1 << ", ";
        //         next = v0 == u0 ? u0 : u1;
        //     }
        // }

        // return os << next << "}";
    }

   private:
    bool updateIsUnknot_() {
        if (tetEdges_.size() <= 1) {
            return isUnknot_ = true;
        }

        bool isAllInOneTet = false;

        for (const auto &[_, edges] : tetEdges_) {
            if (edges.size() == numEdges_) {
                return isUnknot_ = true;
            }
        }

        return false;
    }

    void shrink_() {
        bool isDone = false;
        while (updateIsUnknot_() || !isDone) {
            isDone = true;
            for (const auto &[tet, _] : tetEdges_) {
                // If we find a tetrahedron we can shrink, we need to keep going
                if (shrink_(tet)) {
                    isDone = false;
                    // Since at least one tetrahedron's edge information has
                    // been updated, we must start over
                    break;
                }
            }
        }
    }

    void separate_() {
        bool isDone = false;
        while (updateIsUnknot_() || !isDone) {
            // TODO: Obvious optimizations
            isDone = true;
            for (regina::Tetrahedron<3> *tet : tri_.tetrahedra()) {
                const auto &edges = tetEdges_.at(tet);
                assert(edges.size() <= 2);

                if (edges.size() == 2) {
                    fourOneMove_(tet);
                    isDone = false;
                    break;
                }
            }
        }
    }

    bool shrink_(const regina::Tetrahedron<3> *tet) {
        std::unordered_set<size_t> edges = tetEdges_.at(tet);

        std::unordered_map<size_t, size_t> vertCounts;
        for (size_t e : edges) {
            regina::Perm<4> p = tet->edgeMapping(e);
            ++vertCounts[p[0]];
            ++vertCounts[p[1]];
        }

        if (edges.size() < 2 || (edges.size() == 2 && vertCounts.size() == 4))
            return false;

        for (size_t e : edges) {
            removeEdge_(tet->edge(e));
        }

        /* Add the new edge */
        // In this case, there are exactly two endpoints in this
        // tet
        std::array<size_t, 2> endVerts;
        size_t i = 0;
        for (auto [vert, count] : vertCounts) {
            if (count == 1) endVerts[i++] = vert;
            if (i >= 2) break;
        }

        addEdge_(tet->edge(endVerts[0], endVerts[1]));

        return true;
    }

    void addEdge_(const regina::Edge<3> *edge) {
        for (const regina::EdgeEmbedding<3> &emb : edge->embeddings()) {
            tetEdges_[emb.tetrahedron()].insert(emb.face());
        }
    }

    void removeEdge_(const regina::Edge<3> *edge) {
        for (const regina::EdgeEmbedding<3> &emb : edge->embeddings()) {
            const regina::Tetrahedron<3> *tet = emb.tetrahedron();
            tetEdges_[tet].erase(emb.face());
            if (tetEdges_[tet].empty()) {
                tetEdges_.erase(tet);
            }
        }
    }

    std::array<regina::Tetrahedron<3> *, 4> fourOneTets_() {
        auto tets = tri_.newTetrahedra<4>();

        for (int i = 0; i < 3; ++i) {
            for (int j = i + 1; j < 4; ++j) {
                tets[i]->join(j, tets[j], {i, j});
            }
        }

        // Same as
        // tets[0]->join(1, tets[1], {1, 0, 2, 3});
        // tets[0]->join(2, tets[2], {2, 1, 0, 3});
        // tets[0]->join(3, tets[3], {3, 1, 2, 0});
        // tets[1]->join(2, tets[2], {0, 2, 1, 3});
        // tets[1]->join(3, tets[3], {0, 3, 2, 1});
        // tets[2]->join(3, tets[3], {0, 1, 3, 2});

        return tets;
    }

    void fourOneMove_(regina::Tetrahedron<3> *tet) {
        std::array<std::optional<Gluing<3, 3>>, 4> gluings;

        // For each face of tet, track gluings
        for (int i = 0; i < 4; ++i) {
            regina::Tetrahedron<3> *adjTet = tet->adjacentSimplex(i);
            if (adjTet != nullptr) {
                gluings[i] = {tet, i, adjTet, tet->adjacentGluing(i)};
            }
        }

        std::array<regina::Tetrahedron<3> *, 4> tets = fourOneTets_();
        // Add the edges contained in tet to their appropriate position in tets
        for (size_t e : tetEdges_.at(tet)) {
            regina::Perm<4> p = tet->edgeMapping(e);
            tetEdges_[tets[p[2]]].insert(e);
            tetEdges_[tets[p[3]]].insert(e);
        }

        tri_.removeTetrahedron(tet);
        tetEdges_.erase(tet);

        // Glue in the new tetrahedra
        for (int i = 0; i < 4; ++i) {
            if (auto g = gluings[i]) {
                tets[i]->join(i, g->dst, g->gluing);
            }
        }
    }
};

/** Link Implementation */

// std::pair<EdgeMap, regina::Triangulation<3>> drillDualCurve(
//     const regina::Triangulation<3> &tri, const DualCurve<3> &curve) {
//
// }

class Link {
   private:
    regina::Triangulation<3> &tri_;
    std::vector<Knot> comps_;

   public:
    Link(regina::Triangulation<3> &tri, const std::vector<Curve<3>> &components)
        : tri_(tri) {
        for (const Curve<3> &edges : components) {
            comps_.emplace_back(tri, edges);
        }
    }

    bool isUnlink() const {
        for (const Knot &k : comps_) {
            if (!k.isUnknot()) return false;
        }
        return true;
    }

    void simplify() {
        for (Knot &k : comps_) {
            k.simplify();
            break;
        }
    }

    // regina::Triangulation<3> linkComplement(const
    // regina::Triangulation<3> &t,
    //                                         const Link &l) {
    //     regina::Triangulation<3> tri = t;
    //     Link link = l;
    //
    //     while (!link.empty()) {
    //         const Knot &knot = link.front();
    //         // 1) Push off the current link component
    //         DualCurve<3> curve = knotToDualCurve(knot);
    //
    //         /// 2) Drill, keep track of where the remaining link
    //         components end
    //
    //         const auto &[edgeMap, newTri] = drillDualCurve(tri, curve);
    //
    //         Link newLink;
    //         newLink.reserve(link.size() - 1);
    //         for (const Knot &knot2 : link) {
    //             if (knot2 == knot) continue;  // Ignore the now
    //             drilled-out
    //     ot
    //
    //             Knot newKnot;
    //             newKnot.reserve(knot.size());
    //             for (const regina::Edge<3> *edge : knot2) {
    //                 auto search = edgeMap.find(edge->index());
    //                 if (search == edgeMap.end()) {
    //                     throw regina::InvalidArgument(
    //                         "Could not find corresponding edge in
    //                         edgeMap!");
    //                 }
    //                 newKnot.push_back(newTri.edge(search->second));
    //             }
    //             newLink.push_back(newKnot);
    //         }
    //
    //         link = newLink;
    //         tri = newTri;
    //         // Repeat with newLink in the drilled triangulation
    //     }
    //
    //     // Now tri is a triangulation of the complement of the original
    //     link std::cout << tri.detail() << "\n";
    //
    //     return tri;
    // }

    friend std::ostream &operator<<(std::ostream &os, const Link &l) {
        for (const Knot &k : l.comps_) {
            os << k << "\n";
        }
        return os;
    }
};

template <int dim>
class GluingNode {
   public:
    using AdjList = std::unordered_map<GluingNode *, Gluing<dim, 2>>;

    regina::Triangle<dim> *f;
    AdjList adjList;
    bool visited = false;

    GluingNode(regina::Face<dim, 2> *f) : f(f) {}

    friend std::ostream &operator<<(std::ostream &os, const GluingNode &node) {
        os << "(" << node.f->index() << " { ";
        for (const auto &[_, gluing] : node.adjList) {
            os << gluing << " ";
        }

        os << "})";

        return os;
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
        if (!tri_->isClosed()) bdry_ = tri_->boundaryComponent(0)->build();
    }

    KnottedSurface(const KnottedSurface &other)
        : tri_(other.tri_),
          surface_(other.surface_),
          improperEdges_(other.improperEdges_),
          invariants_(other.invariants_),
          indices_(other.indices_) {
        if (!tri_->isClosed()) bdry_ = tri_->boundaryComponent(0)->build();
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

    bool isProper() const {
        return surface_.isClosed() || improperEdges_.empty();
    }

    bool hasSelfIntersection() const {
        std::unordered_set<const regina::Vertex<dim> *> images;
        for (const regina::Vertex<2> *v : surface_.vertices()) {
            const regina::Vertex<dim> *im = image(v);
            if (images.find(im) != images.end()) return true;
            images.insert(im);
        }
        return false;
    }

    Link boundary() {
        const regina::BoundaryComponent<dim> *b = tri_->boundaryComponent(0);
        std::unordered_map<const regina::Edge<dim> *,
                           const regina::Edge<dim - 1> *>
            bdryMap;
        std::unordered_map<const regina::Edge<2> *,
                           const regina::Edge<dim - 1> *>
            edgeMap;

        for (const regina::Edge<dim - 1> *e : bdry_.edges()) {
            bdryMap[b->edge(e->index())] = e;
        }

        for (const regina::BoundaryComponent<2> *comp :
             surface_.boundaryComponents()) {
            for (const regina::Edge<2> *e : comp->edges()) {
                edgeMap[e] = bdryMap.at(image(e));
            }
        }

        std::vector<Curve<dim - 1>> comps;

        // TODO: Likely lots of optimizations can be made here
        for (const regina::BoundaryComponent<2> *comp :
             surface_.boundaryComponents()) {
            const regina::Edge<2> *first = comp->edge(0);
            const regina::Edge<2> *next;
            std::unordered_set<const regina::Vertex<2> *> verts = {
                first->vertex(0), first->vertex(1)};

            Curve<dim - 1> edges = {edgeMap.at(first)};

            while (edges.size() < comp->countEdges()) {
                for (const regina::Edge<2> *e : comp->edges()) {
                    if (std::find(edges.begin(), edges.end(), edgeMap.at(e)) ==
                            edges.end() &&
                        (verts.find(e->vertex(0)) != verts.end() ||
                         verts.find(e->vertex(1)) != verts.end())) {
                        next = e;
                        break;
                    }
                }

                verts = {next->vertex(0), next->vertex(1)};
                edges.push_back(edgeMap.at(next));
            }

            comps.push_back(edges);
        }

        return {bdry_, comps};
    }

    template <int facedim>
    regina::Face<dim, facedim> *image(const regina::Face<2, facedim> *f) const {
        static_assert(0 <= dim <= 2, "Must have 0 <= facedim <= 2");

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
        if (search == inv_.end()) return nullptr;
        return search->second;
    }

    /** Mutation methods */

    bool addTriangle(const regina::Triangle<dim> *f,
                     const typename GluingNode<dim>::AdjList &adj) {
        regina::Triangle<2> *src = surface_.newSimplex();
        std::array<bool, 3> isBoundaryEdge = {true, true, true};

        emb_[src] = f;
        inv_[f] = src;

        for (const auto &[adjNode, g] : adj) {
            auto search = inv_.find(adjNode->f);
            if (search == inv_.end()) continue;

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
                    ans << "genus " << genus << " torus";

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

template <int dim>
class GluingGraph {
   private:
    size_t calls_ = 0;

    SurfaceCondition cond_;

    const regina::Triangulation<dim> &tri_;
    /**< The given triangulation of a dim-manifold */

    std::map<regina::Triangle<dim> *, GluingNode<dim>> nodes_;
    /**< Guaranteed not to be nodes with the same gluing data */

    std::set<KnottedSurface<dim>> surfaces_;

    KnottedSurface<dim> surface_;

   public:
    GluingGraph(const regina::Triangulation<dim> &tri, SurfaceCondition cond)
        : tri_(tri), surface_(&tri), cond_(cond) {
        buildGluingNodes_();
        std::cout << "Built gluing nodes\n";
        buildGluingEdges_();
        std::cout << "Built gluing edges\n";
        // buildGluingInvalids_();
        // std::cout << "Built gluing invalids\n";
    }

    void pruneSurfaces() {
        // TODO: Do this during DFS
        for (auto it = surfaces_.begin(); it != surfaces_.end();) {
            if ((*it).hasSelfIntersection()) {
                it = surfaces_.erase(it);
                continue;
            }
            ++it;
        }
    }

    std::set<KnottedSurface<dim>> &findSurfaces() {
        if (!surfaces_.empty()) return surfaces_;

        for (regina::Triangle<dim> *triangle : tri_.triangles()) {
            reset_();
            dfs_(&nodes_.at(triangle), 0);
            // std::cout << "TOTAL DFS CALLS = " << calls_ << "\n";
        }

        pruneSurfaces();

        std::cout << "TOTAL DFS CALLS = " << calls_ << "\n\n";

        return surfaces_;
    }

   private:
    void nspaces_(int m) {
        for (int i = 0; i < m; ++i) {
            std::cout << "  ";
        }
    }

    void buildGluingNodes_() {
        for (regina::Triangle<dim> *triangle : tri_.triangles()) {
            nodes_.insert({triangle, triangle});
        }
    }

    void buildGluingEdges_() {
        for (auto &[triangle, node] : nodes_) {
            std::unordered_set<int> selfGluingEdges;
            for (int i = 0; i < 3; ++i) {
                regina::Edge<dim> *edge = triangle->edge(i);
                regina::Perm<dim + 1> edgeToTriangle = triangle->edgeMapping(i);

                // for (const regina::TriangleEmbedding<dim> &emb :
                //      triangle->embeddings()) {
                //     // TODO: GAH
                //     const regina::Simplex<dim> *s = emb.simplex();
                //     for (int j = 0; j < 10; ++j) {
                //         regina::Triangle<dim> *other = s->triangle(j);
                //         GluingNode<dim> &otherNode = nodes_.at(other);

                //        for (int k = 0; k < 3; ++k) {
                //            // Only glue identified edges and ignore identity
                //            // mappings/opposite self gluings
                //            if (other->edge(k) != edge ||
                //                (triangle == other &&
                //                 (i == k || selfGluingEdges.find(i) !=
                //                                selfGluingEdges.end())))
                //                continue;

                //            if (triangle == other) {
                //                selfGluingEdges.insert(k);
                //            }

                //            regina::Perm<dim + 1> edgeToOther =
                //                other->edgeMapping(k);
                //            // regina::Perm is immutable, use std::array
                //            instead std::array<int, 3> p; for (int k = 0; k <
                //            3; ++k) {
                //                p[edgeToTriangle[k]] = edgeToOther[k];
                //            }

                //            node.adjList[&otherNode] = {triangle, i, other,
                //            p};
                //        }
                //    }
                //}

                for (auto &[other, otherNode] : nodes_) {
                    for (int j = 0; j < 3; ++j) {
                        // Only glue identified edges and ignore identity
                        // mappings/opposite self gluings
                        if (other->edge(j) != edge ||
                            (triangle == other &&
                             (i == j || selfGluingEdges.find(i) !=
                                            selfGluingEdges.end())))
                            continue;

                        if (triangle == other) {
                            selfGluingEdges.insert(j);
                        }

                        regina::Perm<dim + 1> edgeToOther =
                            other->edgeMapping(j);
                        // regina::Perm is immutable, use std::array instead
                        std::array<int, 3> p;
                        for (int k = 0; k < 3; ++k) {
                            p[edgeToTriangle[k]] = edgeToOther[k];
                        }

                        node.adjList[&otherNode] = {triangle, i, other, p};
                    }
                }
            }
        }
    }

    // void buildGluingInvalids_() {
    //     for (auto &[f, node] : nodes_) {
    //         std::unordered_set<const regina::Triangle<dim> *> adjTriangles;
    //         std::unordered_set<const regina::Vertex<dim> *> vertices;

    //        for (const auto [adj, _] : node.adjList) {
    //            adjTriangles.insert(adj->f);
    //        }

    //        for (int i = 0; i < 3; ++i) {
    //            vertices.insert(f->vertex(i));
    //        }

    //        // TODO: Optimize
    //        for (auto &[adjTriangle, adjNode] : nodes_) {

    //            if (adjTriangles.find(adjTriangle) != adjTriangles.end())
    //                continue;

    //            bool containsVertex = false;
    //            for (int i = 0; i < 3; ++i) {
    //                if (vertices.find(adjTriangle->vertex(i)) !=
    //                vertices.end())
    //                    containsVertex = true;
    //            }

    //            if (containsVertex) node.invalids.push_back(&adjNode);
    //        }
    //    }
    //}

    void dfs_(GluingNode<dim> *node, int layer) {
        ++calls_;
        if (calls_ % 100000 == 0)
            std::cout << "Call " << calls_ << ", Layer " << layer
                      << ", Surfaces " << surfaces_.size() << "\n";

        if (node->visited) {
            return;
        }

        // nspaces_(layer);
        // std::cout << layer << ": ";
        if (!surface_.addTriangle(node->f, node->adjList)) {
            // std::cout << "Failed to add " << *node << "\n";
            return;
        }

        // std::cout << "Added " << *node << "\n";

        node->visited = true;

        if (surface_.isProper()) {
            //  nspaces_(layer);
            //  std::cout << "IMPROPER EDGES: { ";
            //  for (const regina::Edge<dim> *edge : surface_.improperEdges_) {
            //     std::cout << edge->index() << " ";
            // }
            //  std::cout << "}\n";
            surfaces_.insert(surface_);
        }

        for (auto &[nextNode, g] : node->adjList) {
            ////nspaces_(layer);
            // std::cout << layer << ": Gluing " << node->f->index() << " to "
            //           << nextNode->f->index() << " along (" << g.srcFacet
            //           << ", " << g.gluing[g.srcFacet] << ")\n";
            dfs_(nextNode, layer + 1);
        }

        surface_.removeTriangle(node->f);
    }

    void reset_() {
        surface_ = {&tri_};

        for (auto &[_, node] : nodes_) {
            node.visited = false;
        }
    }
};

void usage(const char *progName, const std::string &error = std::string()) {
    if (!error.empty()) std::cerr << error << "\n\n";

    std::cerr << "Usage:\n";
    std::cerr << "    " << progName
              << " { -a, --all | -b, --boundary | -c, --closed | -l, --links } "
                 "<isosig>\n"
                 "    "
              << progName << " [ -v, --version | -h, --help ]\n\n";
    std::cerr << "    -a, --all      : Find all surfaces, regardless of "
                 "boundary conditions\n";
    std::cerr
        << "    -b, --boundary : Find surfaces such that their boundary is "
           "contained entirely in\n"
           "                     the boundary of the given 4-manifold (Note, "
           "if the 4-manifold\n                     is closed, this is "
           "equivalent to --closed)\n";
    std::cerr << "    -c, --closed   : Find only closed surfaces in the given "
                 "4-manifold\n";
    std::cerr << "    -l, --links    : Same as --boundary, but also gives "
                 "a list of link types which\n"
                 "                     are the boundaries of the surfaces "
                 "we find\n\n";
    std::cerr << "    -v, --version  : Show which version of Regina is "
                 "being used\n";
    std::cerr << "    -h, --help     : Display this help\n";
    exit(1);
}

template <int dim>
void surfacesDetail(std::set<KnottedSurface<dim>> &surfaces,
                    SurfaceCondition cond) {
    std::cout << "--- "
              << (cond == SurfaceCondition::all
                      ? ""
                      : (cond == SurfaceCondition::boundary
                             ? "Proper "
                             : (cond == SurfaceCondition::closed ? "Closed "
                                                                 : "")))
              << "Surfaces ---\n";
    int surfaceCount = 0;
    int closedCount = 0;

    std::map<std::string, int> countMap;

    for (const auto &surface : surfaces) {
        ++surfaceCount;
        if (surface.surface().isClosed()) {
            ++closedCount;
        }

        // For checking if properness works
        // bool isProper = true;
        // for (const regina::BoundaryComponent<2> *comp :
        //     surface.surface().boundaryComponents()) {
        //    for (const regina::Edge<2> *edge : comp->edges()) {
        //        if (edge->isBoundary() && !surface.image(edge)->isBoundary())
        //        {
        //            std::cout << "NOTE: " << surface.detail()
        //                      << " is not proper!!\n";
        //            isProper = false;
        //        }
        //    }
        //}

        // if (isProper) {
        //     std::cout << surface.detail() << " is PROPER!"
        //               << "\n";
        // }

        auto search = countMap.find(surface.detail());
        if (search == countMap.end()) {
            countMap.insert({surface.detail(), 1});
        } else {
            ++search->second;
        }
    }

    std::cout << "\n";

    std::vector<std::string> descriptions;
    descriptions.reserve(countMap.size());
    for (const auto &[description, count] : countMap) {
        descriptions.push_back(description);
    }

    std::sort(descriptions.begin(), descriptions.end(),
              [](const std::string &first, const std::string &second) {
                  return first.size() < second.size() ||
                         (first.size() == second.size() && first < second);
              });

    for (std::string &description : descriptions) {
        int count = countMap.find(description)->second;
        std::cout << "- " << count << (count == 1 ? " copy of " : " copies of ")
                  << description << "\n";
    }

    std::cout << "\n";
    if (cond != SurfaceCondition::closed) {
        std::cout << "Total admissible surfaces = " << surfaceCount << "\n";
    }
    std::cout << "Total closed surfaces = " << closedCount << "\n\n";
}

int main(int argc, char *argv[]) {
    // Check for standard arguments:
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") usage(argv[0]);
        if (arg == "-v" || arg == "--version") {
            if (argc != 2)
                usage(argv[0],
                      "Option --version cannot be used with "
                      "any other arguments.");
            std::cout << PACKAGE_BUILD_STRING << "\n";
            exit(0);
        }
    }

    if (argc != 3) {
        usage(argv[0],
              "Please specify a surface condition and provide an "
              "isomorphism signature.");
    }

    std::string arg = argv[1];
    std::string isoSig = argv[2];
    SurfaceCondition cond;

    if (arg == "-a" || arg == "--all") {
        cond = SurfaceCondition::all;
    } else if (arg == "-b" || arg == "--boundary" || arg == "-l" ||
               arg == "--links") {
        cond = SurfaceCondition::boundary;
    } else if (arg == "-c" || arg == "--closed") {
        cond = SurfaceCondition::closed;
    } else {
        usage(argv[0], "Please specify a valid surface condition.");
    }

    //    regina::Triangulation<4> tri;
    //    tri.newSimplex();
    //    tri.pachner(tri.pentachoron(0));
    //
    // regina::Triangulation<3> tri("caba");
    // for (const regina::Tetrahedron<3> *tet : tri.tetrahedra()) {
    //     std::cout << "Tetrahedron " << tet->index() << "\n";
    //     for (int i = 0; i < 4; ++i) {
    //         std::cout << "Vertex " << i << " = " << tet->vertex(i)->index()
    //                   << "\n";
    //     }
    //     std::cout << "\n";
    // }
    // std::vector<int> edgeIndices = {2, 0, 7, 8, 5};
    // Curve<3> c;
    // for (int i : edgeIndices) {
    //     c.push_back(tri.edge(i));
    // }
    // Knot k(tri, c);
    // std::cout << k << "\n";
    // k.simplify();
    // std::cout << k;

    regina::Triangulation<4> tri;
    tri.newSimplex();

    // regina::Triangulation<2> tri = regina::Example<2>::orientable(6, 0);
    // tri.newSimplex();
    // tri.subdivide();
    // tri.subdivide();
    //    tri.subdivide();

    std::cout << tri.detail() << "\n";

    GluingGraph graph(tri, cond);
    auto &surfaces = graph.findSurfaces();

    surfacesDetail(surfaces, cond);

    for (const auto *comp : tri.boundaryComponents()) {
        for (const auto *v : comp->vertices()) {
            std::cout << v->index() << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";

    int numNonUnlinks = 0;
    for (auto surface : surfaces) {
        if (surface.surface().isClosed()) continue;
        Link l = surface.boundary();
        // std::cout << l << "\n";
        l.simplify();
        if (!l.isUnlink()) {
            ++numNonUnlinks;
            std::cout << "DETAIL: " << surface.detail() << "\n";
            std::cout << l << "\n";
        }
    }

    std::cout << "NUMBER OF NONTRIVIAL LINKS = " << numNonUnlinks;

    return 0;
}

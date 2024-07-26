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

#include <array>
#include <cstring>
#include <iostream>
#include <ostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include "maths/perm.h"
#include "triangulation/forward.h"
#include "triangulation/generic/boundarycomponent.h"
#include "triangulation/generic/face.h"
#include "triangulation/generic/faceembedding.h"
#include "triangulation/generic/triangulation.h"

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

enum SurfaceCondition { closed, boundary };

template <int n>
using Curve = std::vector<regina::Edge<n> *>;

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

    template <int, int>
    friend std::ostream &operator<<(std::ostream &os, const Gluing &gluing);
};

template <int dim, int subdim>
std::ostream &operator<<(std::ostream &os, const Gluing<dim, subdim> &g) {
    return os << "(" << g.src->index() << ", " << g.srcFacet << ", "
              << g.dst->index() << ", " << g.gluing << ")";
}

/** Knot implementation */

class Knot {
   private:
    regina::Triangulation<3> &tri_;
    Curve<3> edges_;
    std::unordered_set<const regina::Edge<3> *> edgeSet_;
    std::unordered_map<const regina::Tetrahedron<3> *, int> tetEdgeCount_;

    bool isUnknot_ = false;

    /**
     * Finds a push off of the given knot onto the dual 1-skeleton of the
     * triangulation it's sitting in. Note that this may change the underlying
     * triangulation.
     */
    DualCurve<3> toDualCurve_(const Knot &knot);

   public:
    Knot(regina::Triangulation<3> &tri, const Curve<3> &edges)
        : tri_(tri), edges_(edges) {
        // Precompute the set of edges and the count of edges in each
        // tetrahedron
        for (const regina::Edge<3> *edge : edges) {
            edgeSet_.insert(edge);
            for (const regina::EdgeEmbedding<3> &emb : edge->embeddings()) {
                auto search = tetEdgeCount_.find(emb.tetrahedron());
                if (search == tetEdgeCount_.end()) {
                    tetEdgeCount_.insert({emb.tetrahedron(), 1});
                } else {
                    ++search->second;
                }
            }
        }
    }

    /**
     * Simplifies a given knot so that each tetrahedron it runs through contains
     * at most one of the knot's edges.
     *
     * Algorithm:
     * 0) If we're the unknot (e.g. all edges run through a single tetrahedron)
     * or if every tetrahedron contains exactly one edge of the knot, then we're
     * done.
     * 1) For each tetrahedron containing an edge, if it contains more than one,
     * we have two cases:
     *    i. The edges are connected. In this case, replace these edges with the
     * single edge running from the two endpoints of the path formed in this
     * tetrahedron. Remember to update the tetrahedron edge counts.
     *     ii. There are exactly two edges using all of the vertices of the
     * tetrahedron. In this case, do a 1-4 move. Remember to update the
     * tetrahedron edge counts, including modifying which tetrahedra edges run
     * through.
     * 2) Repeat from 0). It is very possible that we've reduced to only one
     * terahedron, in which case we've arrived at the unknot.
     */
    void simplify() {
        bool isDone = false;
        while (!isUnknot_ && !isDone) {
            for (const auto &[tet, count] : tetEdgeCount_) {
                if (count == 2 && numVerticesUsed_(tet) == 4) {
                    // Separated case

                } else if (count >= 2) {
                    // Connected case
                }
            }
        }
    }

   private:
    int numVerticesUsed_(const regina::Tetrahedron<3> *tet) {
        std::unordered_set<regina::Vertex<3> *> verts;
        // For each edge in tet
        for (int i = 0; i < 6; ++i) {
            auto search = edgeSet_.find(tet->edge(i));
            if (search != edgeSet_.end()) {
                const regina::Edge<3> *e = *search;
                verts.insert(e->vertex(0));
                verts.insert(e->vertex(1));
            }
        }
        return verts.size();
    }

    std::array<regina::Tetrahedron<3> *, 4> fourOneTet() {
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

    void fourOneMove(regina::Tetrahedron<3> *tet) {
        std::array<std::optional<Gluing<3, 3>>, 4> gluings;
        std::vector<regina::Edge<3> *> edgeMap;
        std::unordered_map<regina::Edge<3> *, std::pair<int, int>> edges;

        int numTets = tri_.countTetrahedra();

        // For each triangle in tet
        for (int i = 0; i < 4; ++i) {
            regina::Tetrahedron<3> *adjTet = tet->adjacentSimplex(i);
            if (adjTet != nullptr) {
                gluings[i] = {tet, i, adjTet, tet->adjacentGluing(i)};
            }
        }

        tri_.removeTetrahedron(tet);
        auto tets = fourOneTet();

        for (int i = 0; i < 4; ++i) {
            if (auto g = gluings[i]) {
                tets[i]->join(i, g->dst, g->gluing);
            }
        }

        // Update edge data
        for (int i = 0; i < numTets - 1; ++i) {
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
    std::vector<Knot> components_;

   public:
    Link(regina::Triangulation<3> &tri, const std::vector<Curve<3>> &components)
        : tri_(tri) {
        for (const Curve<3> &edges : components) {
            components_.emplace_back(tri, edges);
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
};

/** TriangulationEmbedding Implementation */

template <int dim, int subdim>
class TriangulationEmbedding {
   public:
    using EmbeddingMap = std::unordered_map<regina::Simplex<subdim> *,
                                            const regina::Face<dim, subdim> *>;
    using InverseMap = std::unordered_map<const regina::Face<dim, subdim> *,
                                          regina::Simplex<subdim> *>;

    const regina::Triangulation<dim> *tri_;
    /**< The dim-manifold triangulation in which this subdim-manifold is
     * embedded */
    regina::Triangulation<subdim> sub_;
    /**< The triangulation of the sub-triangulation induced by the embedding
     */

    EmbeddingMap emb_;
    /**< A mapping from the top-dimensional simplices of the
     * sub-triangulation into the subdim-simplices of the dim-triangulation
     */
    InverseMap inv_;
    /**< A mapping from the subdim-simplices of the dim-triangulation to the
     * top-dimensional simplices of the sub-triangulation. Note that
     * inv_.at(emb_.at(face)) = face */

    bool isProper_ = true;
    /**< Notes whether an embedding is proper, in the sense that the image
     * of the boundary of the sub-triangulation is entirely contained within
     * the boundary of the dim-manifold triangulation */

    // std::vector<int> cmp_;
    /**< A precomputed list used to compare two embeddings for insertion
     * into e.g. a set. In general, this will consist of some easy
     * invariants like genus for subdim = 2 followed by the list of
     * subdim-simplices appearing in the embedding of the sub-triangulation
     * into the dim-triangulation. */

    TriangulationEmbedding(const regina::Triangulation<dim> *tri) : tri_(tri) {
        static_assert(2 <= subdim <= dim, "Must have 2 <= subdim <= dim");
    }

    TriangulationEmbedding<dim, subdim> &operator=(
        const TriangulationEmbedding<dim, subdim> &src) = default;

    bool isProper() const { return isProper_; }

    template <int facedim>
    regina::Face<dim, facedim> *image(
        const regina::Face<subdim, facedim> *f) const {
        static_assert(0 <= facedim <= subdim,
                      "Must have 0 <= facedim <= subdim");

        const regina::FaceEmbedding<subdim, facedim> &fEmb = f->front();
        return emb_.at(fEmb.simplex())->template face<facedim>(fEmb.face());
    }

    regina::Face<dim, subdim> *image(const regina::Simplex<subdim> *s) const {
        return emb_.at(s);
    }

    regina::Simplex<subdim> *preimage(
        const regina::Face<dim, subdim> *s) const {
        auto search = inv_.find(s);
        if (search == inv_.end()) return nullptr;
        return search->second;
    }

    regina::Simplex<subdim> *addFace(const regina::Face<dim, subdim> *f) {
        auto search = inv_.find(f);

        if (search != inv_.end()) {
            return search->second;
        }

        regina::Simplex<subdim> *s = sub_.newSimplex();
        emb_[s] = f;
        inv_[f] = s;

        return s;
    }

    /** Mutation methods */

    bool addGluing(Gluing<dim, subdim> g) {
        regina::Simplex<subdim> *src = inv_.at(g.src);
        regina::Simplex<subdim> *dst = addFace(g.dst);

        try {
            src->join(g.srcFacet, dst, g.gluing);
        } catch (regina::InvalidArgument &e) {
            if (!sub_.isConnected()) {
                sub_.removeSimplex(dst);
                emb_.erase(dst);
                inv_.erase(g.dst);
            }
            return false;
        }

        if (isProper_) {
            // Quick optimization, boundary containment is entirely
            // determined by the new facets we've made
            for (int i = 0; i < subdim + 1; ++i) {
                if (i == g.gluing[g.srcFacet]) continue;

                const regina::Face<subdim, subdim - 1> *facet =
                    dst->template face<subdim - 1>(i);
                const regina::Face<dim, subdim - 1> *im = image(facet);
                if (facet->isBoundary() && !im->isBoundary()) {
                    isProper_ = false;
                }
            }
        } else {
            updateIsProper_();
        }

        return true;
    }

    bool removeGluing(Gluing<dim, subdim> g) {
        regina::Simplex<subdim> *src = inv_.at(g.src);

        if (src->unjoin(g.srcFacet) == nullptr) {
            return false;
        }

        // If removing this gluing resulted in an isolated simplex
        if (!sub_.isConnected()) {
            regina::Simplex<subdim> *dst = inv_.at(g.dst);
            sub_.removeSimplex(dst);
            emb_.erase(dst);
            inv_.erase(g.dst);
        }

        return true;
    }

    /**
     * Two embeddings are considered equal if they have isomorphic
     * sub-triangulations and have the same embedding data.
     */
    // template <int, int>
    // friend bool operator==(const TriangulationEmbedding &lhs,
    //                        const TriangulationEmbedding &rhs);

    /**
     * Arbitrary tie breaker for sorting in a set by lexicographic ordering
     * on some integer invariants and the sorted set of embedded face
     * indices.
     */
    // template <int, int>
    // friend bool operator<(const TriangulationEmbedding &lhs,
    //                       const TriangulationEmbedding &rhs);

    template <int, int>
    friend std::ostream &operator<<(
        std::ostream &os, const TriangulationEmbedding &TriangulationEmbedding);

   private:
    void updateIsProper_() {
        for (const regina::BoundaryComponent<subdim> *comp :
             sub_.boundaryComponents()) {
            for (const regina::Face<subdim, subdim - 1> *f : comp->facets()) {
                regina::Face<dim, subdim - 1> *im = image(f);
                if (!im->isBoundary()) {
                    isProper_ = false;
                    return;
                }
            }
        }

        isProper_ = true;
    }
};

// template <int dim, int subdim>
// bool operator==(const TriangulationEmbedding<dim, subdim> &lhs,
//                 const TriangulationEmbedding<dim, subdim> &rhs) {
//     return lhs.cmp_ == rhs.cmp_;
// }
//
// template <int dim, int subdim>
// bool operator<(const TriangulationEmbedding<dim, subdim> &lhs,
//                const TriangulationEmbedding<dim, subdim> &rhs) {
//     return lhs.cmp_ < rhs.cmp_;
// }

template <int dim, int subdim>
std::ostream &operator<<(std::ostream &os,
                         const TriangulationEmbedding<dim, subdim> &emb) {
    os << "{ ";
    for (const auto [key, val] : emb.emb_) {
        os << "(" << key->index() << "," << val->index() << ") ";
    }

    return os << "}";
}

/** SurfaceEmbedding Implementation */

template <int dim>
class SurfaceEmbedding : public TriangulationEmbedding<dim, 2> {
   private:
    // std::string detail_;
    /**< A textual description of the underlying surface, e.g.
     * "Non-orientable genus 5 surface, 2 punctures" */
    std::vector<int> invariants_;
    std::set<int> faceIndices_;

   public:
    SurfaceEmbedding(const regina::Triangulation<dim> *tri)
        : TriangulationEmbedding<dim, 2>(tri) {
        invariants_ = {0, 0, 0};
        updateInvariants_();
    }

    const regina::Triangulation<2> &surface() const {
        return TriangulationEmbedding<dim, 2>::sub_;
    }

    std::string detail() const {
        // if (!detail_.empty()) return detail_;

        /* Generate surface details */
        bool isOrientable = invariants_[0];
        int punctures = invariants_[1];
        int genus = invariants_[2];
        std::ostringstream ans;

        if (!TriangulationEmbedding<dim, 2>::sub_.isConnected()) {
            throw regina::InvalidArgument("Surface must be connceted!");
        }

        if (isOrientable) {
            // Special names for surface_s with boundary:
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
            // Special names for surface_s with boundary:
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
        // return detail_ = ans.str();
    }

    regina::Triangle<2> *addFace(const regina::Triangle<dim> *f) {
        int numSimplices =
            TriangulationEmbedding<dim, 2>::sub_.countTriangles();
        regina::Triangle<2> *t = TriangulationEmbedding<dim, 2>::addFace(f);
        if (numSimplices <
            TriangulationEmbedding<dim, 2>::sub_.countTriangles()) {
            faceIndices_.insert(f->index());
        }

        return t;
    }

    bool addGluing(Gluing<dim, 2> g) {
        int numSimplices =
            TriangulationEmbedding<dim, 2>::sub_.countTriangles();
        bool success = TriangulationEmbedding<dim, 2>::addGluing(g);

        if (success) {
            updateInvariants_();
            if (numSimplices <
                TriangulationEmbedding<dim, 2>::sub_.countTriangles())
                faceIndices_.insert(g.dst->index());
        }

        return success;
    }

    bool removeGluing(Gluing<dim, 2> g) {
        int numSimplices =
            TriangulationEmbedding<dim, 2>::sub_.countTriangles();
        bool success = TriangulationEmbedding<dim, 2>::removeGluing(g);

        if (success) {
            updateInvariants_();
            if (TriangulationEmbedding<dim, 2>::sub_.countTriangles() <
                numSimplices)
                faceIndices_.erase(g.dst->index());
        }

        return success;
    }

    friend bool operator==(const SurfaceEmbedding &lhs,
                           const SurfaceEmbedding &rhs) {
        return lhs.invariants_ == rhs.invariants_ &&
               lhs.faceIndices_ == rhs.faceIndices_;
    }

    friend bool operator<(const SurfaceEmbedding &lhs,
                          const SurfaceEmbedding &rhs) {
        return lhs.invariants_ < rhs.invariants_ ||
               (lhs.invariants_ == rhs.invariants_ &&
                lhs.faceIndices_ < rhs.faceIndices_);
    }

   private:
    void updateInvariants_() {
        bool isOrientable = TriangulationEmbedding<dim, 2>::sub_.isOrientable();
        int punctures =
            TriangulationEmbedding<dim, 2>::sub_.countBoundaryComponents();
        int genus =
            isOrientable
                ? (2 - TriangulationEmbedding<dim, 2>::sub_.eulerChar() -
                   punctures) /
                      2
                : 2 - TriangulationEmbedding<dim, 2>::sub_.eulerChar() -
                      punctures;

        invariants_[0] = isOrientable;
        invariants_[1] = punctures;
        invariants_[2] = genus;
    }
};

using KnottedSurface = SurfaceEmbedding<4>;

template <int dim>
class GluingGraph {
   public:
    class GluingNode {
       public:
        const Gluing<dim, 2> gluing_;
        std::unordered_set<GluingNode *> adjList_;
        std::unordered_set<GluingNode *> invalids_;  // Make sure to precompute
        bool valid_ = true;
        bool visited_ = false;

        GluingNode(regina::Triangle<dim> *src, int srcFacet,
                   regina::Triangle<dim> *dst, regina::Perm<3> gluing)
            : gluing_(src, srcFacet, dst, gluing) {}

        /**
         * Comparison operators let us guarantee gluing uniqueness using
         * std::set, since I really don't want to come up with a hash
         * function.
         */
        friend bool operator==(const GluingNode &lhs, const GluingNode &rhs) {
            // We never want two gluing nodes with the same gluing, use this
            // to distinguish the nodes
            return lhs.gluing_.src == rhs.gluing_.src &&
                   lhs.gluing_.srcFacet == rhs.gluing_.srcFacet &&
                   lhs.gluing_.dst == rhs.gluing_.dst &&
                   lhs.gluing_.gluing == rhs.gluing_.gluing;
        }

        friend bool operator<(const GluingNode &lhs, const GluingNode &rhs) {
            // Arbitrary tie breaker. If by some miracle you manage to
            // stick in a triangulation with over LONG_MAX triangles, God
            // help you.
            std::vector<long> lhsLex = {
                lhs.gluing_.src == nullptr
                    ? -1
                    : static_cast<long>(lhs.gluing_.src->index()),
                static_cast<long>(lhs.gluing_.dst->index()),
                lhs.gluing_.srcFacet,
                lhs.gluing_.gluing[0],
                lhs.gluing_.gluing[1],
                lhs.gluing_.gluing[2]};
            std::vector<long> rhsLex = {
                rhs.gluing_.src == nullptr
                    ? -1
                    : static_cast<long>(rhs.gluing_.src->index()),
                static_cast<long>(rhs.gluing_.dst->index()),
                rhs.gluing_.srcFacet,
                rhs.gluing_.gluing[0],
                rhs.gluing_.gluing[1],
                rhs.gluing_.gluing[2]};

            return lhsLex < rhsLex;
        }

        friend std::ostream &operator<<(std::ostream &os,
                                        const GluingNode &node) {
            os << "(" << node.gluing_ << " { ";
            for (const typename GluingGraph<dim>::GluingNode *adj :
                 node.adjList_) {
                os << adj->gluing_ << " ";
            }

            os << "})";

            return os;
        }
    };

   private:
    size_t calls_ = 0;

    SurfaceCondition cond_;

    const regina::Triangulation<dim> *tri_;
    /**< The given triangulation of an n-manifold */

    std::vector<GluingNode> nodes_;
    /**< Guaranteed not to be nodes with the same gluing data */

    SurfaceEmbedding<dim> surface_;

   public:
    GluingGraph(regina::Triangulation<dim> *tri, SurfaceCondition cond)
        : tri_(tri), surface_(tri), cond_(cond) {
        buildGluingNodes_(tri);
        buildGluingEdges_();
        // buildGluingInvalids_();
    }

    std::set<SurfaceEmbedding<dim>> findSurfaces() {
        std::set<SurfaceEmbedding<dim>> ans;

        for (const regina::Triangle<dim> *triangle : tri_->triangles()) {
            auto surfaces = findSurfaces_(triangle);
            for (const SurfaceEmbedding<dim> &surface : surfaces) {
                ans.insert(surface);
            }
            break;
        }

        std::cout << "TOTAL DFS CALLS = " << calls_ << "\n";

        return ans;
    }

   private:
    void nspaces_(int m) {
        for (int i = 0; i < m; ++i) {
            std::cout << "  ";
        }
    }

    std::set<GluingNode> addTriangle_(regina::Triangle<dim> *triangle) {
        std::set<GluingNode> nodes;

        // For each edge of our triangle, find all of the ways adjacent
        // triangles are glued to that edge in the given triangulation
        for (int facet = 0; facet < 3; ++facet) {
            const regina::Edge<dim> *edge = triangle->edge(facet);
            regina::Perm<dim + 1> edgeToTriangle = triangle->edgeMapping(facet);

            for (regina::Triangle<dim> *other : tri_->triangles()) {
                for (int i = 0; i < 3; ++i) {
                    // Only glue identified edges and ignore identity
                    // mappings
                    if (other->edge(i) != edge ||
                        (other->edge(i) == edge && triangle == other))
                        continue;

                    regina::Perm<dim + 1> edgeToOther = other->edgeMapping(i);
                    // regina::Perm is immutable, use std::array instead
                    std::array<int, 3> p;

                    for (int i = 0; i < 3; ++i) {
                        p[edgeToTriangle[i]] = edgeToOther[i];
                    }

                    GluingNode node = {
                        triangle, facet, other, {p[0], p[1], p[2]}};
                    nodes.insert(node);
                }
            }
        }

        return nodes;
    }

    void buildGluingNodes_(regina::Triangulation<dim> *tri) {
        std::set<GluingNode> nodes;

        // TODO: Can we do better than O(|triangles|^2)?
        for (regina::Triangle<dim> *triangle : tri->triangles()) {
            std::set<GluingNode> newNodes = addTriangle_(triangle);

            for (const GluingNode &node : newNodes) {
                nodes.insert(node);
            }
        }

        for (const GluingNode &node : nodes) {
            nodes_.push_back(node);
        }
    }

    void buildGluingEdges_() {
        for (GluingNode &node : nodes_) {
            for (GluingNode &other : nodes_) {
                if (node.gluing_.dst == other.gluing_.src) {
                    node.adjList_.insert(&other);
                }
            }
        }
    }

    void dfs_(GluingNode *node, std::set<SurfaceEmbedding<dim>> &surfaces,
              int layer) {
        ++calls_;

        if (calls_ % 10000 == 0)
            std::cout << "Call " << calls_ << "\n";

        if (!node->valid_) {
            return;
        }

        // std::cout << "GLUING " << node->gluing_ << "\n";

        if (!surface_.addGluing(node->gluing_)) {
            return;
        }

        for (GluingNode *nextNode : node->adjList_) {
            if (nextNode->gluing_.dst == node->gluing_.src) continue;
            regina::Triangle<2> *t = surface_.preimage(nextNode->gluing_.dst);
            if (t != nullptr) {
                if (surface_.addGluing(nextNode->gluing_))
                    nextNode->visited_ = true;
            }
        }

        if (tri_->isClosed() && surface_.surface().isClosed() || surface_.isProper()) {
            surfaces.insert(surface_);
            std::cout << layer << ": " << surface_ << "\n"
                      << surface_.detail() << "\n"
                      << surface_.surface().detail() << "\n";
        }

        // Recursive step: walk along each possible next gluing
        for (GluingNode *nextNode : node->adjList_) {
            if (nextNode->gluing_.dst != node->gluing_.src &&
                !nextNode->visited_)
                dfs_(nextNode, surfaces, layer + 1);
        }

        // Clean up
        for (GluingNode *nextNode : node->adjList_) {
            if (nextNode->gluing_.dst != node->gluing_.src &&
                nextNode->visited_) {
                surface_.removeGluing(nextNode->gluing_);
                nextNode->visited_ = false;
            }
        }

        surface_.removeGluing(node->gluing_);
    }

    void reset_() {
        surface_ = {tri_};

        for (GluingNode &node : nodes_) {
            node.valid_ = true;
            node.visited_ = false;
        }
    }

    std::set<SurfaceEmbedding<dim>> findSurfaces_(
        const regina::Triangle<dim> *triangle) {
        std::set<SurfaceEmbedding<dim>> ans;

        reset_();

        surface_.addFace(triangle);

        for (GluingNode &node : nodes_) {
            if (node.gluing_.src == triangle) {
                dfs_(&node, ans, 1);
            }
        }

        std::cout << "TOTAL DFS CALLS = " << calls_ << "\n";

        return ans;
    }
};

void usage(const char *progName, const std::string &error = std::string()) {
    if (!error.empty()) std::cerr << error << "\n\n";

    std::cerr << "Usage:\n";
    std::cerr << "    " << progName
              << " [ -c, --closed | -b, --boundary | -l, --links ] <isosig>\n"
                 "    "
              << progName << " [ -v, --version | -h, --help ]\n\n";
    std::cerr << "    -c, --closed   : Find closed surfaces in the given "
                 "4-manifold\n";
    std::cerr << "    -b, --boundary : Find surfaces such that its boundary is "
                 "contained entirely in\n"
                 "                     the boundary of the given 4-manifold\n";
    std::cerr << "    -l, --links    : Same as --boundary, but also gives "
                 "a list of "
                 "link types which\n"
                 "                     are the boundaries of the surfaces "
                 "we find\n";
    std::cerr << "    -v, --version  : Show which version of Regina is "
                 "being used\n";
    std::cerr << "    -?, --help     : Display this help\n";
    exit(1);
}

int main(int argc, char *argv[]) {
    // Check for standard arguments:
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0)
            usage(argv[0]);
        if (strcmp(argv[i], "-v") == 0 || strcmp(argv[i], "--version") == 0) {
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

    std::string isoSig = argv[2];
    SurfaceCondition cond;

    if (strcmp(argv[1], "-c") == 0 || strcmp(argv[1], "--closed") == 0) {
        cond = SurfaceCondition::closed;
    } else if (strcmp(argv[1], "-b") == 0 ||
               strcmp(argv[1], "--boundary") == 0) {
        cond = SurfaceCondition::boundary;
    } else if (strcmp(argv[1], "-l") == 0 || strcmp(argv[1], "--links") == 0) {
        cond = SurfaceCondition::boundary;
    } else {
        usage(argv[0], "Please specify a valid surface condition.");
    }

    regina::Triangulation<3> tri(isoSig);
    //tri.subdivide();

    // regina::Triangulation<3> tri;
    // tri.newSimplex();
    // tri.newSimplex();
    // tri.newSimplex();
    // tri.newSimplex();

    // tri.tetrahedron(0)->join(0, tri.tetrahedron(1), {0, 1, 2, 3});
    // tri.tetrahedron(1)->join(1, tri.tetrahedron(2), {0, 1, 2, 3});

    // for (int i = 0; i < tri.countTetrahedra(); ++i) {
    //     std::cout << "Tetrahedron " << i << ": " << tri.tetrahedron(i) <<
    //     "\n"; for (int j = 0; j < 6; ++j) {
    //         auto mapping = tri.tetrahedron(i)->edgeMapping(j);
    //         std::cout << "Edge " << tri.tetrahedron(i)->edge(j)->index() << "
    //         ("
    //                   << mapping[0] << ", " << mapping[1]
    //                   << "): " << tri.tetrahedron(i)->edge(j) << "\n";
    //     }
    //     std::cout << "\n";
    // }

    // tri.removeTetrahedron(tri.tetrahedron(2));
    ////tri.newTetrahedron()->join(3, tri.tetrahedron(0), {0, 1, 2, 3});

    // for (int i = 0; i < tri.countTetrahedra(); ++i) {
    //     std::cout << "Tetrahedron " << i << ": " << tri.tetrahedron(i) <<
    //     "\n"; for (int j = 0; j < 6; ++j) {
    //         auto mapping = tri.tetrahedron(i)->edgeMapping(j);
    //         std::cout << "Edge " << tri.tetrahedron(i)->edge(j)->index() << "
    //         ("
    //                   << mapping[0] << ", " << mapping[1]
    //                   << "): " << tri.tetrahedron(i)->edge(j) << "\n";
    //     }
    //     std::cout << "\n";
    // }

    // tri.subdivide();
    ////  auto tri = regina::Example<4>::fourSphere();
    //// regina::Triangulation<4> tri(isoSig);
    //// regina::Triangulation<3> tri("caba");
    std::cout << tri.detail() << "\n\n";
    // Knot<3> k;

    // knotComplement(tri, k);

    std::cout << "Boundary = ";
    for (const auto comp : tri.boundaryComponents()) {
        for (const auto edge : comp->edges()) {
            std::cout << edge->index() << " ";
        }
    }
    std::cout << "\n\n";

    std::cout << "--- Closed Surfaces ---\n";
    int surfaceCount = 0;
    int closedCount = 0;
    GluingGraph graph(&tri, cond);

    auto surfaces = graph.findSurfaces();
    std::map<std::string, int> countMap;

    for (const auto &surface : surfaces) {
        ++surfaceCount;
        if (surface.surface().isClosed()) {
            ++closedCount;
        }

        // std::cout << " - " << surface << "\n";
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
    std::cout << "Total admissible surfaces = " << surfaceCount
              << "\nTotal closed surfaces = " << closedCount << "\n\n";

    return 0;
}

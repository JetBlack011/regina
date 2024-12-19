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

#include <iostream>
#include <ostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "triangulation/forward.h"
#include "triangulation/generic/face.h"
#include "triangulation/generic/triangulation.h"

/* Custom types */
enum SurfaceCondition { closed, boundary };

template <int n>
using Curve = std::vector<regina::Edge<n> *>;

template <int n>
using DualCurve = std::vector<regina::Triangle<n> *>;

class Knot {
   private:
    using EdgeSet = std::unordered_set<const regina::Edge<3> *>;
    using TetEdgeCount =
        std::unordered_map<const regina::Tetrahedron<3> *, int>;

    const regina::Triangulation<3> &tri_;
    Curve<3> edges_;
    EdgeSet edgeSet_;
    TetEdgeCount tetEdgeCount_;

    bool isUnknot_ = false;

    /**
     * Finds a push off of the given knot onto the dual 1-skeleton of the
     * triangulation it's sitting in. Note that this may change the underlying
     * triangulation.
     */
    DualCurve<3> toDualCurve_(const Knot &knot);

   public:
    Knot(const regina::Triangulation<3> &tri, const Curve<3> &edges);

    /**
     * Simplifies the knot so that every tetrahedron in the triangulation
     * contains at most one edge of the knot.
     *
     * Outline:
     * 1) Base case: if every edge is contained in one tetrahedron, this is the
     * unknot.
     * 2) For each tetrahedron we know contains an edge, if it contains more
     * than one edge, there are two cases: either these edges are a connected
     * path or there are exactly two edges touching every vertex in the
     * tetrahedron. In the first case, replace all of these edges with the edge
     * going from the incoming vertex to the outgoing vertex, and remember to
     * update the counts of other tetrahedra containing the deleted edges. In
     * the second, perform a 1-4 move, and update the counts.
     * 3) Repeat 1)
     */
    void simplify();

    regina::Triangulation<3> complement();
};

class Link {
   private:
    const regina::Triangulation<3> tri_;
    std::vector<Knot> components_;

   public:
    Link(const regina::Triangulation<3> &tri,
         const std::vector<Curve<3>> &components);

    /**
     * Identifies a link in a given 3-manifold triangulation.
     *
     * Outline:
     * For each component of the link,
     * 1) Push off the given link into the dual 1-skeleton
     * 2) Drill the dual curve to obtain the link complement, keeping track of
     * where the other components of the link end up in the complement.
     * 3) Identify the link using the complement. SnapPy has a routine for this,
     * might end up running Python code in here.
     */
    regina::Triangulation<3> complement(const regina::Triangulation<3> &t,
                                        const Link &l);
};

template <int dim>
struct Gluing {
    const regina::Simplex<dim> *src;
    const int srcFacet;
    const regina::Simplex<dim> *dst;
    const regina::Perm<dim + 1> gluing;

    template <int>
    friend std::ostream &operator<<(std::ostream &os, const Gluing &gluing);
};

/**
 * Represents a surface with a piecewise-linear embedding into a triangulated
 * dim-manifold.
 */
template <int dim, int subdim>
class TriangulationEmbedding {
   public:
    using EmbeddingMap = std::unordered_map<regina::Simplex<subdim> *,
                                            regina::Face<dim, subdim> *>;
    using InverseMap = std::unordered_map<regina::Face<dim, subdim> *,
                                          regina::Simplex<subdim> *>;

   protected:
    const regina::Triangulation<dim> &tri_;
    /**< The dim-manifold triangulation in which this subdim-manifold is
     * embedded */
    regina::Triangulation<subdim> sub_;
    /**< The triangulation of the sub-triangulation induced by the embedding */

    EmbeddingMap emb_;
    /**< A mapping from the top-dimensional simplices of the
     * sub-triangulation into the subdim-simplices of the dim-triangulation */
    InverseMap inv_;
    /**< A mapping from the subdim-simplices of the dim-triangulation to the
     * top-dimensional simplices of the sub-triangulation. Note that
     * inv_.at(emb_.at(face)) = face */

    bool isProper_ = true;
    /**< Notes whether an embedding is proper, in the sense that the image of
     * the boundary of the sub-triangulation is entirely contained within the
     * boundary of the dim-manifold triangulation */

    std::vector<int> cmp_;
    /**< A precomputed list used to compare two embeddings for insertion into
     * e.g. a set. In general, this will consist of some easy invariants like
     * genus for subdim = 2 followed by the list of subdim-simplices appearing
     * in the embedding of the sub-triangulation into the dim-triangulation. */

   public:
    TriangulationEmbedding(const regina::Triangulation<dim> &tri);

    /** Const methods */

    /**
     * Checks whether the boundary of the embedded surface is contained entirely
     * within the boundary of the triangulation it's embedded in.
     */
    bool isProper() const;

    void updateIsProper();

    template <int facedim>
    regina::Face<dim, facedim> *image(
        const regina::Face<subdim, facedim> *f) const;

    regina::Face<dim, subdim> *image(const regina::Simplex<subdim> *s) const;

    regina::Face<subdim, subdim> *preimage(
        const regina::Face<dim, subdim> *f) const;

    const regina::Triangulation<2> &subtriangulation() const;

    const EmbeddingMap &embedding() const;

    /** Mutation methods for use when building an embedded triangulation */

    regina::Simplex<subdim> *addFace(const regina::Face<dim, subdim> *f);

    bool addGluing(Gluing<dim> g);

    void removeGluing(Gluing<dim> g);

    /**
     * Two embeddings are considered equal if they have isomorphic
     * sub-triangulations and have the same embedding data.
     */
    template <int, int>
    friend bool operator==(const TriangulationEmbedding &lhs,
                           const TriangulationEmbedding &rhs);

    /**
     * Arbitrary tie breaker for sorting in a set by lexicographic ordering on
     * some integer invariants and the sorted set of embedded face indices.
     */
    template <int, int>
    friend bool operator<(const TriangulationEmbedding &lhs,
                          const TriangulationEmbedding &rhs);

    template <int, int>
    friend std::ostream &operator<<(
        std::ostream &os, const TriangulationEmbedding &TriangulationEmbedding);
};

class KnottedSurface : public TriangulationEmbedding<4, 2> {
   private:
    //std::string detail_;
    /**< A textual description of the underlying surface, e.g. "Non-orientable
     * genus 5 surface, 2 punctures" */

    std::vector<int> invariants_;

   public:
    KnottedSurface(const regina::Triangulation<4> &tri);

    const regina::Triangulation<2> &surface() const;

    std::string detail() const;

private:
    void computeInvariants_();
};

#endif


/**************************************************************************
 *                                                                        *
 *  Regina - A Normal Surface Theory Calculator                           *
 *  Computational Engine                                                  *
 *                                                                        *
 *  Copyright (c) 1999-2014, Ben Burton                                   *
 *  For further details contact Ben Burton (bab@debian.org).              *
 *                                                                        *
 *  This program is free software; you can redistribute it and/or         *
 *  modify it under the terms of the GNU General Public License as        *
 *  published by the Free Software Foundation; either version 2 of the    *
 *  License, or (at your option) any later version.                       *
 *                                                                        *
 *  As an exception, when this program is distributed through (i) the     *
 *  App Store by Apple Inc.; (ii) the Mac App Store by Apple Inc.; or     *
 *  (iii) Google Play by Google Inc., then that store may impose any      *
 *  digital rights management, device limits and/or redistribution        *
 *  restrictions that are required by its terms of service.               *
 *                                                                        *
 *  This program is distributed in the hope that it will be useful, but   *
 *  WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU     *
 *  General Public License for more details.                              *
 *                                                                        *
 *  You should have received a copy of the GNU General Public             *
 *  License along with this program; if not, write to the Free            *
 *  Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston,       *
 *  MA 02110-1301, USA.                                                   *
 *                                                                        *
 **************************************************************************/

/* end stub */

/*! \file dim2/dim2edge.h
 *  \brief Deals with edges in the 1-skeleton of a 2-manifold triangulation.
 */

#ifndef __DIM2EDGE_H
#ifndef __DOXYGEN
#define __DIM2EDGE_H
#endif

#include "regina-core.h"
#include "output.h"
#include "generic/face.h"
#include "maths/nperm3.h"
// NOTE: More #includes follow after the class declarations.

namespace regina {

class Dim2BoundaryComponent;

template <int> class Component;
template <int> class Simplex;
template <int> class Triangulation;
typedef Component<2> Dim2Component;
typedef Simplex<2> Dim2Triangle;
typedef Triangulation<2> Dim2Triangulation;
typedef Face<2, 0> Dim2Vertex;

/**
 * \weakgroup dim2
 * @{
 */

/**
 * Details how an edge of a 2-manifold triangulation appears within each
 * triangle.
 *
 * This is a specialisation of the generic FaceEmbedding class template;
 * see the documentation for FaceEmbedding (and also Face) for a general
 * overview of how these face-related classes work.
 *
 * This 2-dimensional specialisation of FaceEmbedding offers additional
 * dimension-specific aliases of some member functions.
 */
template <>
class REGINA_API FaceEmbedding<2, 1> : public FaceEmbeddingBase<2, 1> {
    public:
        /**
         * Default constructor.  This object is unusable until it has
         * some data assigned to it using <tt>operator =</tt>.
         *
         * \ifacespython Not present.
         */
        FaceEmbedding();

        /**
         * Creates a new object containing the given data.
         *
         * @param tri the triangle in which the underlying edge
         * of the triangulation is contained.
         * @param edge the corresponding edge number of \a tri.
         * This must be between 0 and 2 inclusive.
         */
        FaceEmbedding(Dim2Triangle* tri, int edge);

        /**
         * Creates a new copy of the given object.
         *
         * @param cloneMe the object to copy.
         */
        FaceEmbedding(const FaceEmbedding& cloneMe);

        /**
         * A dimension-specific alias for getSimplex().
         *
         * See getSimplex() for further information.
         */
        Dim2Triangle* getTriangle() const;

        /**
         * A dimension-specific alias for getFace().
         *
         * See getFace() for further information.
         */
        int getEdge() const;
};

/**
 * A convenience typedef for FaceEmbedding<2, 1>.
 */
typedef FaceEmbedding<2, 1> Dim2EdgeEmbedding;

/**
 * Represents an edge in the skeleton of a 2-manifold triangulation.
 *
 * This is a specialisation of the generic Face class template; see the
 * documentation for Face for a general overview of how this class works.
 *
 * These specialisations for Regina's \ref stddim "standard dimensions",
 * offer significant extra functionality.
 */
template <>
class REGINA_API Face<2, 1> : public FaceBase<2, 1>, public Output<Face<2, 1>> {
    private:
        /**
         * An array that hard-codes the results of ordering().
         *
         * See ordering() for further details.
         */
        static const NPerm3 ordering_[3];

        Dim2BoundaryComponent* boundaryComponent_;
            /**< The boundary component that this edge is a part of,
                 or 0 if this edge is internal. */

    public:
        /**
         * Returns the boundary component of the triangulation to which
         * this edge belongs.
         *
         * @return the boundary component containing this edge, or 0
         * if this edge does not lie entirely within the boundary of
         * the triangulation.
         */
        Dim2BoundaryComponent* getBoundaryComponent() const;

        /**
         * Returns the vertex of the 2-manifold triangulation corresponding
         * to the given vertex of this edge.
         *
         * @param vertex the vertex of this edge to examine.  This
         * should be either 0 or 1.
         * @return the corresponding vertex of the 2-manifold triangulation.
         */
        Dim2Vertex* getVertex(int vertex) const;

        /**
         * Determines if this edge lies entirely on the boundary of the
         * triangulation.
         *
         * @return \c true if and only if this edge lies on the boundary.
         */
        bool isBoundary() const;

        /**
         * Determines whether this edge represents a dual edge in the
         * maximal forest that has been chosen for the dual 1-skeleton of the
         * triangulation.
         *
         * When the skeletal structure of a triangulation is first computed,
         * a maximal forest in the dual 1-skeleton of the triangulation is
         * also constructed.  Each dual edge in this maximal forest
         * represents a (transverse) edge in the primal skeleton of the
         * triangulation.
         *
         * This maximal forest will remain fixed until the triangulation
         * changes, at which point it will be recomputed (as will all
         * other skeletal objects, such as connected components and so on).
         * There is no guarantee that, when it is recomputed, the
         * maximal forest will use the same dual edges as before.
         *
         * This routine identifies whether this edge corresponds to a
         * member of this dual forest.  In this sense it performs a similar
         * role to Simplex::facetInMaximalForest(), but this routine is
         * typically easier to use.
         *
         * If the skeleton has already been computed, then this routine is
         * very fast (since it just returns a precomputed answer).
         *
         * @return \c true if and only if this edge represents a
         * dual edge in the maximal forest.
         */
        bool inMaximalForest() const;

        /**
         * Given an edge number within a triangle, returns the canonical
         * ordering of those triangle vertices that make up the given edge.
         *
         * This means that the vertices of edge \a i in a triangle are,
         * in canonical order, <tt>ordering[i][0,1]</tt>.
         *
         * Regina defines canonical order to be \e increasing order.
         * That is, <tt>ordering[i][0] &lt; ordering[i][1]</tt>.
         *
         * This routine does \e not describe the mapping from edges
         * of the triangulation into individual triangles; for that, see
         * the routine Dim2Triangle::getEdgeMapping().  Instead, this routine
         * just provides a neat and consistent way of listing the vertices of
         * any given edge of any given triangle.
         *
         * @param edge identifies which edge of a triangle to query.
         * This must be between 0 and 2 inclusive.
         * @return the canonical ordering of the triangle vertices that
         * make up the given triangle edge.
         */
        static NPerm<3> ordering(unsigned edge);

        /**
         * Writes a short text representation of this object to the
         * given output stream.
         *
         * \ifacespython Not present.
         *
         * @param out the output stream to which to write.
         */
        void writeTextShort(std::ostream& out) const;
        /**
         * Writes a detailed text representation of this object to the
         * given output stream.
         *
         * \ifacespython Not present.
         *
         * @param out the output stream to which to write.
         */
        void writeTextLong(std::ostream& out) const;

    private:
        /**
         * Creates a new edge and marks it as belonging to the
         * given triangulation component.
         *
         * @param component the triangulation component to which this
         * edge belongs.
         */
        Face(Dim2Component* component);

    friend class Triangulation<2>;
    friend class TriangulationBase<2>;
};

/**
 * A convenience typedef for Face<2, 1>.
 */
typedef Face<2, 1> Dim2Edge;

/*@}*/

} // namespace regina
// Some more headers that are required for inline functions:
#include "dim2/dim2triangle.h"
namespace regina {

// Inline functions for Dim2EdgeEmbedding

inline FaceEmbedding<2, 1>::FaceEmbedding() :
        FaceEmbeddingBase<2, 1>() {
}

inline FaceEmbedding<2, 1>::FaceEmbedding(
        Dim2Triangle* tri, int edge) :
        FaceEmbeddingBase<2, 1>(tri, edge) {
}

inline FaceEmbedding<2, 1>::FaceEmbedding(
        const Dim2EdgeEmbedding& cloneMe) :
        FaceEmbeddingBase<2, 1>(cloneMe) {
}

inline Dim2Triangle* FaceEmbedding<2, 1>::getTriangle() const {
    return getSimplex();
}

inline int FaceEmbedding<2, 1>::getEdge() const {
    return getFace();
}

// Inline functions for Dim2Edge

inline Face<2, 1>::Face(Dim2Component* component) :
        FaceBase<2, 1>(component), boundaryComponent_(0) {
}

inline Dim2BoundaryComponent* Face<2, 1>::getBoundaryComponent() const {
    return boundaryComponent_;
}

inline Dim2Vertex* Face<2, 1>::getVertex(int vertex) const {
    return front().getTriangle()->getVertex(front().getVertices()[vertex]);
}

inline bool Face<2, 1>::isBoundary() const {
    return (boundaryComponent_ != 0);
}

inline bool Face<2, 1>::inMaximalForest() const {
    return front().getTriangle()->facetInMaximalForest(front().getEdge());
}

inline NPerm<3> Face<2, 1>::ordering(unsigned edge) {
    return ordering_[edge];
}

inline void Face<2, 1>::writeTextShort(std::ostream& out) const {
    out << (boundaryComponent_ ? "Boundary " : "Internal ") << "edge";
}

} // namespace regina

#endif


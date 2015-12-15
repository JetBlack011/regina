
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

#include <algorithm>
#include <queue>
#include <utility>
#include "dim2/dim2triangulation.h"

namespace regina {

void Triangulation<2>::calculateSkeleton() {
    TriangulationBase<2>::calculateSkeleton();

    // Get rid of the empty triangulation now, so that all the helper routines
    // can happily assume at least one triangle.
    if (simplices_.empty())
        return;

    // Off we go!
    calculateVertices();
        // Sets:
        // - vertices_
        // - Dim2Component::vertices_
        // - Dim2Triangle::vertex_
        // - Dim2Triangle::vertexMapping_
        // - all Dim2Vertex members except boundaryComponent_

    calculateBoundary();
        // Sets:
        // - boundaryComponents_
        // - Dim2Component::boundaryComponents_
        // - Dim2 [ Edge, Vertex ]::boundaryComponent_
        // - all Dim2BoundaryComponent members

    // Flesh out the details of each component.
    for (auto e : getEdges())
        e->component()->edges_.push_back(e);
}

void Triangulation<2>::calculateVertices() {
    TriangleIterator it;
    int loopVtx;
    for (it = simplices_.begin(); it != simplices_.end(); ++it)
        for (loopVtx = 0; loopVtx < 3; ++loopVtx)
            (*it)->vertex_[loopVtx] = 0;

    Dim2Vertex* label;
    Dim2Triangle *loopTri, *tri, *adjTri;
    int adjVertex;
    NPerm3 map, adjMap;
    int dir, exitEdge;
    for (it = simplices_.begin(); it != simplices_.end(); ++it) {
        loopTri = *it;
        for (loopVtx = 2; loopVtx >= 0; --loopVtx) {
            if (loopTri->vertex_[loopVtx])
                continue;

            label = new Dim2Vertex(loopTri->component_);
            FaceList<2, 0>::push_back(label);
            loopTri->component_->vertices_.push_back(label);

            // Since triangle vertices are joined together in a loop, the
            // depth-first search is really just a straight line in either
            // direction.  We therefore do away with the usual stack and
            // just keep track of the next triangle to process in the current
            // direction.
            loopTri->vertex_[loopVtx] = label;
            loopTri->vertexMapping_[loopVtx] =
                NPerm3(loopVtx, (loopVtx + 1) % 3, (loopVtx + 2) % 3);
            label->push_back(Dim2VertexEmbedding(loopTri, loopVtx));

            for (dir = 0; dir < 2; ++dir) {
                // Start at the start and walk in one particular direction.
                tri = loopTri;
                map = tri->vertexMapping_[loopVtx];

                while (true) {
                    // Move through to the next triangle.
                    exitEdge = map[dir == 0 ? 1 : 2];
                    adjTri = tri->adjacentTriangle(exitEdge);
                    if (! adjTri)
                        break;

                    adjMap = tri->adjacentGluing(exitEdge) * map *
                        NPerm3(0, 2, 1);
                    adjVertex = adjMap[0];

                    if (adjTri->vertex_[adjVertex]) {
                        // We looped right around.
                        break;
                    }

                    // We have not yet seen this triangle vertex.  Label it.
                    adjTri->vertex_[adjVertex] = label;
                    adjTri->vertexMapping_[adjVertex] = adjMap;

                    if (dir == 0)
                        label->push_back(Dim2VertexEmbedding(
                            adjTri, adjVertex));
                    else
                        label->push_front(Dim2VertexEmbedding(
                            adjTri, adjVertex));

                    tri = adjTri;
                    map = adjMap;
                }
            }
        }
    }
}

void Triangulation<2>::calculateBoundary() {
    // Are there any boundary edges at all?
    long nBdry = 2 * getNumberOfEdges() - 3 * simplices_.size();
    if (nBdry == 0)
        return;

    Dim2BoundaryComponent* label;
    Dim2Triangle *tri, *adjTri;
    int edgeId, adjEdgeId;
    int vertexId, adjVertexId;
    Dim2Edge *adjEdge;
    Dim2Vertex* vertex;
    Dim2VertexEmbedding vertexEmb;

    for (Dim2Edge* edge : getEdges()) {
        // We only care about boundary edges that we haven't yet seen..
        if (edge->getDegree() == 2 || edge->boundaryComponent_)
            continue;

        label = new Dim2BoundaryComponent();
        boundaryComponents_.push_back(label);
        edge->getComponent()->boundaryComponents_.push_back(label);

        // Loop around from this boundary edge to
        // completely enumerate all edges in this boundary component.

        tri = edge->front().getTriangle();
        edgeId = edge->front().getEdge();
        vertexId = edge->front().getVertices()[0];
        vertex = tri->vertex_[vertexId];
        while (true) {
            if (! edge->boundaryComponent_) {
                edge->boundaryComponent_ = label;
                label->edges_.push_back(edge);

                vertex->boundaryComponent_ = label;
                label->vertices_.push_back(vertex);
            } else {
                // We've looped right around.
                break;
            }

            // Find the next edge along the boundary.

            // We can be clever about this.  The current
            // boundary edge is one end of the vertex link; the
            // *adjacent* boundary edge must be at the other.
            vertexEmb = vertex->front();
            if (vertexEmb.getTriangle() == tri &&
                    vertexEmb.getVertices()[0] == vertexId &&
                    vertexEmb.getVertices()[2] == edgeId) {
                // We are currently looking at the embedding at the
                // front of the list.  Take the one at the back.
                vertexEmb = vertex->back();

                adjTri = vertexEmb.getTriangle();
                adjEdgeId = vertexEmb.getVertices()[1];
                adjEdge = adjTri->SimplexFaces<2, 1>::face_[adjEdgeId];
                adjVertexId = vertexEmb.getVertices()[2];
            } else {
                // We must be looking at the embedding at the back
                // of the list.  Take the one at the front (which is
                // already stored in vertexEmb).
                adjTri = vertexEmb.getTriangle();
                adjEdgeId = vertexEmb.getVertices()[2];
                adjEdge = adjTri->SimplexFaces<2, 1>::face_[adjEdgeId];
                adjVertexId = vertexEmb.getVertices()[1];

                // TODO: Sanity checking; remove this eventually.
                vertexEmb = vertex->back();
                if (! (vertexEmb.getTriangle() == tri &&
                        vertexEmb.getVertices()[0] == vertexId &&
                        vertexEmb.getVertices()[1] == edgeId)) {
                    std::cerr << "ERROR: Something has gone terribly "
                        "wrong in computeBoundaryComponents()."
                        << std::endl;
                    ::exit(1);
                }
            }

            edge = adjEdge;
            tri = adjTri;
            edgeId = adjEdgeId;
            vertexId = adjVertexId;
            vertex = tri->vertex_[vertexId];
        }
    }
}

} // namespace regina


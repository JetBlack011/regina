
/**************************************************************************
 *                                                                        *
 *  Regina - A Normal Surface Theory Calculator                           *
 *  Computational Engine                                                  *
 *                                                                        *
 *  Copyright (c) 1999-2013, Ben Burton                                   *
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

#include <cassert>
#include <iomanip>
#include <iostream>
#include <sstream>
#include "dim2/dim2triangulation.h"
#include "utilities/xmlutils.h"

namespace regina {

Dim2Triangulation::Dim2Triangulation(const std::string& description) :
        calculatedSkeleton_(false) {
    Dim2Triangulation* attempt;

    if ((attempt = fromIsoSig(description))) {
        cloneFrom(*attempt);
        setPacketLabel(description);
    }

    delete attempt;
}

bool Dim2Triangulation::isIdenticalTo(const Dim2Triangulation& other) const {
    if (triangles_.size() != other.triangles_.size())
        return false;

    unsigned long i;
    unsigned j;
    for (i = 0; i < triangles_.size(); ++i)
        for (j = 0; j < 3; ++j) {
            if (triangles_[i]->adjacentTriangle(j)) {
                if (! other.triangles_[i]->adjacentTriangle(j))
                    return false;
                if (triangles_[i]->adjacentTriangle(j)->markedIndex() !=
                        other.triangles_[i]->adjacentTriangle(j)->markedIndex())
                    return false;
                if (triangles_[i]->adjacentGluing(j) !=
                        other.triangles_[i]->adjacentGluing(j))
                    return false;
            } else {
                if (other.triangles_[i]->adjacentTriangle(j))
                    return false;
            }
        }

    return true;
}

bool Dim2Triangulation::isMinimal() const {
    // 2-sphere:
    if (getEulerChar() == 2)
        return (triangles_.size() == 2);

    // Projective plane and disc:
    if (getEulerChar() == 1)
        return (triangles_.size() == (isClosed() ? 2 : 1));

    // All other closed manifolds:
    if (isClosed())
        return (vertices_.size() == 1);

    // All other bounded manifolds:
    return (vertices_.size() == boundaryComponents_.size());
}

void Dim2Triangulation::swapContents(Dim2Triangulation& other) {
    ChangeEventSpan span1(this);
    ChangeEventSpan span2(&other);

    clearAllProperties();
    other.clearAllProperties();

    triangles_.swap(other.triangles_);

    TriangleIterator it;
    for (it = triangles_.begin(); it != triangles_.end(); ++it)
        (*it)->tri_ = this;
    for (it = other.triangles_.begin(); it != other.triangles_.end(); ++it)
        (*it)->tri_ = &other;
}

void Dim2Triangulation::moveContentsTo(Dim2Triangulation& dest) {
    ChangeEventSpan span1(this);
    ChangeEventSpan span2(&dest);

    clearAllProperties();
    dest.clearAllProperties();

    TriangleIterator it;
    for (it = triangles_.begin(); it != triangles_.end(); ++it) {
        // This is an abuse of NMarkedVector, since for a brief moment
        // each triangle belongs to both vectors triangles_ and dest.triangles_.
        // However, the subsequent clear() operation does not touch the
        // triangle markings (indices), and so we end up with the
        // correct result (i.e., the markings are correct for dest).
        (*it)->tri_ = &dest;
        dest.triangles_.push_back(*it);
    }
    triangles_.clear();
}

void Dim2Triangulation::writeTextLong(std::ostream& out) const {
    if (! calculatedSkeleton_)
        calculateSkeleton();

    out << "Size of the skeleton:\n";
    out << "  Triangles: " << triangles_.size() << '\n';
    out << "  Edges: " << edges_.size() << '\n';
    out << "  Vertices: " << vertices_.size() << '\n';
    out << '\n';

    Dim2Triangle* tri;
    Dim2Triangle* adjTri;
    unsigned triPos;
    int i, j;
    NPerm3 adjPerm;

    out << "Triangle gluing:\n";
    out << "  Triangle  |  glued to:     (01)     (02)     (12)\n";
    out << "  ----------+--------------------------------------\n";
    for (triPos=0; triPos < triangles_.size(); triPos++) {
        tri = triangles_[triPos];
        out << "      " << std::setw(4) << triPos << "  |           ";
        for (i = 2; i >= 0; --i) {
            out << " ";
            adjTri = tri->adjacentTriangle(i);
            if (! adjTri)
                out << "boundary";
            else {
                adjPerm = tri->adjacentGluing(i);
                out << std::setw(3) << triangleIndex(adjTri) << " (";
                for (j = 0; j < 3; ++j) {
                    if (j == i) continue;
                    out << adjPerm[j];
                }
                out << ")";
            }
        }
        out << '\n';
    }
    out << '\n';

    out << "Vertices:\n";
    out << "  Triangle  |  vertex:    0   1   2\n";
    out << "  ----------+----------------------\n";
    for (triPos = 0; triPos < triangles_.size(); ++triPos) {
        tri = triangles_[triPos];
        out << "      " << std::setw(4) << triPos << "  |          ";
        for (i = 0; i < 3; ++i)
            out << ' ' << std::setw(3) <<
                vertexIndex(tri->getVertex(i));
        out << '\n';
    }
    out << '\n';

    out << "Edges:\n";
    out << "  Triangle  |  edge:   01  02  12\n";
    out << "  ----------+--------------------\n";
    for (triPos = 0; triPos < triangles_.size(); ++triPos) {
        tri = triangles_[triPos];
        out << "      " << std::setw(4) << triPos << "  |        ";
        for (i = 2; i >= 0; --i)
            out << ' ' << std::setw(3) << edgeIndex(tri->getEdge(i));
        out << '\n';
    }
    out << '\n';
}

void Dim2Triangulation::insertTriangulation(const Dim2Triangulation& X) {
    ChangeEventSpan span(this);

    unsigned long nOrig = getNumberOfTriangles();
    unsigned long nX = X.getNumberOfTriangles();

    unsigned long triPos;
    for (triPos = 0; triPos < nX; ++triPos)
        newTriangle(X.triangles_[triPos]->getDescription());

    // Make the gluings.
    unsigned long adjPos;
    Dim2Triangle* tri;
    Dim2Triangle* adjTri;
    NPerm3 adjPerm;
    int edge;
    for (triPos = 0; triPos < nX; ++triPos) {
        tri = X.triangles_[triPos];
        for (edge = 0; edge < 3; ++edge) {
            adjTri = tri->adjacentTriangle(edge);
            if (adjTri) {
                adjPos = X.triangleIndex(adjTri);
                adjPerm = tri->adjacentGluing(edge);
                if (adjPos > triPos ||
                        (adjPos == triPos && adjPerm[edge] > edge)) {
                    triangles_[nOrig + triPos]->joinTo(edge,
                        triangles_[nOrig + adjPos], adjPerm);
                }
            }
        }
    }
}

void Dim2Triangulation::insertConstruction(unsigned long nTriangles,
        const int adjacencies[][3], const int gluings[][3][3]) {
    if (nTriangles == 0)
        return;

    Dim2Triangle** tri = new Dim2Triangle*[nTriangles];

    unsigned i, j;
    NPerm3 p;

    ChangeEventSpan span(this);

    for (i = 0; i < nTriangles; ++i)
        tri[i] = newTriangle();

    for (i = 0; i < nTriangles; ++i)
        for (j = 0; j < 3; ++j)
            if (adjacencies[i][j] >= 0 &&
                    ! tri[i]->adjacentTriangle(j)) {
                p = NPerm3(gluings[i][j][0], gluings[i][j][1],
                    gluings[i][j][2]);
                tri[i]->joinTo(j, tri[adjacencies[i][j]], p);
            }

    delete[] tri;
}

std::string Dim2Triangulation::dumpConstruction() const {
    std::ostringstream ans;
    ans <<
"/**\n";
    if (! getPacketLabel().empty())
        ans <<
" * 2-manifold triangulation: " << getPacketLabel() << "\n";
    ans <<
" * Code automatically generated by dumpConstruction().\n"
" */\n"
"\n";

    if (triangles_.empty()) {
        ans <<
"/* This triangulation is empty.  No code is being generated. */\n";
        return ans.str();
    }

    ans <<
"/**\n"
" * The following arrays describe the individual gluings of\n"
" * triangle edges.\n"
" */\n"
"\n";

    unsigned long nTriangles = triangles_.size();
    Dim2Triangle* tri;
    NPerm3 perm;
    unsigned long p;
    int e, i;

    ans << "const int adjacencies[" << nTriangles << "][3] = {\n";
    for (p = 0; p < nTriangles; ++p) {
        tri = triangles_[p];

        ans << "    { ";
        for (e = 0; e < 3; ++e) {
            if (tri->adjacentTriangle(e)) {
                ans << triangleIndex(tri->adjacentTriangle(e));
            } else
                ans << "-1";

            if (e < 2)
                ans << ", ";
            else if (p != nTriangles - 1)
                ans << "},\n";
            else
                ans << "}\n";
        }
    }
    ans << "};\n\n";

    ans << "const int gluings[" << nTriangles << "][3][3] = {\n";
    for (p = 0; p < nTriangles; ++p) {
        tri = triangles_[p];

        ans << "    { ";
        for (e = 0; e < 3; ++e) {
            if (tri->adjacentTriangle(e)) {
                perm = tri->adjacentGluing(e);
                ans << "{ ";
                for (i = 0; i < 3; ++i) {
                    ans << perm[i];
                    if (i < 2)
                        ans << ", ";
                    else
                        ans << " }";
                }
            } else
                ans << "{ 0, 0, 0 }";

            if (e < 2)
                ans << ", ";
            else if (p != nTriangles - 1)
                ans << " },\n";
            else
                ans << " }\n";
        }
    }
    ans << "};\n\n";

    ans <<
"/**\n"
" * The following code actually constructs a 2-manifold triangulation\n"
" * based on the information stored in the arrays above.\n"
" */\n"
"\n"
"Dim2Triangulation tri;\n"
"tri.insertConstruction(" << nTriangles << ", adjacencies, gluings);\n"
"\n";

    return ans.str();
}

void Dim2Triangulation::writeXMLPacketData(std::ostream& out) const {
    using regina::xml::xmlEncodeSpecialChars;
    using regina::xml::xmlValueTag;

    // Write the triangle gluings.
    TriangleIterator it;
    Dim2Triangle* adjTri;
    int edge;

    out << "  <triangles ntriangles=\"" << triangles_.size() << "\">\n";
    for (it = triangles_.begin(); it != triangles_.end(); ++it) {
        out << "    <triangle desc=\"" <<
            xmlEncodeSpecialChars((*it)->getDescription()) << "\"> ";
        for (edge = 0; edge < 3; ++edge) {
            adjTri = (*it)->adjacentTriangle(edge);
            if (adjTri) {
                out << triangleIndex(adjTri) << ' '
                    << static_cast<int>((*it)->
                        adjacentGluing(edge).getPermCode()) << ' ';
            } else
                out << "-1 -1 ";
        }
        out << "</triangle>\n";
    }
    out << "  </triangles>\n";
}

void Dim2Triangulation::cloneFrom(const Dim2Triangulation& X) {
    ChangeEventSpan span(this);

    removeAllTriangles();

    TriangleIterator it;
    for (it = X.triangles_.begin(); it != X.triangles_.end(); ++it)
        newTriangle((*it)->getDescription());

    // Make the gluings.
    long triPos, adjPos;
    Dim2Triangle* tri;
    Dim2Triangle* adjTri;
    NPerm3 adjPerm;
    int edge;
    triPos = 0;
    for (it = X.triangles_.begin(); it != X.triangles_.end(); ++it) {
        tri = *it;
        for (edge = 0; edge < 3; ++edge) {
            adjTri = tri->adjacentTriangle(edge);
            if (adjTri) {
                adjPos = X.triangleIndex(adjTri);
                adjPerm = tri->adjacentGluing(edge);
                if (adjPos > triPos ||
                        (adjPos == triPos && adjPerm[edge] > edge)) {
                    triangles_[triPos]->joinTo(edge,
                        triangles_[adjPos], adjPerm);
                }
            }
        }
        ++triPos;
    }

    // Properties:
    // None yet for 2-manifold triangulations.
}

void Dim2Triangulation::deleteTriangles() {
    for (TriangleIterator it = triangles_.begin(); it != triangles_.end(); ++it)
        delete *it;
    triangles_.clear();
}

void Dim2Triangulation::deleteSkeleton() {
    for (VertexIterator it = vertices_.begin(); it != vertices_.end(); ++it)
        delete *it;
    for (EdgeIterator it = edges_.begin(); it != edges_.end(); ++it)
        delete *it;
    for (ComponentIterator it = components_.begin();
            it != components_.end(); ++it)
        delete *it;
    for (BoundaryComponentIterator it = boundaryComponents_.begin();
            it != boundaryComponents_.end(); ++it)
        delete *it;

    vertices_.clear();
    edges_.clear();
    components_.clear();
    boundaryComponents_.clear();

    calculatedSkeleton_ = false;
}

void Dim2Triangulation::clearAllProperties() {
    if (calculatedSkeleton_)
        deleteSkeleton();
}


bool Dim2Triangulation::oneThreeMove( Dim2Triangle *t, bool check, bool perform )
{
    // not yet okay-ed by bab
    #ifdef DEBUG
    std::cerr << "Performing 1-3 move\n";
    #endif

    if (!perform) return true; // this move can always be done.

    // two calls to newTriangle( std::string desc )
    Dim2Triangle *newTu = newTriangle();
    Dim2Triangle *newTl = newTriangle();
    // glue these two triangles together
    newTu->joinTo(1, newTl, NPerm3( 0, 2, 1 ) );

    // get info on the triangles across edges 1 and 2 from t before we disconnect those gluings
    Dim2Triangle *triu( t->adjacentTriangle(1) );
    Dim2Triangle *tril( t->adjacentTriangle(2) );
    NPerm3 glueU(t->adjacentGluing(1));
    NPerm3 glueL(t->adjacentGluing(2));

    // unjoin twice, 
    t->unjoin(1); t->unjoin(2);

    // 4 joins with the pre-existing triangles and new triangles
    t->joinTo(1, newTl, NPerm3( 1, 0, 2 ) );
    t->joinTo(2, newTu, NPerm3( 2, 1, 0 ) );
    newTu->joinTo(2, triu, NPerm3( glueU[0], glueU[1], glueU[2] ) );
    newTl->joinTo(1, tril, NPerm3( glueL[0], glueL[1], glueL[2] ) );

    // TODO: put in the labels of the new triangles which triangles, edges and vertex
    //       indices are "new" vs "old" so that we can keep track in the circle packing
    //       code which introductions were artificial. 

    return true; 
}

} // namespace regina

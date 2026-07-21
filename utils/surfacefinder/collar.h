//
//  collar.h
//
//  Traces a knot/link's edges through CobordismBuilder<3>'s thicken()
//  layers, producing the collar of triangles that connects the edges'
//  original position to wherever thickening has pushed them -- a fixed
//  base suitable for SurfaceFinder<4>::findSurfaces(startingTriangles),
//  so the search only has to find a cap rather than rediscovering the
//  collar itself from scratch among everything else in the cobordism.
//

#ifndef COLLAR_H

#define COLLAR_H

#include <array>
#include <unordered_set>
#include <vector>

#include <triangulation/dim3.h>
#include <triangulation/dim4.h>

#include "cobordismbuilder.h"

// IMPORTANT (both learned the hard way -- see class usage notes below):
//
// 1. CobordismBuilder<3>'s constructor takes a *copy* of whatever
//    triangulation it's given (and may reorder it), so a Tetrahedron<3>*/
//    Edge<3>* obtained from the caller's own triangulation is a dangling
//    reference to a different object as far as the CobordismBuilder is
//    concerned. Only indices survive the copy (and any subsequent order())
//    unchanged; CollarBuilder therefore takes edge *indices*, to be
//    resolved against CobordismBuilder::baseTriangulation() rather than
//    against whatever triangulation the caller originally built the knot
//    in.
//
// 2. Face<dim,subdim>* objects (Triangle<4>*, Edge<4>*, ...) -- as opposed
//    to Simplex<dim>* itself -- go dangling the moment the triangulation
//    is further modified and its skeleton gets recomputed, even though the
//    Simplex<dim>* they came from stays valid. So CollarBuilder never
//    stores a Triangle<4>*; it stores (Simplex<4>*, local vertex triple)
//    descriptors during construction and only resolves them to actual
//    Triangle<4>* pointers in resolve(), which must be called after every
//    thicken()/cone() call on the CobordismBuilder is done.
class CollarBuilder {
  private:
    struct FaceDesc {
        regina::Simplex<4> *simplex;
        std::array<int, 3> verts;
    };

    std::vector<int> edgeIndices_;
    std::vector<FaceDesc> faces_;

  public:
    // edgeIndices: indices (within cob.baseTriangulation(), NOT within
    // whatever triangulation the caller used to build the knot/link) of
    // the base edges to trace.
    explicit CollarBuilder(std::vector<int> edgeIndices)
        : edgeIndices_(std::move(edgeIndices)) {}

    // Call once per thicken() layer, immediately after that thicken() call
    // and before the next one -- CobordismBuilder only retains the most
    // recently built layer's prisms (currentTopSimplex()), so a layer's
    // sweep must be extracted before it's superseded.
    void addLayer(const CobordismBuilder<3> &cob);

    // Resolves every accumulated face descriptor into an actual
    // Triangle<4>*. Must only be called after all thicken()/cone() calls
    // on `cob` are complete (see class-level note #2).
    std::unordered_set<regina::Triangle<4> *> resolve() const;

  private:
    // Local vertex index (within a SimplicialPrism<4> sub-simplex) holding
    // (v, isTop) -- matches SimplicialPrism<4>::localVertex(), duplicated
    // here since that's a static method on a class this header would
    // otherwise need to instantiate a whole prism just to call.
    static int encodeVertex_(int v, bool isTop);

    // The two triangles tracing base edge {a,b} (a<b, tetrahedron-local
    // vertex indices) as it sweeps through the most recently built
    // thickening layer. See SimplicialPrism's decode_/encode_ for the
    // staircase encoding this is derived from: for a<b, the sweep is
    // {(a,bottom),(b,bottom),(a,top)} in the layer's a-th sub-simplex, and
    // {(b,bottom),(a,top),(b,top)} in its b-th sub-simplex.
    void addEdgeSweep_(const CobordismBuilder<3> &cob,
                       regina::Tetrahedron<3> *tet, int a, int b);

    static regina::Triangle<4> *resolveFace_(const FaceDesc &f);
};

#endif // COLLAR_H

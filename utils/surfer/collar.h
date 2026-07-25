//
//  collar.h
//
//  Created by John Teague on 07/05/2026.
//

#ifndef COLLAR_H

#define COLLAR_H

#include <array>
#include <unordered_set>
#include <vector>

#include <triangulation/dim3.h>
#include <triangulation/dim4.h>

#include "cobordismbuilder.h"

/*! \file utils/surfer/collar.h
 *  \brief Traces a knot/link through a cobordism's thickening layers.
 */

/**
 * Traces a knot/link's edges through CobordismBuilder<3>'s thicken()
 * layers, producing the collar of triangles connecting the edges'
 * original position to wherever thickening has pushed them.
 *
 * This is a fixed base for seeding SurfaceSearch (see embeddingsearch.h),
 * so the search only has to find a cap rather than rediscovering the
 * collar itself from scratch among everything else in the cobordism.
 *
 * \warning See CobordismBuilder::baseTriangulation()'s documentation:
 * its constructor takes a *copy* of the given triangulation, so pointers
 * obtained from the caller's own triangulation go stale. CollarBuilder
 * therefore takes edge *indices*, resolved against
 * `cob.baseTriangulation()` rather than the triangulation the caller
 * originally built the knot in.
 *
 * \warning Face<dim,subdim>* objects (Triangle<4>*, Edge<4>*, ...), unlike
 * Simplex<dim>* itself, go dangling the moment the triangulation is
 * further modified and its skeleton recomputed. So CollarBuilder never
 * stores a Triangle<4>*; it stores (Simplex<4>*, local vertex triple)
 * descriptors, and only resolves them to actual Triangle<4>* pointers in
 * resolve(), which must be called after every thicken()/cone() call on
 * the CobordismBuilder is done.
 */
class CollarBuilder {
  private:
    /** A triangle of the cobordism, described positionally (see the class \warning above). */
    struct FaceDesc {
        regina::Simplex<4> *simplex;
        std::array<int, 3> verts;
    };

    std::vector<int> edgeIndices_;
        /**< Indices, within cob.baseTriangulation(), of the base edges being traced. */
    std::vector<FaceDesc> faces_; /**< Accumulated face descriptors, from every addLayer() call so far. */

  public:
    /**
     * Traces the base edges at `edgeIndices` (indices within
     * cob.baseTriangulation(), not the triangulation the caller used to
     * build the knot/link).
     */
    explicit CollarBuilder(std::vector<int> edgeIndices)
        : edgeIndices_(std::move(edgeIndices)) {}

    /**
     * Extracts the sweep of the tracked edges through the most recently
     * built thickening layer.
     *
     * \pre Called once per thicken() layer, immediately after that
     * thicken() call and before the next -- CobordismBuilder only retains
     * the most recently built layer's prisms (currentTopSimplex()), so a
     * layer's sweep must be extracted before it's superseded.
     */
    void addLayer(const CobordismBuilder<3> &cob);

    /**
     * Resolves every accumulated face descriptor into an actual
     * Triangle<4>*.
     *
     * \pre All thicken()/cone() calls on `cob` are complete (see the
     * second class-level \warning above).
     */
    std::unordered_set<regina::Triangle<4> *> resolve() const;

  private:
    /**
     * The local vertex index (within a SimplicialPrism<4> sub-simplex)
     * holding (v, isTop) -- matches SimplicialPrism<4>::localVertex(),
     * duplicated here to avoid instantiating a whole prism just for this.
     */
    static int encodeVertex_(int v, bool isTop);

    /**
     * Records the two triangles tracing base edge {a,b} (a < b,
     * tetrahedron-local vertex indices) as it sweeps through the most
     * recently built thickening layer.
     *
     * See SimplicialPrism's decode_()/encode_() for the staircase
     * encoding this is derived from: for a < b, the sweep is
     * {(a,bottom),(b,bottom),(a,top)} in the layer's a-th sub-simplex,
     * and {(b,bottom),(a,top),(b,top)} in its b-th.
     */
    void addEdgeSweep_(const CobordismBuilder<3> &cob,
                       regina::Tetrahedron<3> *tet, int a, int b);

    /** Resolves a single FaceDesc into its actual Triangle<4>*. */
    static regina::Triangle<4> *resolveFace_(const FaceDesc &f);
};

#endif // COLLAR_H

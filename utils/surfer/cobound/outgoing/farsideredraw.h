//
//  farsideredraw.h
//
//  A stored witness's outgoing link, recovered from its pair signature.
//

/*! \file utils/surfer/cobound/outgoing/farsideredraw.h
 *  \brief Redraws witnesses of one row from their pair signatures, exactly
 *  as the search itself would have seen them.
 *
 *  The row is thickened as verifyslicegenus thickens it
 *  (rowsearch::buildRow()). A witness's decoded pair is
 *  carried onto that thickening by an isomorphism sending its incoming curve
 *  onto the row's L x {0}; from there the outgoing side is read through the
 *  thickening's own product structure (farside::OutgoingMap), so no
 *  automorphism of T can mirror it, and oriented against the row
 *  (farside::orientedOutgoingLink()). The isomorphism is the only search,
 *  and L x {0} pins it.
 *
 *  Used by farsidediagram (diagrams) and farsidename (exact names).
 */

#ifndef SURFER_FARSIDEREDRAW_H
#define SURFER_FARSIDEREDRAW_H

#include <memory>
#include <optional>
#include <unordered_map>
#include <string>
#include <vector>

#include <triangulation/dim3.h>
#include <triangulation/dim4.h>

#include "diagramtriangulation/thickening/thickening.h"
#include "cobound/search/incoming.h"
#include "cobound/outgoing/outgoinglink.h"
#include "diagramtriangulation/todiagram.h"
#include "diagramtriangulation/fromdiagram.h"
#include "cobound/search/rowsearch.h"
#include "surfer/submanifold/skeleton.h"

namespace farside {

class WitnessRedrawer {
  public:
    /**
     * \param rowPD the row's PD code, as the tables write it.
     * \param layers the witnesses' thicken_layers (cobordisms.csv).
     */
    WitnessRedrawer(const std::string &rowPD, int layers);
    WitnessRedrawer(const WitnessRedrawer &) = delete;
    WitnessRedrawer &operator=(const WitnessRedrawer &) = delete;

    /**
     * The witness's surface as triangle indices of the thickening, or
     * nullopt with `why` set when no isomorphism carries its incoming curve
     * onto L x {0}.
     */
    std::optional<std::vector<int>> carry(const std::string &pairsig, std::string &why) const;

    /** The oriented outgoing link of a witness, or nullopt with `why`. */
    std::optional<OutgoingLink> outgoingLink(const std::string &pairsig, std::string &why) const;

    /**
     * As outgoingLink(), without rebuilding what a stored witness no longer
     * needs checked -- and what made the reference path ~1 s per witness:
     *
     *   - A pair signature is the AMBIENT's own isomorphism signature, a
     *     delimiter, and the surface's faces in that ambient's canonical
     *     reconstruction (pairsig.h). Every witness of a row therefore has
     *     the same ambient part: it is decoded, and its isomorphisms onto
     *     the thickening found, once per row; per witness only the face
     *     suffix is decoded, and the cached isomorphisms tried until one
     *     carries the incoming curve onto L x {0}.
     *   - The surface is not rebuilt as a KnottedSurface, whose addFaces()
     *     re-runs the search's embeddedness and local-flatness checks with a
     *     cold cache (~290 ms): a stored witness passed them when found. Its
     *     triangles are glued into a plain regina::Triangulation<2>, each
     *     component oriented, and the boundary read off exactly as
     *     KnottedSurface::orientedBoundaryLinks() does.
     *
     * farsidename checks this against outgoingLink() (--reference).
     */
    std::optional<OutgoingLink> outgoingLinkFast(const std::string &pairsig,
                                                 std::string &why) const;

    /**
     * Rebuilds a surface given by its triangles in thickening() -- as an
     * in-process cascade hop keeps them -- face by face into `surface`, which
     * must be empty, over skeleton(), made with resolveUnlinked on. It runs
     * the search's own checks: every face must add, and the whole must
     * satisfy `proper`, be acceptable (embedded or resolvable, and smooth at
     * the boundary), and meet the incoming boundary in exactly L x {0}.
     * False, with `why`, otherwise.
     */
    bool rebuild(const std::vector<int> &faces, KnottedSurface &surface,
                 std::string &why) const;

    /** As outgoingLink(), for a surface given by its faces (rebuild()). */
    std::optional<OutgoingLink> outgoingLinkFromFaces(const std::vector<int> &faces,
                                                      std::string &why) const;

    /**
     * A digest of thickening(): its size, every gluing, and which pentachoron
     * face each triangle is. Faces recorded against one build are read only
     * in a build with the same digest, so a change to the construction (or to
     * Regina's skeleton numbering) is refused rather than misread.
     */
    std::string buildChecksum() const;

    const regina::Triangulation<3> &knotT() const { return rb_.link.tri; }
    const regina::Triangulation<4> &thickening() const { return rb_.tri; }
    const Skeleton<4, 2> &skeleton() const { return *skeleton_; }
    const OutgoingMap &outgoing() const { return *outgoing_; }
    const cobordismgraph::RowOrientation &row() const { return *rb_.orientation; }
    size_t incomingBC() const { return rb_.searchSideBC; }
    const std::vector<size_t> &rowEdges() const { return rb_.searchEdges; }
    const knotbuilder::DiagramDrawer &drawer() const { return *drawer_; }
    /** The row's own components, in DiagramDrawer::cyclesOf() order. */
    const std::vector<knotbuilder::EdgeCycle> &rowCycles() const { return rowCycles_; }
    /** Which of rowCycles() the search-side edge `edgeIndex` (an index into
     *  the incoming boundary component's built triangulation) lies on.
     *  \throws std::out_of_range for an edge that is not one of L's. */
    size_t rowComponentOf(size_t edgeIndex) const { return rowComponentOf_.at(edgeIndex); }
    const knotbuilder::TriangulationWithLink &built() const { return rb_.link; }
    /** The whole row build: a search run in thickening() (cascadesearch's
     *  in-process hops) sees exactly what this redrawer reads. */
    const rowsearch::RowBuild &rowBuild() const { return rb_; }

    /** Cumulative milliseconds spent decoding pair signatures, and searching
     *  for the isomorphism onto the thickening (carry()). */
    double msDecode() const { return msDecode_; }
    double msIsomorphism() const { return msIso_; }
    /** Of outgoingLink()'s own time: building the KnottedSurface (its
     *  boundary-component builds, then addFaces()), and reading the oriented
     *  outgoing link off it. */
    double msSurfaceBuild() const { return msSurface_; }
    double msBoundaryRead() const { return msRead_; }
    /** Of msSurfaceBuild(): the KnottedSurface constructor alone (it builds
     *  every boundary component of the thickening as a 3-triangulation). */
    double msBoundaryBuild() const { return msBoundaryBuild_; }

  private:
    rowsearch::RowBuild rb_; /**< T, the thickening, its collar and row map. */
    std::unique_ptr<OutgoingMap> outgoing_;
    std::unique_ptr<knotbuilder::DiagramDrawer> drawer_;
    std::unique_ptr<Skeleton<4, 2>> skeleton_;
    std::vector<knotbuilder::EdgeCycle> rowCycles_;
    std::unordered_map<size_t, size_t> rowComponentOf_; /**< see rowComponentOf() */
    mutable double msDecode_ = 0, msIso_ = 0, msSurface_ = 0, msRead_ = 0, msBoundaryBuild_ = 0;

    // The fast path's per-row state (outgoingLinkFast()).
    mutable std::string ambientSig_;
    mutable std::unique_ptr<regina::Triangulation<4>> ambient_;
    mutable std::vector<regina::Isomorphism<4>> isos_;
    mutable int faceWidth_ = 0;
    std::vector<regina::Triangulation<3>> boundaries_; /**< W's boundary components, built */
    std::unordered_map<const regina::Edge<4> *, std::pair<size_t, size_t>> boundaryEdge_;
    /**< an edge of W in its boundary -> (component, local edge index) */
};

} // namespace farside

#endif

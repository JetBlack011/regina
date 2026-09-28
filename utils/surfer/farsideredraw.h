//
//  farsideredraw.h
//
//  A stored witness's outgoing link, recovered from its pair signature.
//

/*! \file utils/surfer/farsideredraw.h
 *  \brief Redraws witnesses of one row from their pair signatures, exactly
 *  as the search itself would have seen them.
 *
 *  The row is thickened as verifyslicegenus thickens it (knotbuilder,
 *  CobordismBuilder x layers, CollarBuilder). A witness's decoded pair is
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

#include "cobordismbuilder.h"
#include "cobordismgraph.h"
#include "farsidecurves.h"
#include "knotbuilder/diagramdrawer.h"
#include "knotbuilder/knotbuilder.h"
#include "skeleton.h"

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

    const regina::Triangulation<3> &knotT() const { return built_.tri; }
    const regina::Triangulation<4> &thickening() const { return cob_->getCobordism(); }
    const Skeleton<4, 2> &skeleton() const { return *skeleton_; }
    const OutgoingMap &outgoing() const { return *outgoing_; }
    const cobordismgraph::RowOrientation &row() const { return row_; }
    size_t incomingBC() const { return incomingBC_; }
    const std::vector<size_t> &rowEdges() const { return rowEdges_; }
    const knotbuilder::DiagramDrawer &drawer() const { return *drawer_; }
    /** The row's own components, in DiagramDrawer::cyclesOf() order. */
    const std::vector<knotbuilder::EdgeCycle> &rowCycles() const { return rowCycles_; }
    const knotbuilder::TriangulationWithLink &built() const { return built_; }

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
    knotbuilder::PDCode pd_;
    knotbuilder::TriangulationWithLink built_;
    std::unique_ptr<CobordismBuilder<3>> cob_;
    size_t incomingBC_ = 0;
    std::vector<size_t> rowEdges_;
    cobordismgraph::RowOrientation row_;
    std::unique_ptr<OutgoingMap> outgoing_;
    std::unique_ptr<knotbuilder::DiagramDrawer> drawer_;
    std::unique_ptr<Skeleton<4, 2>> skeleton_;
    std::vector<knotbuilder::EdgeCycle> rowCycles_;
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

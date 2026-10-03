//
//  fromdatabase.h
//
//  A stored cobordism's outgoing link: read back from its pair signature
//  onto its row's thickening (WitnessRedrawer), and kept across runs
//  (RowReadBacks).
//

/*! \file utils/surfer/cobound/outgoing/fromdatabase.h
 *  \brief Reading stored cobordisms' outgoing links, and keeping them.
 *
 *  **WitnessRedrawer** redraws witnesses of one row from their pair
 *  signatures, exactly as the search itself would have seen them. The row is
 *  thickened as verifyslicegenus thickens it (search::buildRow()). A
 *  witness's decoded pair is carried onto that thickening by an isomorphism
 *  sending its incoming curve onto the row's L x {0}; from there the outgoing
 *  side is read through the thickening's own product structure
 *  (OutgoingMap), so no automorphism of T can mirror it, and oriented against
 *  the row (outgoing::orientedOutgoingLink()). The isomorphism is the only
 *  search, and L x {0} pins it. Used by farsidediagram (diagrams),
 *  farsidename (exact names) and the cascade (stored cobordisms as graph
 *  edges).
 *
 *  **RowReadBacks** keeps those read-backs across runs. Reading a master
 *  witness back (WitnessRedrawer::outgoingLinkFast()) enumerates the
 *  isomorphisms from its pair signature's ambient onto the row's thickening
 *  -- most of a run's single-threaded driver time on 2026-09-30, and the same
 *  witnesses of the same popular rows (L6a5{0;1}, L7a7{0;0}, L10a147{0;0},
 *  ...) were read again by nearly every target's run. The result is fixed by
 *  the witness and the row's exact thickening, so it is kept: one file per
 *  row, one line per witness.
 *
 *    <dir>/<rowkey>.readback       rowkey = witnessKey(row PD + "|" + layers)
 *      readback 1 <build digest>   the row's WitnessRedrawer::buildChecksum()
 *      <witness key>\tok\t<link>   a read-back (serialiseLink())
 *      <witness key>\tfail\t<why>  a read-back that fails (it always will)
 *
 *  A file whose digest is not the row's as built now is ignored and replaced:
 *  knotbuilder or the thickening changed, so its edge numbers mean nothing.
 *  Appends hold an flock on the file; a torn last line is ignored.
 */

#ifndef CASCADE_READBACKCACHE_H
#define CASCADE_READBACKCACHE_H

#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

#include <triangulation/dim3.h>
#include <triangulation/dim4.h>

#include "cobound/outgoing/outgoinglink.h"
#include "cobound/search/incoming.h"
#include "diagramtriangulation/fromdiagram.h"
#include "diagramtriangulation/thickening/thickening.h"
#include "diagramtriangulation/todiagram.h"
#include "surfer/submanifold/skeleton.h"

namespace outgoing {

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
    const search::RowOrientation &row() const { return *rb_.orientation; }
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
    const search::RowBuild &rowBuild() const { return rb_; }

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
    /** Calls its argument on isomorphisms (ambient -> thickening) until it
     *  returns true. */
    using IsoVisitor = std::function<bool(const regina::Isomorphism<4> &)>;
    using IsoSource = std::function<void(const IsoVisitor &)>;

    /**
     * The one carry: `faces` (triangles of `ambient`) carried onto
     * thickening() by the first isomorphism `isos` offers whose image meets
     * the incoming boundary in exactly L x {0}, or nullopt. carry() decodes
     * the whole pair signature and enumerates the isomorphisms as it goes;
     * outgoingLinkFast() decodes only the face suffix and offers the row's
     * isomorphisms, found once.
     */
    std::optional<std::vector<int>> pinned_(const regina::Triangulation<4> &ambient,
                                            const std::vector<int> &faces,
                                            const IsoSource &isos) const;

    search::RowBuild rb_; /**< T, the thickening, its collar and row map. */
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

} // namespace outgoing

namespace outgoing {

/// One cached read-back: the link, or why reading it back fails.
struct CachedReadBack {
  std::optional<outgoing::OutgoingLink> link;
  std::string why;
};

std::string serialiseLink(const outgoing::OutgoingLink &link);
/// Inverse of serialiseLink(); nullopt if the text is malformed.
std::optional<outgoing::OutgoingLink> parseLink(const std::string &text);

/// One row's read-backs: loaded from its file, appended as new ones are
/// computed. Not shared between threads (one per row being read).
class RowReadBacks {
public:
  /// An empty dir means no cache (get() finds nothing, put() keeps nothing).
  RowReadBacks(const std::string &dir, const std::string &rowPD, int layers,
               const std::string &buildDigest);

  const CachedReadBack *get(const std::string &witnessKey) const;
  /// Records one read-back, in memory and (buffered) for the file.
  void put(const std::string &witnessKey, const CachedReadBack &r);
  /// Appends what put() recorded since the last flush, under the file's lock.
  void flush();

  size_t hits() const { return hits_; }
  size_t loaded() const { return entries_.size(); }

private:
  std::string path_, digest_;
  bool rewrite_ = false; // the file's digest was stale (or absent): start it afresh
  std::unordered_map<std::string, CachedReadBack> entries_;
  std::string pending_;
  mutable size_t hits_ = 0;
};

} // namespace outgoing

#endif

//
//  outgoingnamer.h
//
//  Naming outgoing links from their diagrams, in the search.
//

/*! \file utils/surfer/cobound/outgoing/outgoingnamer.h
 *  \brief Names a surface's outgoing curves by drawing them
 *  (knotbuilder::DiagramDrawer) instead of drilling their complement.
 *
 *  The outgoing boundary of the search's thickening is a copy of
 *  knotbuilder's triangulation T of the incoming link (OutgoingMap, thickening.h), so
 *  its curves draw straight into a diagram: microseconds, where drilling,
 *  simplifying and naming a complement -- with a Pachner search behind
 *  a census miss -- took tens of milliseconds and was nearly all of the
 *  drain. OutgoingNamer draws, and linknaming's DiagramNamer names the drawing
 *  with proof (linknamer.h), falling back to the complement route for what
 *  a drawing cannot name.
 */

#ifndef SURFER_OUTGOINGNAMER_H
#define SURFER_OUTGOINGNAMER_H

#include <map>
#include <memory>
#include <mutex>
#include <optional>
#include <string>
#include <unordered_map>

#include <triangulation/dim3.h>

#include "diagramtriangulation/thickening/thickening.h"
#include "linknaming/linknamer.h"
#include "cobound/outgoing/outgoinglink.h"
#include "diagramtriangulation/todiagram.h"
#include "linknaming/complement/linkcomplement.h"
#include "surfer/enumeration/surfacesearch.h"

namespace outgoing {

/**
 * Names every boundary curve by its complement (census::nameComplement(), the
 * census route) -- the search's namer where a search has no OutgoingNamer, and
 * OutgoingNamer's own fallback.
 */
class ComplementNamer : public ComplementBoundaryNamer {
  public:
    ComplementNamer();
};

/**
 * The search's BoundaryNamer: the outgoing boundary's curves by their
 * drawing (name()), every other boundary component's, and each curve of a
 * multi-curve component on its own, by the complement route as
 * ComplementNamer names them. One per search; thread-safe.
 */
class OutgoingNamer : public BoundaryNamer {
  public:
    /**
     * \param knotT knotbuilder's triangulation of the incoming link, unmodified.
     * \param crossings the incoming diagram's crossing count.
     * \param cob the incoming link's thickening, after its last thicken() and not coned.
     * \param table outlives this namer.
     */
    OutgoingNamer(const regina::Triangulation<3> &knotT, size_t crossings,
                 const CobordismBuilder<3> &cob, const linknaming::SignatureTable &table);

    /** Whether boundary component `bc` is the outgoing one, which this
     *  namer draws. */
    bool handles(size_t bc) const { return bc == map_.boundaryComponent(); }
    /** All the curves of the outgoing boundary component, named together
     *  from their drawing (DiagramNamer). */
    std::string name(const Link &curves) const;

    std::string nameLink(size_t bc, const Link &curves) const override;
    std::string nameCurve(size_t bc, const Knot &curve) const override;
    const linknaming::NamingStats &stats() const { return diagramNamer_.stats(); }

    /**
     * Turns on orientedName() (outgoing_names; formerly --exact-far-side-names).
     * \param tables outlives this namer.
     * \param caches what naming learns about `tables`, shared with other
     *        namers over them (a goal run's searches share one); new when null.
     */
    void enableOrientedNames(const linknaming::Tables &tables,
                          std::shared_ptr<linknaming::TableCaches> caches = nullptr);
    bool orientedNamesOn() const { return linkNamer_ != nullptr; }

    /**
     * The name (linknaming/) of ONE surface's outgoing curves,
     * oriented as a cobordism from the incoming link: each curve reversed when
     * `flips` says its surface component runs against the incoming link
     * (outgoing::incomingFlips()). name() cannot give this -- it is asked once
     * per edge set, and an edge set's orientation depends on the surface --
     * so the search asks it per cobordism, for deduplicating by the ORIENTED
     * outgoing link (two surfaces whose outgoing links are different orientation
     * variants of one link are two cobordisms). Only the link namer's fast path
     * runs here (diagram matches and visible cuts, no Reidemeister search);
     * a piece it cannot match is written as its exact oriented signature,
     * which `cobound name` refines offline from the stored pair signature.
     * Cached by the drawn diagram's exact signature. nullopt when names
     * are off or the drawing fails.
     */
    std::optional<std::string> orientedName(
        const std::vector<OrientedCurve> &outgoing,
        const std::map<const regina::Edge<3> *, size_t> &surfaceOf,
        const std::map<size_t, int> &flips) const;

  private:
    /** This namer's drawing of `curves`, as namer_ asks for it. */
    linknaming::DrawnCurves draw(const Link &curves) const;

    ComplementNamer complement_; ///< everything name() does not draw
    OutgoingMap map_;
    knotbuilder::DiagramDrawer drawer_;
    linknaming::DiagramNamer diagramNamer_;
    std::unique_ptr<linknaming::LinkNamer> linkNamer_;
    mutable std::mutex orientedMutex_;
    mutable std::unordered_map<std::string, std::string> orientedCache_;
    /**< drawn diagram's exact signature -> its name. */
};

} // namespace outgoing

#endif

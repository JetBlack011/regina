//
//  outgoingnamer.h
//
//  Naming outgoing links from their diagrams, in the search.
//

/*! \file utils/surfer/cobound/outgoing/outgoingnamer.h
 *  \brief Names a surface's outgoing curves by drawing them
 *  (diagramtriangulation::DiagramDrawer) and naming the drawing with the one
 *  link namer (linknaming::LinkNamer::nameDrawing()).
 *
 *  The outgoing boundary of the search's thickening is a copy of
 *  the diagram's triangulation T of the incoming link (OutgoingMap, thickening.h), so
 *  its curves draw straight into a diagram: microseconds, where drilling,
 *  simplifying and naming a complement -- with a Pachner search behind
 *  a census miss -- took tens of milliseconds and was nearly all of the
 *  drain. What a drawing cannot name goes to the complement route, which the
 *  link namer runs on the curves' edges (census::nameComplement()).
 */

#ifndef SURFER_OUTGOINGNAMER_H
#define SURFER_OUTGOINGNAMER_H

#include <map>
#include <memory>
#include <optional>
#include <string>

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
 * census route): OutgoingNamer's for every boundary component it does not
 * draw, and for each curve of a multi-curve component on its own.
 */
class ComplementNamer : public ComplementBoundaryNamer {
  public:
    ComplementNamer();
};

/**
 * The search's BoundaryNamer: the outgoing boundary's curves drawn and named
 * by the link namer (name(), and per surface orientedName()); every other
 * boundary component's, and each curve of a multi-curve component on its
 * own, by the complement route as ComplementNamer names them. One per search;
 * thread-safe.
 */
class OutgoingNamer : public BoundaryNamer {
  public:
    /**
     * \param knotT the diagram's triangulation T of the incoming link, unmodified.
     * \param crossings the incoming diagram's crossing count.
     * \param cob the incoming link's thickening, after its last thicken() and not coned.
     * \param tables outlives this namer.
     * \param caches what naming learns about `tables`, shared with other
     *        namers over them (a goal run's searches share one); new when null.
     */
    OutgoingNamer(const regina::Triangulation<3> &knotT, size_t crossings,
                  const CobordismBuilder<3> &cob, const linknaming::Tables &tables,
                  std::shared_ptr<linknaming::TableCaches> caches = nullptr);

    /**
     * The link namer's limits in a search: its diagram steps and isometry, no
     * Reidemeister search (forward or from the table's side: on an outgoing
     * link that is none of the shortlisted entries it runs to the end, a
     * minute or more per name -- a 50k-surface search spent 57,000
     * thread-seconds there, 2026-09-29; `cobound name` refines such names
     * offline from the pair signature). The attempts each route made before
     * the routes were merged are kept: a knot one simplify() (no extra tries),
     * a link three (two extra).
     */
    static linknaming::NamerLimits limits();

    /** Whether boundary component `bc` is the outgoing one, which this
     *  namer draws. */
    bool handles(size_t bc) const { return bc == map_.boundaryComponent(); }
    /** All the curves of the outgoing boundary component, named together
     *  from their drawing, each curve oriented as the search carries it. */
    std::string name(const Link &curves) const;

    std::string nameLink(size_t bc, const Link &curves) const override;
    std::string nameCurve(size_t bc, const Knot &curve) const override;
    const linknaming::NamingStats &stats() const { return stats_; }

    /**
     * The name -- or description -- of ONE surface's outgoing curves,
     * oriented as a cobordism from the incoming link: each curve reversed
     * when `flips` says its surface component runs against the incoming link
     * (outgoing::incomingFlips()). name() cannot give this -- it is asked once
     * per edge set, and an edge set's orientation depends on the surface --
     * so the search asks it per cobordism, for deduplicating by the ORIENTED
     * outgoing link (two surfaces whose outgoing links are different
     * orientation variants of one link are two cobordisms). Never empty: a
     * drawing that fails, or a surface whose orientation cannot be read, is
     * named by the complement route, which only describes a link
     * (`complement:`, names.h).
     */
    std::string orientedName(const std::vector<OrientedCurve> &outgoing,
                             const std::map<const regina::Edge<3> *, size_t> &surfaceOf,
                             const std::map<size_t, int> &flips) const;

  private:
    /** This namer's drawing of `cycles`, as LinkNamer::nameDrawing() asks for it. */
    linknaming::DrawnCurves draw(const std::vector<diagramtriangulation::EdgeCycle> &cycles) const;

    ComplementNamer complement_; ///< everything name() does not draw
    OutgoingMap map_;
    diagramtriangulation::DiagramDrawer drawer_;
    linknaming::LinkNamer namer_;
    mutable linknaming::NamingStats stats_;
};

} // namespace outgoing

#endif

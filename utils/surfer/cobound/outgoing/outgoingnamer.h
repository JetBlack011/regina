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
 *  knotbuilder's triangulation T of the row (OutgoingMap, thickening.h), so
 *  its curves draw straight into a diagram: microseconds, where drilling,
 *  simplifying and recognising a complement -- with a Pachner search behind
 *  a census miss -- took tens of milliseconds and was nearly all of the
 *  drain. DiagramNamer draws, and linknaming's LinkNamer names the drawing
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

namespace farside {

/**
 * The search's BoundaryNamer for the outgoing boundary. One per row; name()
 * is thread-safe.
 */
class DiagramNamer : public BoundaryNamer {
  public:
    /**
     * \param knotT knotbuilder's triangulation of the row, unmodified.
     * \param crossings the row's crossing count.
     * \param cob the row's thickening, after its last thicken() and not coned.
     * \param table outlives this namer.
     */
    DiagramNamer(const regina::Triangulation<3> &knotT, size_t crossings,
                 const CobordismBuilder<3> &cob, const SignatureTable &table);

    bool handles(size_t bc) const override { return bc == map_.boundaryComponent(); }
    std::string name(const Link &curves) const override;
    const NamingStats &stats() const { return namer_.stats(); }

    /**
     * Turns on orientedName() (verifyslicegenus --exact-far-side-names).
     * \param tables outlives this namer.
     * \param caches what naming learns about `tables`, shared with other
     *        namers over them (a cascade's hops share one); new when null.
     */
    void enableExactNames(const exactnaming::ExactTables &tables,
                          std::shared_ptr<exactnaming::TableCaches> caches = nullptr);
    bool exactNamesOn() const { return exact_ != nullptr; }

    /**
     * The exact name (exactnaming/) of ONE surface's outgoing curves,
     * oriented as a cobordism from the row: each curve reversed when
     * `flips` says its surface component runs against the row
     * (farside::incomingFlips()). name() cannot give this -- it is asked once
     * per edge set, and an edge set's orientation depends on the surface --
     * so the search asks it per witness, for deduplicating by the ORIENTED
     * far side (two surfaces whose far sides are different orientation
     * variants of one link are two witnesses). Only exactnaming's fast path
     * runs here (diagram matches and visible cuts, no Reidemeister search);
     * a piece it cannot match is written as its exact oriented signature,
     * which farsidename refines offline from the stored pair signature.
     * Cached by the drawn diagram's exact signature. nullopt when exact names
     * are off or the drawing fails.
     */
    std::optional<std::string> orientedName(
        const std::vector<OrientedCurve> &outgoing,
        const std::map<const regina::Edge<3> *, size_t> &surfaceOf,
        const std::map<size_t, int> &flips) const;

  private:
    /** This namer's drawing of `curves`, as namer_ asks for it. */
    DrawnCurves draw(const Link &curves) const;

    OutgoingMap map_;
    knotbuilder::DiagramDrawer drawer_;
    LinkNamer namer_;
    std::unique_ptr<exactnaming::ExactNamer> exact_;
    mutable std::mutex exactMutex_;
    mutable std::unordered_map<std::string, std::string> exactCache_;
    /**< drawn diagram's exact signature -> its exact name. */
};

} // namespace farside

#endif

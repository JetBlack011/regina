//
//  farsidenaming.h
//
//  Naming far sides from their diagrams, in the search.
//

/*! \file utils/surfer/farsidenaming.h
 *  \brief Names a surface's outgoing curves by drawing them
 *  (knotbuilder::DiagramDrawer) instead of drilling their complement.
 *
 *  The outgoing boundary of the search's thickening is a copy of
 *  knotbuilder's triangulation T of the row (farsidecurves.h), so its curves
 *  draw straight into a diagram: microseconds, where drilling, simplifying
 *  and recognising a complement -- with a Pachner search behind a census
 *  miss -- took tens of milliseconds and was nearly all of the drain.
 *  DiagramNamer draws, lets Regina's Link::simplify() reduce the diagram,
 *  and names it with proof:
 *
 *    - no crossings left: "Unknot", or "<n>-component unlink" (a diagram
 *      without crossings is the unlink);
 *    - a knot whose simplified diagram is exactly a table knot's diagram
 *      (knotSig, mirror and reversal allowed -- a slice genus sees neither):
 *      that table name;
 *    - a knot seen before under this diagram: the name proved then;
 *    - a link whose diagram is exactly one of the link table's (any of its
 *      orientation variants, so up to orientation): that link's base name.
 *      A link name bears no bound in the solvers, as before;
 *    - any other link that is provably not an unlink -- a nonzero linking
 *      number, or a Jones polynomial other than the unlink's:
 *      "diagram:<signature>", a name that bears nothing but tells distinct
 *      links apart exactly (the unlink is the one link name that would
 *      bear a bound, so that is what must be ruled out first).
 *
 *  Everything else -- a knot the table does not know, a link the Jones
 *  polynomial cannot tell from an unlink, a degenerate drawing, a drawing
 *  the drawer refuses as not planar (a drawer defect, counted in
 *  NamingStats::nonPlanar) -- falls back to identify::identify(), the
 *  complement route, exactly as before.
 *  Whatever it returns is remembered against the diagram's signature (a
 *  diagram determines its link), so each distinct diagram costs at most one
 *  fallback per row, and repeats of it get the same name.
 */

#ifndef SURFER_FARSIDENAMING_H
#define SURFER_FARSIDENAMING_H

#include <atomic>
#include <map>
#include <memory>
#include <mutex>
#include <optional>
#include <string>
#include <unordered_map>
#include <unordered_set>

#include <triangulation/dim3.h>

#include "cobordismbuilder.h"
#include "exactnaming/exactnamer.h"
#include "farsidecurves.h"
#include "knotbuilder/diagramdrawer.h"
#include "linkcomplement.h"
#include "surfacesearch.h"

namespace farside {

/**
 * Diagram signatures of the knot and link tables: exact diagrams, so a hit
 * is a proof. Built once per process from the tables' PD codes.
 */
class SignatureTable {
  public:
    /**
     * \param knotTable, linkTable CSV files of "Name,PD,..." rows (either
     *        may be empty to skip it). Link names are oriented
     *        ("L6a3{1}"); the table maps every variant to its base name.
     */
    static SignatureTable fromTables(const std::string &knotTable,
                                     const std::string &linkTable);

    const std::string *knot(const std::string &knotSig) const;
    const std::string *link(const std::string &linkSig) const;
    bool isKnotName(const std::string &name) const { return knotNames_.contains(name); }
    size_t knots() const { return knots_.size(); }
    size_t links() const { return links_.size(); }

  private:
    std::unordered_map<std::string, std::string> knots_, links_;
    std::unordered_set<std::string> knotNames_;
};

/** How DiagramNamer named what it was asked to (cumulative). */
struct NamingStats {
    std::atomic<long long> calls{0}, unknots{0}, unlinks{0}, tableKnots{0},
        learnedKnots{0}, tableLinks{0}, diagramLinks{0}, jonesLinks{0},
        learnedLinks{0}, fallbacks{0}, learned{0}, microsDiagram{0},
        microsFallback{0};
    std::atomic<long long> nonPlanar{0};
    /**< Drawings the drawer refused as not planar (knotbuilder::NonPlanar):
         each is a drawer defect, named by the complement route instead. */
    std::atomic<long long> exactNamed{0}, exactCacheHits{0}, exactFailed{0};
    /**< orientedName(): names computed, answered from the cache, and
         drawings that failed (the witness then keeps its unoriented name). */
};

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
    const NamingStats &stats() const { return stats_; }

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
    std::string nameOnce(const Link &curves) const;

    OutgoingMap map_;
    knotbuilder::DiagramDrawer drawer_;
    const SignatureTable &table_;
    mutable std::mutex learnedMutex_;
    mutable std::unordered_map<std::string, std::string> learned_;
    /**< "K" + knotSig or "L" + unoriented link signature -> the name the
         complement route gave that diagram. */
    mutable std::unordered_map<size_t, regina::Laurent<regina::Integer>> unlinkJones_;
    /**< n -> the Jones polynomial of the n-component unlink. */
    std::unique_ptr<exactnaming::ExactNamer> exact_;
    mutable std::mutex exactMutex_;
    mutable std::unordered_map<std::string, std::string> exactCache_;
    /**< drawn diagram's exact signature -> its exact name. */
    mutable NamingStats stats_;
};

} // namespace farside

#endif

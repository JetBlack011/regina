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
 *  polynomial cannot tell from an unlink, a degenerate drawing -- falls
 *  back to identify::identify(), the complement route, exactly as before.
 *  Whatever it returns is remembered against the diagram's signature (a
 *  diagram determines its link), so each distinct diagram costs at most one
 *  fallback per row, and repeats of it get the same name.
 */

#ifndef SURFER_FARSIDENAMING_H
#define SURFER_FARSIDENAMING_H

#include <atomic>
#include <mutex>
#include <string>
#include <unordered_map>
#include <unordered_set>

#include <triangulation/dim3.h>

#include "cobordismbuilder.h"
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
    mutable NamingStats stats_;
};

} // namespace farside

#endif

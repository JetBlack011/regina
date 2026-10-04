//
//  exactnamer.h
//
//  Exact names for oriented far-side diagrams.
//

/*! \file utils/surfer/linknaming/linknamer.h
 *  \brief Names an oriented link diagram -- a far side drawn with the
 *  orientation its surface induces -- with a proof, as a table entry or a
 *  split union / connected sum of table entries, in the atlas's syntax.
 *
 *  1. Simplify (Regina's simplify() and simplifyExhaustive(), which never
 *     reflect or reverse), then cut the diagram into split pieces and
 *     visible connected summands (gaussdiagram.h), recursively. Components
 *     that end up in no piece are split unknots.
 *  2. Name each prime-looking piece:
 *     - exactly: its diagram, after further simplify() tries, IS one version
 *       (as written, mirrored, reversed, both) of a table entry
 *       (exacttables.h) -- this pins the oriented variant, and the mirror
 *       and reversal relative to the drawing;
 *     - by isometry, for hyperbolic pieces: among the table entries sharing
 *       the piece's HOMFLY polynomial (a hint), one whose complement is
 *       isometric to the piece's by an isometry carrying meridians to
 *       meridians (snappeaisometry.h: the SnapPea kernel's own test, whose
 *       positive answer is a combinatorial isomorphism, hence exact). Filling
 *       along the meridians, the piece is that link up to mirror and
 *       orientations; its variant is then pinned as for a search;
 *     - by search: a bounded Reidemeister search (regina::Link::rewrite(),
 *       which may reflect and reverse) reaches a diagram of a table base, so
 *       the piece is that link up to mirror and orientations; then every
 *       orientation variant of the base, and its mirror, whose HOMFLY
 *       polynomial or linking numbers differ from the piece's (computed on
 *       the drawing's own orientation) is ruled out. What survives is proved
 *       to contain the piece;
 *     - from the table's side, when that search finds nothing: table entries
 *       whose HOMFLY polynomial (in either mirror) equals the piece's, and
 *       whose table diagram is no bigger, are candidates -- a hint only.
 *       The proof is a diagram: for an alternating piece, its planar graph
 *       is in the flype orbit of a candidate's alternating table diagram
 *       (a flype is an isotopy, and an alternating diagram is its graph up
 *       to mirror); otherwise a bounded rewrite() outward FROM the
 *       candidate's table diagram reaches the piece in some orientation of
 *       its components. Then the variant is pinned as above. By the
 *       Menasco-Thistlethwaite flyping theorem, every reduced alternating
 *       diagram of a prime alternating table link is in the orbit, so this
 *       names the minimal diagrams the forward search cannot reach; the
 *       table-side rewrite covers non-alternating and non-minimal ones;
 *     - otherwise the piece is untabulated: `diagram:<sig>`, the signature
 *       of its oriented diagram (reversal disallowed, the lesser of the
 *       piece and its global reverse), which determines the link.
 *  3. Compose the name (the atlas grammar): `A u B` split factors (Unknot
 *     last), `A#mB` composite knots, `K #_? L` knots summed into a link,
 *     `#{A[?] # B[?] ; ...}` links summed along components. `?` marks a
 *     component index this namer does not compute.
 *
 *  A name is EXACT when it is an identity: it determines the far side up to
 *  mirror and global reversal, neither of which a slice genus sees, so it
 *  may stand as a node of the cobordism graph and receive bounds. That is:
 *  the unknot and unlinks; one tabulated or untabulated piece whose variant
 *  is pinned, possibly with split unknots; and composite knots whose every
 *  relative mirror and reversal that matters is pinned. Every other name is
 *  a proved DESCRIPTION: each piece is exactly the named table link, but how
 *  the pieces are joined is not all recorded, which is enough for rules
 *  that bound its slice genus from its pieces' but not for an identity.
 */

#ifndef SURFER_EXACTNAMING_EXACTNAMER_H
#define SURFER_EXACTNAMING_EXACTNAMER_H

#include <atomic>
#include <functional>
#include <map>
#include <memory>
#include <mutex>
#include <optional>
#include <set>
#include <string>
#include <vector>

#include <link/link.h>

#include "linknaming/tables.h"
#include "linknaming/diagrams/gaussdiagram.h"
#include "linknaming/isometry/isometry.h"

namespace linknaming {

struct NamerLimits {
    /** Name a piece by its own diagram when it is a version of a table
     *  entry. Always on in use; validation turns it off to force the other
     *  steps (tests/isometry_validation.cpp). */
    bool exactDiagram = true;
    /** Hyperbolic pieces: an isometry of complements carrying meridians to
     *  meridians (snappeaisometry.h), against a HOMFLY shortlist, tried
     *  before any Reidemeister search. */
    bool isometry = true;
    /** Check every variant the isometry pins against the HOMFLY polynomial
     *  and linking numbers; a disagreement is a bug and throws. */
    bool checkIsometryPins = true;
    int simplifyTries = 24;         /**< extra simplify() tries per piece */
    int exhaustiveHeight = 1;       /**< simplifyExhaustive() height (0: none) */
    int searchHeight = 2;           /**< rewrite() height (-1: no search) */
    size_t searchVisits = 20000;    /**< rewrite() diagrams per piece */
    size_t maxSearchCrossings = 16; /**< no search above this many crossings */
    /** A second, deeper search for small pieces the first misses (the old
     *  diagram pipeline's second round): height, visits, crossings. */
    int deepHeight = 3;
    size_t deepVisits = 200000;
    size_t maxDeepCrossings = 12;
    /** From the table's side (see step 2): a HOMFLY shortlist, then the
     *  flype orbit (alternating) and a rewrite() outward from each
     *  candidate's table diagram, heights 1..tableSideHeight in turn
     *  across candidates (-1: none). A table diagram of 8 crossings needed
     *  height 4 to reach a 10-crossing drawing of it (~30 s). */
    int tableSideHeight = 4;
    size_t tableSideVisits = 3000000; /**< rewrite() diagrams per candidate and height */
    size_t maxTableSideCrossings = 12; /**< no table-side step above this */
};

/** How one prime-looking piece was named. */
struct PieceName {
    /** How the name was proved: the diagram itself; an isometry of
     *  complements carrying meridians to meridians, which also pins the
     *  variant (and a knot's mirror and reversal) from its action on the
     *  oriented meridians; a Reidemeister search (forward, flype orbit or
     *  from the table's side), which proves the base, the variant then
     *  pinned by invariants. */
    enum class By { exactDiagram, isometry, searchAndInvariants, untabulated };
    By by = By::untabulated;
    size_t components = 0;
    size_t crossings = 0;           /**< of its simplified diagram */
    std::vector<size_t> origin;     /**< far-side component of each of its components */
    std::string base;               /**< table base, when tabulated */
    /** Canonical names (exacttables.h) the piece may be; exactly one when
     *  its oriented variant is pinned. */
    std::vector<std::string> names;
    std::optional<bool> mirror;     /**< relative to the drawing, when pinned */
    std::optional<bool> reverse;
    std::string sig;                /**< untabulated: its oriented signature */

    bool pinned() const { return by == By::untabulated || names.size() == 1; }
    /** The name alone: the table name or `diagram:<sig>`; alternatives joined by '|'. */
    std::string display() const;
};

/** A far side's name, with how it was proved. */
struct LinkName {
    std::string name;
    bool exact = false;             /**< an identity; see the file comment */
    bool pinned = false;            /**< every piece's variant pinned */
    size_t factors = 0;             /**< split factors, split unknots included */
    size_t splitUnknots = 0;
    std::vector<PieceName> pieces;
    /** e.g. "diagram", "search+homfly", "untabulated", joined per piece. */
    std::string proof() const;
};

/**
 * What an ExactNamer learns about its tables, whatever its limits: the
 * entries' HOMFLY polynomials and the index over them, their flype orbits,
 * their complements in the SnapPea kernel, and their canonical names. Each
 * is built on first use -- the HOMFLY index alone is ~34,000 polynomials --
 * so namers over the same tables can share one (ExactNamer::caches()), as a
 * cascade's node namer and every hop's far-side namer do. Guarded by its
 * own locks.
 */
struct TableCaches {
    explicit TableCaches(const ExactTables &t) : tables(&t) {}
    const ExactTables *tables; /**< the tables these entries belong to */

    std::mutex cacheMutex; /**< homfly and flypeOrbits */
    std::map<std::pair<const TableEntry *, bool>, regina::Laurent2<regina::Integer>> homfly;
    /** HOMFLY polynomial (either mirror, as a string) -> the bases having it;
     *  built on first use of the table-side step. */
    std::once_flag homflyIndexOnce;
    std::map<std::string, std::vector<std::string>> homflyIndex;
    /** Flype orbits by entry, built on first use. */
    std::map<const TableEntry *, std::set<std::string>> flypeOrbits;
    /** Table entries' complements in the SnapPea kernel, built on first use
     *  (guarded by their own mutex: building one takes the kernel's). */
    std::mutex kernelMutex;
    std::map<const TableEntry *, std::unique_ptr<KernelLink>> kernelLinks;
    /** canonicalName(), by entry, filled a whole base at a time. */
    std::mutex classMutex;
    std::map<const TableEntry *, std::string> canonicalOf;
    std::vector<std::string> classConflicts;
};

class ExactNamer {
  public:
    /** \param caches shared with other namers over the same `tables`
     *         (caches()); a new one when null. */
    explicit ExactNamer(const ExactTables &tables, NamerLimits limits = {},
                        std::shared_ptr<TableCaches> caches = nullptr);

    /** This namer's table caches, to share with another namer over the same
     *  tables. */
    const std::shared_ptr<TableCaches> &caches() const { return caches_; }

    /**
     * \param drawn an oriented, planar diagram of the far side, component i
     *        the i-th far-side curve (knotbuilder::Diagram::link()).
     */
    LinkName name(const regina::Link &drawn) const;

    /** Piece identification alone (exposed for tests). */
    PieceName identify(const GaussDiagram &piece) const;

    /** An entry's canonical name: one per class of variants of its base that
     *  are the same oriented link up to mirror and global reversal, whether
     *  their table diagrams coincide (ExactTables::canonical()) or an
     *  isometry carries meridians to meridians with a uniform orientation
     *  sign. Built per base on first use. */
    const std::string &canonicalName(const TableEntry &e) const;

    /** Cuts a diagram at its visible connected-sum spheres (after
     *  simplifying it) into prime pieces, each keeping its components'
     *  origins; a diagram with no visible sphere gives itself. Exposed for
     *  the cascade, which turns the pieces into summand nodes. */
    void decompose(const GaussDiagram &g, std::vector<GaussDiagram> &primes) const;

  private:
    const regina::Laurent2<regina::Integer> &homfly(const TableEntry &e, bool mirror) const;
    /** Table entries whose HOMFLY polynomial (either mirror) is `h` and whose
     *  table diagram (minimal) is no bigger than `l`, smallest first; the
     *  index is built on first use. A shortlist only: it proves nothing. */
    std::vector<const TableEntry *> homflyCandidates(
        const regina::Link &l, const regina::Laurent2<regina::Integer> &h) const;
    /** What the isometry step of identify() proves: the base of a
     *  shortlisted table entry whose complement is isometric to `l`'s by an
     *  isometry carrying meridians to meridians, and the canonical names of
     *  the variants some such isometry reaches with a uniform orientation
     *  sign (the same oriented link up to mirror and global reversal), with
     *  the mirror and reversal relative to the named entry when pinned. */
    struct IsometryMatch {
        std::string base;
        std::set<std::string> names;
        std::optional<bool> mirror, reverse;
    };
    std::optional<IsometryMatch> isometryMatch(const regina::Link &l,
                                               const regina::Laurent2<regina::Integer> &h) const;
    /** The variants of `base` (canonical names) whose HOMFLY polynomial and
     *  linking numbers, in one mirror or the other, are the piece's; the
     *  mirrors that survive go to `mirrors`. */
    std::set<std::string> invariantSurvivors(const GaussDiagram &piece,
                                             const regina::Laurent2<regina::Integer> &h,
                                             const std::string &base, std::set<bool> *mirrors) const;
    /** The kernel's complement of an entry's table diagram, built on first use. */
    const KernelLink &kernelLinkOf(const TableEntry &e) const;
    /** The table-side step of identify(): the base of a table entry proved
     *  to be the link of `l` up to mirror and orientations, if any. */
    std::optional<std::string> tableSideBase(const regina::Link &l,
                                             const regina::Laurent2<regina::Integer> &h) const;
    /** Whether an alternating graph (canonicalPlantri, reflection allowed) is
     *  in the flype orbit of an entry's table diagram's graph. */
    bool inFlypeOrbit(const TableEntry &e, const std::string &graph) const;

    const ExactTables &tables_;
    NamerLimits limits_;
    std::shared_ptr<TableCaches> caches_; /**< never null */

  public:
    /** Variants merged by isometry whose literature 4-genera differ: a
     *  table error or a bug, and never silently. */
    std::vector<std::string> classConflicts() const;
};

/** Pairwise linking numbers of a diagram's components, sorted. */
std::vector<long> linkingNumbers(const GaussDiagram &g);

} // namespace linknaming

/** The curves of one boundary component, as linkcomplement.h holds them
 *  (edges of a triangulation), for the complement route. */
class Link;

namespace linknaming {

/** How LinkNamer (and DiagramNamer's exact names) named what it was asked
 *  to (cumulative). */
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
    std::atomic<long long> microsExact{0};
    /**< Time in the exact namer itself (computed names only, not cache hits). */

    /** The slowest single naming so far: which route, what it named, how
        long. One slow name can hold a whole drain's last thread. */
    void noteDuration(long long micros, const char *route, const std::string &name);
    long long slowestMicros() const { return slowestMicros_.load(); }
    std::string slowest() const; ///< "<route> <name>", or empty

    /** The `diagram naming:` body, as verifyslicegenus and each cascade hop
        print it: counts by outcome, times by route, and the slowest name. */
    std::string summary() const;

  private:
    std::atomic<long long> slowestMicros_{0};
    mutable std::mutex slowestMutex_;
    std::string slowest_;
};

/** A drawing of one boundary component's curves, as LinkNamer::name()
 *  asks for it. */
struct DrawnCurves {
    enum class Outcome {
        drawn,     /**< `diagram` is the drawing. */
        nonPlanar, /**< The drawer refused it as not planar: a drawer defect. */
        failed     /**< Degenerate, or not a drawable set of curves. */
    };
    Outcome outcome = Outcome::failed;
    regina::Link diagram; /**< When drawn: component i is the i-th curve. */
    bool someLinking = false;
    /**< When drawn: whether some pair of curves has a nonzero linking number. */
};

/**
 * Names the curves of one boundary component from their drawing, with
 * proof: lets Regina's Link::simplify() reduce the diagram, and then
 *
 *   - no crossings left: "Unknot", or "<n>-component unlink" (a diagram
 *     without crossings is the unlink);
 *   - a knot whose simplified diagram is exactly a table knot's diagram
 *     (knotSig, mirror and reversal allowed -- a slice genus sees neither):
 *     that table name;
 *   - a knot seen before under this diagram: the name proved then;
 *   - a link whose diagram is exactly one of the link table's (any of its
 *     orientation variants, so up to orientation): that link's base name.
 *     A link name bears no bound in the solvers;
 *   - any other link that is provably not an unlink -- a nonzero linking
 *     number, or a Jones polynomial other than the unlink's:
 *     "diagram:<signature>", a name that bears nothing but tells distinct
 *     links apart exactly (the unlink is the one link name that would bear
 *     a bound, so that is what must be ruled out first).
 *
 * Everything else -- a knot the table does not know, a link the Jones
 * polynomial cannot tell from an unlink, a drawing that failed or that the
 * drawer refused as not planar (counted in NamingStats::nonPlanar) -- falls
 * back to census::identify(), the complement route. Whatever it returns is
 * remembered against the diagram's signature (a diagram determines its
 * link), so each distinct diagram costs at most one fallback, and repeats
 * of it get the same name. Names are perturbed as identify()'s are under
 * census::perturbNamesForTesting.
 */
class LinkNamer {
  public:
    /** \param table outlives this namer. */
    explicit LinkNamer(const SignatureTable &table);

    /**
     * The name of `curves` -- all the curves of one boundary component --
     * from the drawing `draw` makes of them. `draw` reports a drawing that
     * failed rather than throwing. Thread-safe.
     */
    std::string name(const Link &curves,
                     const std::function<DrawnCurves()> &draw) const;

    /** What this namer has named (cumulative). A caller naming by other
     *  routes too (DiagramNamer::orientedName()) adds its counts here. */
    NamingStats &stats() const { return stats_; }

  private:
    std::string nameOnce(const Link &curves,
                         const std::function<DrawnCurves()> &draw) const;

    const SignatureTable &table_;
    mutable std::mutex learnedMutex_;
    mutable std::unordered_map<std::string, std::string> learned_;
    /**< "K" + knotSig or "L" + unoriented link signature -> the name the
         complement route gave that diagram. */
    mutable std::unordered_map<size_t, regina::Laurent<regina::Integer>> unlinkJones_;
    /**< n -> the Jones polynomial of the n-component unlink. */
    mutable NamingStats stats_;
};

} // namespace linknaming

#endif

//
//  linknamer.h
//
//  THE namer: names, with proof, for oriented link diagrams, and the route a
//  search names its outgoing links by.
//

/*! \file utils/surfer/linknaming/linknamer.h
 *  \brief Names an oriented link diagram -- an outgoing link drawn with the
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
 *       (tables.h) -- this pins the oriented variant, and the mirror
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
 *  A NAME is an identity: it determines the outgoing link up to mirror and
 *  global reversal, neither of which a slice genus sees, so it may stand as
 *  a link of the cobordism graph and receive bounds (LinkName::isName). That is:
 *  the unknot and unlinks; one tabulated or untabulated piece whose variant
 *  is pinned, possibly with split unknots; and composite knots whose every
 *  relative mirror and reversal that matters is pinned. Anything weaker is
 *  a proved DESCRIPTION: each piece is exactly the named table link, but how
 *  the pieces are joined is not all recorded (`?`), which is enough for rules
 *  that bound its slice genus from its pieces' but not for an identity; or a
 *  link of more than one component named only by its complement
 *  (`complement:`, names.h), which bears nothing.
 *
 *  A search names its outgoing links by the full route, nameDrawing(): a
 *  drawing (todiagram.h) by the steps above, under the search's limits;
 *  a knot no table piece names, and every drawing that failed, by its
 *  complement (census::nameComplement()), which names a knot (Gordon-Luecke)
 *  and only describes a link. See nameDrawing().
 */

#ifndef SURFER_LINKNAMING_LINKNAMER_H
#define SURFER_LINKNAMING_LINKNAMER_H

#include <atomic>
#include <functional>
#include <map>
#include <memory>
#include <mutex>
#include <optional>
#include <set>
#include <string>
#include <unordered_map>
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
    /** Extra simplify() tries for a one-component piece; -1: simplifyTries.
     *  A search's namer keeps the attempts each of its routes made before
     *  they were merged (outgoing::OutgoingNamer::limits()). */
    int knotSimplifyTries = -1;
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
    std::vector<size_t> origin;     /**< outgoing component of each of its components */
    std::string base;               /**< table base, when tabulated */
    /** Canonical names (tables.h) the piece may be; exactly one when
     *  its oriented variant is pinned. */
    std::vector<std::string> names;
    std::optional<bool> mirror;     /**< relative to the drawing, when pinned */
    std::optional<bool> reverse;
    std::string sig;                /**< untabulated: its oriented signature */

    bool pinned() const { return by == By::untabulated || names.size() == 1; }
    /** The name alone: the table name or `diagram:<sig>`; alternatives joined by '|'. */
    std::string display() const;
};

/** A link's name, with how it was proved. */
struct LinkName {
    std::string name;
    bool isName = false;             /**< an identity; see the file comment */
    bool pinned = false;            /**< every piece's variant pinned */
    size_t factors = 0;             /**< split factors, split unknots included */
    size_t splitUnknots = 0;
    std::vector<PieceName> pieces;
    /** e.g. "diagram", "search+homfly", "untabulated", joined per piece. */
    std::string proof() const;
};

/** A drawing of the curves of one boundary component, as LinkNamer::nameDrawing()
 *  names it. */
struct DrawnCurves {
    enum class Outcome {
        drawn,     /**< `diagram` is the drawing. */
        nonPlanar, /**< The drawer refused it as not planar: a drawer defect. */
        failed     /**< Degenerate, or not a drawable set of curves. */
    };
    Outcome outcome = Outcome::failed;
    regina::Link diagram; /**< When drawn: component i is the i-th curve. */
};

/** How LinkNamer::nameDrawing() named a drawing (NamingStats counts each). */
enum class NamedBy {
    unknot,      /**< one curve, no crossing left */
    unlink,      /**< several curves, no crossing left; or proved by the complement */
    tableKnot,   /**< a knot of table pieces: its table signature, a piece's
                      diagram or isometry, a sum of them */
    learnedKnot, /**< a knot whose simplified diagram the complement named before */
    tableLink,   /**< a link of table pieces only */
    otherLink,   /**< a link with an untabulated piece (`diagram:`) */
    learnedLink, /**< a link the complement proved an unlink before */
    complement,  /**< the complement route: an untabulated knot, a failed drawing */
};

/** What LinkNamer::nameDrawing() made of one drawing. */
struct RoutedName {
    std::string name; /**< a name, or a description (`?`, `|`, `complement:`) */
    bool isName = false; /**< an identity (see the file comment) */
    NamedBy by = NamedBy::complement;
    bool cached = false;           /**< this drawing was named before */
    bool complementCalled = false; /**< the complement route ran for this call */
    bool nonPlanar = false;        /**< the drawer refused the drawing as not planar */
    bool learned = false;          /**< a complement answer was kept against its diagram */
    long long microsDiagram = 0;   /**< time in the diagram route */
    long long microsComplement = 0; /**< time in the complement route */
};

/**
 * What a LinkNamer learns about its tables, whatever its limits: the
 * entries' HOMFLY polynomials and the index over them, their flype orbits,
 * their complements in the SnapPea kernel, and their canonical names. Each
 * is built on first use -- the HOMFLY index alone is ~34,000 polynomials --
 * so namers over the same tables can share one (LinkNamer::caches()), as a
 * goal run's link namer and every search's outgoing namer do. Guarded by its
 * own locks.
 */
struct TableCaches {
    explicit TableCaches(const Tables &t) : tables(&t) {}
    const Tables *tables; /**< the tables these entries belong to */

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

class LinkNamer {
  public:
    /** \param caches shared with other namers over the same `tables`
     *         (caches()); a new one when null. */
    explicit LinkNamer(const Tables &tables, NamerLimits limits = {},
                        std::shared_ptr<TableCaches> caches = nullptr);

    /** This namer's table caches, to share with another namer over the same
     *  tables. */
    const std::shared_ptr<TableCaches> &caches() const { return caches_; }

    /**
     * \param drawn an oriented, planar diagram of the link, component i its
     *        i-th curve (diagramtriangulation::Diagram::link()).
     * \param simplified `drawn` was simplified already (once): it is not
     *        simplified again before it is cut into pieces.
     */
    LinkName name(const regina::Link &drawn, bool simplified = false) const;

    /** The complement route: the drawn curves named by their complement
     *  (census::nameComplement() of their edges). */
    using ComplementRoute = std::function<std::string()>;

    /**
     * THE route a search names an outgoing link by (plan, phase 7.1): one
     * namer for every drawing, `components` curves drawn as `drawing`.
     *
     *   - A drawing that failed (degenerate, or refused as not planar) goes
     *     to the complement route: a knot is named by its complement
     *     (Gordon-Luecke); a link only described, `complement:<name>`
     *     (names.h), unless the complement proves it an unlink.
     *   - A knot: simplified once; no crossing left, the unknot; its table
     *     signature (the simplified diagram IS a version of a table knot);
     *     a complement answer learned for that diagram; its pieces (cut at
     *     visible connected-sum spheres; each by its diagram, with
     *     NamerLimits::knotSimplifyTries more attempts, its isometry, its
     *     search, as the limits allow): a table knot or a sum of them. A
     *     knot some piece of which no table names is named by its complement,
     *     the answer learned against its simplified diagram.
     *   - A link: name() under this namer's limits (its pieces, oriented as
     *     drawn). One whose untabulated piece might be an unlink in disguise
     *     -- every linking number 0 and the unlink's Jones polynomial -- is
     *     tried on the complement route, which only ever replaces the name by
     *     a proved unlink.
     *
     * Every drawing's answer is kept against its exact signature, so a drawing
     * seen again is answered the same, and costs one signature. Names and
     * descriptions alike are perturbed under census::perturbNamesForTesting
     * (SURFER_TEST_PERTURB_NAMES). Thread-safe.
     */
    RoutedName nameDrawing(const DrawnCurves &drawing, size_t components,
                           const ComplementRoute &complement) const;

    /** Naming one piece alone (exposed for tests). */
    PieceName namePiece(const GaussDiagram &piece) const;

    /** An entry's canonical name: one per class of variants of its base that
     *  are the same oriented link up to mirror and global reversal, whether
     *  their table diagrams coincide (Tables::canonical()) or an
     *  isometry carries meridians to meridians with a uniform orientation
     *  sign. Built per base on first use. */
    const std::string &canonicalName(const TableEntry &e) const;

    /** Cuts a diagram at its visible connected-sum spheres (after
     *  simplifying it, unless `simplified`) into prime pieces, each keeping
     *  its components' origins; a diagram with no visible sphere gives
     *  itself. Exposed for goal runs, which turn the pieces into summand
     *  links. */
    void decompose(const GaussDiagram &g, std::vector<GaussDiagram> &primes,
                   bool simplified = false) const;

  private:
    const regina::Laurent2<regina::Integer> &homfly(const TableEntry &e, bool mirror) const;
    /** Table entries whose HOMFLY polynomial (either mirror) is `h` and whose
     *  table diagram (minimal) is no bigger than `l`, smallest first; the
     *  index is built on first use. A shortlist only: it proves nothing. */
    std::vector<const TableEntry *> homflyCandidates(
        const regina::Link &l, const regina::Laurent2<regina::Integer> &h) const;
    /** What the isometry step of namePiece() proves: the base of a
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
    /** The table-side step of namePiece(): the base of a table entry proved
     *  to be the link of `l` up to mirror and orientations, if any. */
    std::optional<std::string> tableSideBase(const regina::Link &l,
                                             const regina::Laurent2<regina::Integer> &h) const;
    /** Whether an alternating graph (canonicalPlantri, reflection allowed) is
     *  in the flype orbit of an entry's table diagram's graph. */
    bool inFlypeOrbit(const TableEntry &e, const std::string &graph) const;
    /** nameDrawing() before the perturbation, and without the drawing cache. */
    RoutedName nameDrawingUncached(const DrawnCurves &drawing, size_t components,
                                   const ComplementRoute &complement) const;

    const Tables &tables_;
    NamerLimits limits_;
    std::shared_ptr<TableCaches> caches_; /**< never null */

    /** What nameDrawing() keeps: complement answers learned per simplified
     *  diagram, every drawing's answer, the unlinks' Jones polynomials. */
    struct RouteMemory {
        std::mutex mutex;
        std::unordered_map<std::string, RoutedName> learned;
        /**< "K" + knotSig, or "L" + the unoriented signature -> the complement's answer. */
        std::unordered_map<std::string, RoutedName> drawings;
        /**< a drawing's exact signature (Link::sig<2>(false, false, true)) -> its answer. */
        std::unordered_map<size_t, regina::Laurent<regina::Integer>> unlinkJones;
        /**< n -> the Jones polynomial of the n-component unlink. */
    };
    /** Drawings remembered before the memory is cleared and starts again. */
    static constexpr size_t kDrawingMemoryLimit = 1'000'000;
    std::unique_ptr<RouteMemory> memory_; /**< never null */

  public:
    /** Variants merged by isometry whose literature 4-genera differ: a
     *  table error or a bug, and never silently. */
    std::vector<std::string> classConflicts() const;
};

/** Pairwise linking numbers of a diagram's components, sorted. */
std::vector<long> linkingNumbers(const GaussDiagram &g);

/** What a search's outgoing namer named, by route (cumulative): nameDrawing()'s
 *  answers, counted per entry point (outgoing::OutgoingNamer). */
struct NamingStats {
    /** Per edge set ("far sides drawn"): every answer, by how it was named. A
     *  drawing named before counts as it was named then, but a complement
     *  answer reused counts as learned. */
    std::atomic<long long> calls{0}, unknots{0}, unlinks{0}, tableKnots{0},
        learnedKnots{0}, tableLinks{0}, otherLinks{0}, learnedLinks{0};
    /** The complement route, both entry points: how often it ran, answers kept
     *  for their diagrams, and drawings refused as not planar (each a drawer
     *  defect: diagramtriangulation::NonPlanar). */
    std::atomic<long long> fallbacks{0}, learned{0}, nonPlanar{0};
    /** Per surface: oriented names computed, answered from the drawing memory,
     *  and the complement-route answers among them (drawings that failed). */
    std::atomic<long long> orientedNamed{0}, orientedCacheHits{0}, orientedByComplement{0};
    std::atomic<long long> microsDiagram{0}, microsFallback{0}, microsOriented{0};
    /**< Time in the diagram route per edge set, in the complement route, and in
         the diagram route per surface. */

    /** Counts one answer: `oriented` for a surface's oriented name, else one
     *  edge set's (with its component count). */
    void count(const RoutedName &r, bool oriented, size_t components);

    /** The slowest single naming so far: which route, what it named, how
        long. One slow name can hold a whole drain's last thread. */
    void noteDuration(long long micros, const char *route, const std::string &name);
    long long slowestMicros() const { return slowestMicros_.load(); }
    std::string slowest() const; ///< "<route> <name>", or empty

    /** The `diagram naming:` body, as every search prints it: counts by
        outcome, times by route, and the slowest name. */
    std::string summary() const;

  private:
    std::atomic<long long> slowestMicros_{0};
    mutable std::mutex slowestMutex_;
    std::string slowest_;
};

} // namespace linknaming

#endif

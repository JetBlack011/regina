//
//  exactnamer.h
//
//  Exact names for oriented far-side diagrams.
//

/*! \file utils/surfer/exactnaming/exactnamer.h
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

#include <map>
#include <memory>
#include <mutex>
#include <optional>
#include <set>
#include <string>
#include <vector>

#include <link/link.h>

#include "exacttables.h"
#include "gaussdiagram.h"
#include "snappeaisometry.h"

namespace exactnaming {

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
struct FarSideName {
    std::string name;
    bool exact = false;             /**< an identity; see the file comment */
    bool pinned = false;            /**< every piece's variant pinned */
    size_t factors = 0;             /**< split factors, split unknots included */
    size_t splitUnknots = 0;
    std::vector<PieceName> pieces;
    /** e.g. "diagram", "search+homfly", "untabulated", joined per piece. */
    std::string proof() const;
};

class ExactNamer {
  public:
    explicit ExactNamer(const ExactTables &tables, NamerLimits limits = {})
        : tables_(tables), limits_(limits) {}

    /**
     * \param drawn an oriented, planar diagram of the far side, component i
     *        the i-th far-side curve (knotbuilder::Diagram::link()).
     */
    FarSideName name(const regina::Link &drawn) const;

    /** Piece identification alone (exposed for tests). */
    PieceName identify(const GaussDiagram &piece) const;

    /** An entry's canonical name: one per class of variants of its base that
     *  are the same oriented link up to mirror and global reversal, whether
     *  their table diagrams coincide (ExactTables::canonical()) or an
     *  isometry carries meridians to meridians with a uniform orientation
     *  sign. Built per base on first use. */
    const std::string &canonicalName(const TableEntry &e) const;

  private:
    void decompose(const GaussDiagram &g, std::vector<GaussDiagram> &primes) const;
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
    mutable std::mutex cacheMutex_;
    mutable std::map<std::pair<const TableEntry *, bool>, regina::Laurent2<regina::Integer>> homfly_;
    /** HOMFLY polynomial (either mirror, as a string) -> the bases having it;
     *  built on first use of the table-side step. */
    mutable std::once_flag homflyIndexOnce_;
    mutable std::map<std::string, std::vector<std::string>> homflyIndex_;
    /** Flype orbits by entry, built on first use. */
    mutable std::map<const TableEntry *, std::set<std::string>> flypeOrbits_;
    /** Table entries' complements in the SnapPea kernel, built on first use
     *  (guarded by their own mutex: building one takes the kernel's). */
    mutable std::mutex kernelCacheMutex_;
    mutable std::map<const TableEntry *, std::unique_ptr<KernelLink>> kernelLinks_;
    /** canonicalName(), by entry, filled a whole base at a time. */
    mutable std::mutex classMutex_;
    mutable std::map<const TableEntry *, std::string> canonicalOf_;

  public:
    /** Variants merged by isometry whose literature 4-genera differ: a
     *  table error or a bug, and never silently. */
    std::vector<std::string> classConflicts() const;

  private:
    mutable std::vector<std::string> classConflicts_;
};

/** Pairwise linking numbers of a diagram's components, sorted. */
std::vector<long> linkingNumbers(const GaussDiagram &g);

} // namespace exactnaming

#endif

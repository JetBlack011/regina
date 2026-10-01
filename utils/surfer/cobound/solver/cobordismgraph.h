//
//  cobordismgraph.h
//
//  Created by John Teague on 08/02/2026.
//

#ifndef COBORDISMGRAPH_H

#define COBORDISMGRAPH_H

#include <climits>
#include <functional>
#include <map>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

#include "surfer/enumeration/surfacesearch.h"

/*! \file utils/surfer/cobound/solver/cobordismgraph.h
 *  \brief The name/genus resolution graph verifyslicegenus.cpp builds up
 *  across its --input rows: given a set of witnessed surfaces and cobordisms
 *  between named knots/links, works out what each row's slice genus must be.
 *
 *  \section cg_math The inequality this is all built on
 *
 *  Throughout, `g_4(L)` is the slice genus of a link in the following sense:
 * the minimum genus over **connected** orientable properly embedded surfaces in
 * `B^4` bounded by `L`. Let `Sigma_1` be a connected genus-`g` cobordism in
 * `S^3 x I` from `L_0` (with `n_0` components) to `L_1` (with `n_1`
 * components), and let `Sigma_2` be a connected genus-`h` surface in `B^4`
 * bounded by `L_1`. Gluing them along `L_1` gives a connected surface bounded
 * by `L_0` with
 *  \f[
 *    \chi = \chi(\Sigma_1) + \chi(\Sigma_2)
 *         = (2 - 2g - n_0 - n_1) + (2 - 2h - n_1),
 *  \f]
 *  and matching that against `chi = 2 - 2G - n_0` gives `G = g + h + n_1 - 1`.
 *  Hence
 *  \f[
 *    g_4(L_0) \le g_4(L_1) + g + n_1 - 1, \qquad
 *    g_4(L_0) \ge g_4(L_1) - g - n_0 + 1,
 *  \f]
 *  the second being the first with the roles of `L_0` and `L_1` swapped.
 *
 *  For knots, both collapse to the familiar `|g_4(K_0) - g_4(K_1)| <= g`
 * exactly when `n_0 = n_1 = 1`. Note that the component-count terms are not
 * optional for links: a genus-0 cobordism from a knot `K` to a 2-component
 * unlink bounds `g_4(K)` by `0 + 0 + 2 - 1 = 1`, not by 0, so dropping the
 * correction would "prove" `K` slice on wrong evidence.
 *
 *  \section cg_farside Which far sides may carry a bound at all
 *
 *  identify::identify() names a far side by its COMPLEMENT. Whether that
 *  name may feed either inequality above depends entirely on how many
 *  components the far side has, and the rule is decided by farSideBearsBound():
 *
 *  - **One component.** Gordon-Luecke: a knot is determined by its
 *    complement up to mirroring, and `g_4` is mirror-invariant. The name is
 *    the knot, and it bounds.
 *  - **An `n`-component unlink.** identify() emits that name only from a
 *    structural proof (free pi_1 => split unlink), which forces the link
 *    itself, not just its exterior. It bounds, with no component penalty
 *    (see propagate()).
 *  - **Anything else with two or more components: it bounds NOTHING.** One
 *    link complement belongs to infinitely many non-isotopic links (Rolfsen
 *    twisting along an unknotted component), most of them not orientation
 *    variants of each other and with different slice genera. Our own
 *    peripheral tests found one census name standing for 18 different links
 *    with `g_4` ranging over 0..2. So "the complement matched L6a3" does not
 *    mean "this is some orientation of L6a3", and a max/min over L6a3's
 *    oriented variants is not a bound over the actual possibilities. Such a
 *    far side is still RECORDED (it is the honest observation, and the input
 *    a later peripheral resolution needs) but the solver skips it.
 *
 *  Orientation is the lesser half of the same problem: `L4a1{0}` (slice genus
 *  0) and `L4a1{1}` (slice genus 1) are the same manifold, which is why
 *  NameTable::candidates() still widens a base name to its oriented variants.
 *  That widening is correct as far as it goes; it is simply not sufficient
 *  on its own, and so it is never the thing that licenses a bound.
 */

namespace cobordismgraph {

/** One --input row: a name to verify/bound the slice genus of, plus its
 * literature bounds and the PD code to search from. */
struct InputRow {
    std::string name;
    std::string pdNotation;
    int lo = 0, hi = 0; // literature genus bounds; lo == hi except for a
                        // handful of [lo;hi] range rows
    int crossings = 0;
};

/** Sentinels for "no bound derived yet" */
constexpr int NO_UPPER_BOUND = INT_MAX;
constexpr int NO_LOWER_BOUND = INT_MIN;

/* Names */

/**
 * The number of components of whatever `name` names, worked out from the
 * name itself.
 *
 * \note This is a fallback for names we have no better information about.
 * When a component count is *observed* (the number of boundary curves a
 * found surface actually put on an ambient boundary component), prefer the
 * observation: it is a fact about this surface, whereas this is an
 * inference from a string.
 */
int componentsFromName(const std::string &name);

/** `name` with any trailing orientation tag removed: `L6a3{0}` -> `L6a3`. */
std::string baseName(const std::string &name);

/**
 * A COMPOSITE far-side name, "K #_c L": the knot K connect-summed into
 * component c of the link L (L a base name, c in L's PD component order).
 *
 * The component is part of the name because K # L is only well defined once
 * it is given. It does not enter the bound: g_4(K #_c L) <= g_4(K) + g_4(L)
 * for every c, by boundary-connect-summing the two minimal surfaces.
 * \a knot keeps its chirality ("m3_1"); g_4 is mirror-invariant, so the
 * solver looks it up without the "m".
 */
struct CompositeName {
    std::string knot;
    int component = 0;  /**< -1 for "?": the namer did not compute it */
    std::string link;   /**< a base, or a tagged variant ("L7n1{1}") from an exact name */
};

/** The parts of a composite name, or nullopt if `name` is not one. */
std::optional<CompositeName> compositeParts(const std::string &name);

/**
 * The summands of a composite KNOT name, "A#B#..." (each "K", "mK" or
 * "Unknot"), or an empty vector if `name` is not one. A census name
 * ("m129 : #2"), a knot-into-link sum ("3_1 #_0 L2a1") and a split name
 * ("A u B") are all refused.
 */
std::vector<std::string> knotSummands(const std::string &name);

/**
 * A (possibly marked) prime knot name without its marks: "mr8_17" -> "8_17".
 * An exact name marks a summand "m" (mirrored) and/or "r" (reversed)
 * relative to the table's diagram, where its symmetry type makes the mark
 * matter; g_4 sees neither. Any other name is returned unchanged.
 */
std::string stripKnotMarks(const std::string &name);

/**
 * The pieces of a sum along components, each with its component count:
 * "K1 #_? K2 #_? L" (knots summed into components of the link L), or
 * "#{A[?] # B[?] ; ...}" (sum sites, one per shared component). nullopt for
 * anything else, or when a piece's name does not state its component count.
 *
 * In "#{...}" a piece summed at two sites is written at both, so every
 * occurrence counts as a piece of its own. That can only over-count, which
 * weakens both bounds built from the pieces and never makes one unsound.
 */
std::optional<std::vector<std::pair<std::string, int>>> sumPieces(const std::string &name);

/**
 * What a far side named EXACTLY (--far-side-exact) may be: the name itself,
 * or -- where the namer proved it one of a few orientation variants it could
 * not tell apart -- those alternatives "A|B". Never widened to a base's
 * variants. A split or a sum is one candidate: its alternatives live inside
 * its factors and pieces.
 */
std::vector<std::string> exactCandidates(const std::string &name);

/**
 * A prime knot's symmetry type, as data/knot_symmetry.csv records it
 * (KnotInfo, cross-checked against SnapPy through 10 crossings). It fixes the
 * concordance inverse -K = m(K^r):
 *   reversible            -K = mK
 *   fullyAmphicheiral     -K = K = mK
 *   negativeAmphicheiral  -K = K   (and K^r = mK, so mK is its own inverse too)
 *   positiveAmphicheiral  -K = K^r    } not expressible without a reversal
 *   chiral                -K = m(K^r) } marker, which our names lack
 */
enum class SymmetryType {
    chiral,
    reversible,
    positiveAmphicheiral,
    negativeAmphicheiral,
    fullyAmphicheiral
};

/** Parses KnotInfo's spelling ("negative amphicheiral"), or nullopt. */
std::optional<SymmetryType> parseSymmetryType(const std::string &text);

/**
 * An identify() result reduced to the name the graph should key on.
 *
 * Strips the trailing parenthetical if there is one; leaves everything else
 * (`"Unknot"`, `"3-component unlink"`, a bare isoSig, an undecorated census
 * name) untouched.
 */
std::string normalizeIdentifiedName(const std::string &name);

/** What we know about one name independently of any surface we have found. */
struct NameInfo {
    int components = 1;
    bool haveLiterature = false;
    int litLo = 0;
    int litHi = 0;
};

/**
 * Populated from the --input tables. A links run and a knots run pointed at
 * the same table set see the same NameTable, which is what lets a knot row's
 * search use a link far side and vice versa.
 */
class NameTable {
  public:
    /** Registers one literature row. Safe to call twice for the same name. */
    void addLiterature(const std::string &name, int lo, int hi);

    /** Returns what's known about `name`, or nullptr if it was never seen. */
    const NameInfo *find(const std::string &name) const;

    /**
     * `name`'s component count: the registered value if we have one, else
     * componentsFromName().
     */
    int components(const std::string &name) const;

    /**
     * The oriented variants sharing `name`'s base.
     *
     * Returns every registered name sharing `name`'s base, filtered to those
     * with `observedComponents` components when that is given (a variant with
     * the wrong number of components simply isn't what was found). Falls back
     * to `{name}` when the base is unregistered: a knot name, `"Unknot"`, an
     * `"<n>-component unlink"`, or a bare isoSig, none of which have oriented
     * variants to disambiguate between.
     *
     * This is a description of a NAME, not a claim about a far side: for a
     * multi-component far side the true candidate set is not enumerable from
     * a complement at all (see \ref cg_farside), so callers must never treat
     * this list as licensing a bound. farSideBearsBound() decides that.
     */
    std::vector<std::string>
    candidates(const std::string &name,
               std::optional<int> observedComponents = std::nullopt) const;

    size_t size() const { return info_.size(); }

    /** Records a knot's symmetry type (--knot-symmetry). */
    void setSymmetry(const std::string &knot, SymmetryType type) {
        symmetry_[knot] = type;
    }

    /** A knot's symmetry type, or nullptr if unknown. */
    const SymmetryType *symmetry(const std::string &knot) const {
        auto it = symmetry_.find(knot);
        return it == symmetry_.end() ? nullptr : &it->second;
    }

    /**
     * Whether upperOf()/lowerOf() bound sums along components and splits
     * with link factors from their pieces (--sum-rules): the additive upper
     * bound, and the lower bounds g_4(F_i) - sum_{j != i} (g_4(F_j) + n(F_j)
     * - 1) for a split and g_4(P_i) - sum_{j != i} (g_4(P_j) + n(P_j) - 1)
     * for a sum. Off by default, like every other widening of what bounds.
     */
    void setSumRules(bool on) { sumRules_ = on; }
    bool sumRules() const { return sumRules_; }

  private:
    bool sumRules_ = false;
    std::unordered_map<std::string, SymmetryType> symmetry_;
    std::unordered_map<std::string, NameInfo> info_;
    std::unordered_map<std::string, std::vector<std::string>> byBase_;
};

/**
 * Whether the knot `name` is slice by ELEMENTARY concordance-group reasoning,
 * making it an anchor exactly like the unknot.
 *
 * `K # -K` bounds a ribbon disc for every K, where `-K = m(K^r)` is the
 * concordance inverse, and a sum of slice knots is slice. So a composite knot
 * is elementarily slice when its summands pair off into inverse pairs. The
 * inverse depends on the summand's symmetry: for an invertible K, `-K = mK`;
 * if K is also amphicheiral, `-K = K`. Summands without a certified symmetry
 * (NameTable::symmetry), and every NON-invertible summand, are refused: our
 * names record chirality but not reversal, so for a non-invertible K the name
 * cannot say whether its neighbour is -K or its reverse (8_17 is the trap).
 *
 * The two long-standing anchors "3_1#m3_1" and "4_1#4_1" are accepted even
 * with no symmetry data loaded, so a run without --knot-symmetry loses
 * nothing it had before.
 */
bool isElementarySlice(const std::string &name, const NameTable &names);

/* Witnesses */

/** Whether a witness bounds its subject on its own, or only relative to
 * another name. */
enum class WitnessKind { direct, cobordism };

/**
 * One surface we actually found.
 */
struct Witness {
    WitnessKind kind = WitnessKind::cobordism;

    std::string
        subject; /**< The name whose own boundary this surface realizes. */
    int subjectComponents = 1; /**< Observed curve count on `subject`'s side. */

    std::string other; /**< The far side's identified name; empty iff `kind ==
                          direct`. */
    std::vector<std::string> otherCandidates;
    /**< The oriented variants `other` could be (see NameTable::candidates).
         Empty iff `kind == direct`. */
    int otherComponents = 1; /**< Observed curve count on the far side. */

    int genus = 0;
    /**< The genus of the connected surface this witness provides. For a
         disconnected find, the tubed genus (SurfaceFoundInfo::tubedGenus),
         not the meaningless whole-complex figure. */
    bool tubed = false;
    /**< Whether the surface as found was disconnected and had to be tubed
         to give `genus`. Recorded so a reader can tell that the pairSig
         names a disconnected complex, with the tubing as the step from it
         to the surface the bound is about. */

    std::string pairSig;
    /**< The found surface's pair signature, if captured -- held in memory
         only until the witness is on disk (see fileOffset). */
    long long fileOffset = -1;
    /**< Byte offset of this witness's line in the witness file it was loaded
         from or appended to, or -1 while it exists only in memory. A loaded
         witness never holds its ~11 KB pair signature in memory; anything
         that needs it reads it back through this offset. */
    std::string pairSigKey;
    /**< witnesskey::witnessKey(pairSig), sha1(pairsig)[:12] -- the key the
         per-witness far-side resolution file is keyed on. Filled at load
         only when resolutions are in use, and for every new witness. */
    /**
     * Whether this witness's far side has been PROVED, per witness, to be
     * the link named in `other` -- by the peripheral system (an isometry or
     * curve-respecting isomorphism carrying our meridians), or by a split
     * decomposition whose every factor is a knot (Gordon-Luecke). Set only
     * on the solver-side copy by applyFarSideResolutions(); never recorded.
     *
     * This is what lets farSideBearsBound() accept a multi-component far
     * side: the complement alone never determines a link, but complement
     * plus meridians does, and then the oriented variants of `other` ARE a
     * complete candidate set, so the max/min over them is sound.
     */
    bool farSideProved = false;
    /**
     * Whether the far side was named EXACTLY (--far-side-exact): redrawn
     * from this witness's own pair signature, oriented by its surface, and
     * named with a proof that is an identity (one link up to mirror and
     * global reversal). Such a far side may also RECEIVE a bound from the
     * subject, whatever its component count. Solver-side only; never
     * recorded.
     */
    bool farSideExact = false;

    // Provenance: which search produced this, and under what budget.
    std::string sourceRow;
    int thickenLayers = 1;
    long long maxFaces = 0; /**< 0 means the search was unbounded. */

    int resolvedVertices = 0;
    /**< How many ambient vertices the found surface meets itself at. 0 for
         an embedded surface. Positive only under --resolve-unlinked, where
         each such vertex is an unlinked self-intersection, and a
         perturbation near those vertices turns the surface into an embedded
         one of the same topology and boundary (paper §4.5), so the witness
         bounds exactly as an embedded one would. The solver ignores it, but
         whether it is zero is part of witnessIdentity(): an embedded
         witness must never be discarded as a duplicate of one that rests on
         the resolution theorem. */
};

/** The dedup identity of a witness: kind, subject and its component count,
 * far-side name, far-side component count, genus, tubed, and whether it is
 * resolved. Only the first witness of each identity is recorded.
 *
 * For a far side that bears a bound (farSideBearsBound()) equal identity
 * means the witnesses bound exactly the same thing. For any other far side
 * the name comes from its complement, which several links can share, so two
 * surfaces reaching DIFFERENT links can collapse here. That is deliberate:
 * keying such far sides by their own edges instead (tried 2026-09-26) made
 * nearly every surface its own witness -- 19,690 witnesses from 23,672
 * surfaces for L9a43{1;0} at cap 3, against 109 by name -- which no store
 * can hold at a million surfaces per row, while the per-witness pipeline
 * has found such a collapse in 21 of 228,580 witnesses, and those far sides
 * bear nothing in the solvers. A cheap invariant of the far-side LINK
 * (per-curve knot types, linking numbers) would separate them boundedly. */
std::string witnessIdentity(const Witness &w);

/** Whether `witnesses` already contains a witness with `w`'s
 * witnessIdentity(). Linear; the search keeps a hashed set of identities
 * instead, and this remains for tests and small callers. */
bool haveWitness(const std::vector<Witness> &witnesses, const Witness &w);

/**
 * Whether `w`'s far side is allowed to supply or receive a slice-genus bound.
 *
 * True exactly when the far side has one observed component (a knot; sound by
 * Gordon-Luecke) or is a structurally recognised `"Unknot"` /
 * `"<n>-component unlink"`. False for every other multi-component far side,
 * whatever it is called: a Thistlethwaite name, a census name, a bare
 * isoSig -- none of these determines a link from a complement alone (see
 * \ref cg_farside). Gated on the OBSERVED component count rather than the
 * spelling of the name, so an alias that renames a two-component far side to
 * something knot-shaped cannot slip through.
 */
bool farSideBearsBound(const Witness &w);

/* Solving */

/** How a derived upper bound was arrived at. */
enum class Basis {
    constructive,
    /**< Obtained using the surfaces we've found (assuming unlinks are genus 0)
     */
    literatureAssisted,
    /**< Obtained by appealing to a known value. */
};

/** The bounds propagate() derived for one name. */
struct Bounds {
    int hi = NO_UPPER_BOUND; /**< Derived upper bound; see Basis. */

    std::vector<std::string> support;
    std::vector<std::string> lowerSupport;

    int lo = NO_LOWER_BOUND;
    Basis basis = Basis::constructive;

    // Provenance of whichever witness last improved `hi`.
    WitnessKind kind = WitnessKind::direct;
    std::string viaName;
    int viaGenus = 0;
    std::string pairSig;
    long long pairSigOffset = -1; /**< The witness's Witness::fileOffset, for
                                       reading pairSig back when it is not in
                                       memory. */
    bool tubed = false;

    bool haveUpper() const { return hi != NO_UPPER_BOUND; }
    bool haveLower() const { return lo != NO_LOWER_BOUND; }
};

/**
 * Derives every bound the witness set supports, by relaxing to a fixpoint.
 *
 * For each witness, with `g` its genus, `n_a` the subject's component count
 * and `n_b` the far side's:
 * \f[
 *   hi[a] \leftarrow \min\bigl(hi[a],\ \max_{c \in cand(b)} hi(c) + g + n_b -
 * 1\bigr),
 * \f]
 * \f[
 *   lo[a] \leftarrow \max\bigl(lo[a],\ \min_{c \in cand(b)} lo(c) - g - n_a +
 * 1\bigr),
 * \f]
 * and `hi[a] <- min(hi[a], g)` for a direct witness. Every edge is relaxed in
 * both directions, since a cobordism is symmetric.
 *
 * `hi(c)` prefers the derived bound and falls back to `c`'s literature upper
 * bound, marking the result Basis::literatureAssisted when it does. `lo(c)`
 * uses the literature lower bound, or `c`'s own derived one.
 *
 * Both relaxations are monotone (`hi` only falls, `lo` only rises) and clamped
 * to `[0, ...)`, since a genus is never negative; each propagation step
 * subtracts a non-negative `g + n - 1`, so no cycle can move a value upward
 * indefinitely.
 *
 * Seeded with the two facts that need no literature and no search: the
 * unknot bounds a disk, and an `n`-component unlink bounds a connected
 * planar surface (tube the `n` discs together -- genus 0, `n` boundary
 * circles), so both get `hi = lo = 0` constructively.
 */
/**
 * An upper bound proved outside this solver's witnesses: a certified
 * cascadesearch proof (regina-john cascade/, checked by cascade_check.py),
 * read from the atlas's data/cascade_proofs.csv (--cascade-proofs).
 * `support` names the literature values the proof's leaves use; empty for a
 * constructive proof. Seeded like a direct witness, so it grounds chains
 * exactly as one does.
 */
struct ExternalProof {
    std::string name;
    int genus = 0;
    std::vector<std::string> support;
    std::string source; ///< e.g. "cascade:2026-09-28_knots8/8_16"
};

std::unordered_map<std::string, Bounds>
propagate(const std::vector<Witness> &witnesses, const NameTable &names,
          const std::vector<ExternalProof> &external = {});

/* Reporting */

/** What a name's derived bounds amount to, judged against the literature. */
enum class Status {
    verified,
    /**< Verified using surfaces we found */
    verifiedAssisted,
    /**< Verified by appealing to known values */
    improved,
    /**< Our upper bound beats the literature's but hasn't reached its
         lower bound (yippee!!) */
    pinned,
    /**< Derived upper and lower bounds agree, pinning the value, without
         the literature having an exact value to match. */
    bounded, /**< Some derived bound exists, but it doesn't improve anything. */
    unresolved, /**< Nothing derived. */
    contradiction,
    /**< Derived bounds are inconsistent with the literature. Mathematically
       impossible, so it means a bug -- see verifyslicegenus.cpp's fatal-bug
       halt. */
};

/** One name's final verdict, for writing out and for printing. */
struct Verdict {
    Status status = Status::unresolved;
    Bounds bounds;
    bool haveLiterature = false;
    int litLo = 0, litHi = 0;
    int value = 0;      /**< The established genus, when status pins one. */
    std::string reason; /**< Human-readable one-liner, for a contradiction. */
};

/** Judges `bounds` for `name` against whatever literature `names` holds. */
Verdict judge(const std::string &name, const Bounds &bounds,
              const NameTable &names);

/** Walks `via` back through `bounds`' own provenance chain to whatever it
 * ultimately rests on (a direct witness, the Unknot, or an unlink), joining
 * the path with ';'. No cycles. */
std::string
buildDependsOn(const std::string &via,
               const std::unordered_map<std::string, Bounds> &bounds);

/* Boundary classification */

/**
 * One non-search-side ambient boundary component's identity, as classified
 * by splitBoundary().
 */
struct BoundarySide {
    std::string name;
    int components = 1;
};

/**
 * Splits a SurfaceBoundaryInfo::boundaryComponents grouping into "the
 * search side" (this row's own side) and every other ambient boundary
 * component the surface touches.
 */
struct BoundarySplit {
    size_t searchCurveCount = 0; // 0 if the search side has no boundary here
    bool searchSideRejected = false;
    /**< Set when `requiredSearchEdges` was given and the curves on
         searchSideBC are not exactly those edges -- the surface's boundary
         there is some other link, so it says nothing about this row. */
    bool unnamedSide = false;
    /**< Set when a non-search-side component carried no name at all, which
         describeBoundary_() never produces; the caller treats it as a bug. */
    std::vector<BoundarySide> otherSides;
};

/**
 * The search side is identified by geometry, never by name.
 *
 * Component `searchSideBC` is the search side. In a seeded search that is
 * all there is to it: the seed is L x {0} and no other triangle with an edge
 * in that boundary component is ever searchable, so its curves are L by
 * construction (verifyslicegenus asserts this once per row). An unseeded
 * search has no such guarantee, and passes `requiredSearchEdges` -- the
 * sorted boundary-triangulation edge indices of the row's own link -- so
 * that a surface whose search-side curves are any other edge set is
 * rejected. Identified names are deliberately not consulted: they are not
 * canonical (a census hit's "#N" varies from one identification to the
 * next), and comparing them once silently discarded whole rows.
 */
BoundarySplit
splitBoundary(const std::vector<BoundaryComponentNames> &boundaryComponents,
              size_t searchSideBC,
              const std::vector<size_t> *requiredSearchEdges = nullptr);

/* Orientation matching */

/**
 * A row's own PD-tagged diagram edges (knotbuilder::TriangulationWithLink's
 * `edges`/`reversed`), translated into directed pairs of *vertex indices*
 * within some triangulation combinatorially isomorphic to the row's own
 * knotbuilder triangulation -- normally the ambient search-side boundary
 * component, rebuilt via `BoundaryComponent<4>::build()` (see
 * buildRowOrientation()).
 */
struct RowOrientation {
    std::unordered_map<size_t, size_t> tailOf;
    /**< Edge index of L in the search-side triangulation -> index of that
         edge's tail vertex under the row's PD orientation. Keyed by edge,
         not by vertex pair, so two edges joining the same pair of vertices
         can never be confused. */
    std::vector<size_t> edges; /**< Sorted keys of tailOf: L's edge set. */
    std::unordered_map<size_t, size_t> rowIndexOf;
    /**< Edge index of L in the search-side triangulation -> that edge's
         position in the rowEdges given to buildRowOrientation(), i.e. which
         edge of the row's own link it is. Lets a caller tell which
         component of L a search-side curve is. */
    size_t components = 0;
    /**< How many closed curves `edges` forms, each checked to chain head to
         tail under the PD orientation (buildRowOrientation() throws
         otherwise). */
    bool divergedFromDefaultIsomorphism = false;
    /**< Whether the isomorphism isIsomorphicTo() would have returned maps L
         differently -- onto other edges, or with other directions. That is
         the map this code used to trust; see buildRowOrientation(). */
};

/**
 * Builds `rowEdges`/`rowReversed`'s RowOrientation against `searchSideTri`.
 *
 * The row's own triangulation and `searchSideTri` are isomorphic, but PD
 * triangulations have automorphisms, so "an isomorphism" does not determine
 * where L goes. When `requiredEdges` is given (a seeded search: the seed's
 * own edges in that boundary component, i.e. L x {0} exactly), the
 * isomorphism used is one taking L's edges onto exactly that set. Any such
 * choice is sound: two of them differ by an automorphism of the row's
 * triangulation preserving L, a PL homeomorphism of S^3 taking L to itself,
 * and the slice genus is invariant under homeomorphism and mirroring.
 *
 * \throws regina::InvalidArgument if `rowEdges` is empty, if the two
 * triangulations are not isomorphic, if no isomorphism takes L onto
 * `requiredEdges`, or if the image of L fails to chain into closed directed
 * curves.
 */
RowOrientation
buildRowOrientation(const std::vector<const regina::Edge<3> *> &rowEdges,
                    const std::vector<bool> &rowReversed,
                    const regina::Triangulation<3> &searchSideTri,
                    const std::vector<size_t> *requiredEdges = nullptr);

/** How a surface's search-side boundary compares with the row's orientation;
 * see classifyRowOrientation(). */
enum class OrientationVerdict {
    match,        /**< Some choice of orientation on each surface component
                       induces the row's own orientation. */
    mismatch,     /**< Some surface component's curves induce an orientation
                       pattern no flip of that component can fix: the surface
                       witnesses a different oriented variant of the link. */
    foreignEdge,  /**< A search-side edge is not one of the row's own. */
    incoherentCurve, /**< A single curve's edges disagree in direction, or a
                          curve's surface component is unknown. */
};

/**
 * Compares one found surface's induced boundary orientation on the
 * search side with `row`.
 *
 * `curves` come from KnottedSurface::orientedBoundaryLinks(), which orients
 * each CONNECTED COMPONENT of the surface independently and arbitrarily;
 * `surfaceComponentOf` (KnottedSurface::boundaryEdgeSurfaceComponent())
 * says which component each edge belongs to. So the curves are grouped by
 * surface component, and each group must match `row` all at once or be
 * reversed all at once; different groups are independent. That is exactly
 * the freedom the surface has: its components can be oriented separately
 * and still tubed into one oriented surface, because a split two-component
 * unlink bounds an oriented annulus for either relative orientation.
 *
 * `foreignEdge` and `incoherentCurve` cannot happen for a correctly built
 * row in a seeded search; the caller treats them as bugs.
 */
OrientationVerdict classifyRowOrientation(
    const RowOrientation &row, const std::vector<OrientedCurve> &curves,
    const std::map<const regina::Edge<3> *, size_t> &surfaceComponentOf);

} // namespace cobordismgraph

#endif // COBORDISMGRAPH_H

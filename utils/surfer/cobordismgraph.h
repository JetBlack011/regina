//
//  cobordismgraph.h
//
//  Created by John Teague on 08/02/2026.
//

#ifndef COBORDISMGRAPH_H

#define COBORDISMGRAPH_H

#include <climits>
#include <functional>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

#include "surfacesearch.h"

/*! \file utils/surfer/cobordismgraph.h
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

  private:
    std::unordered_map<std::string, NameInfo> info_;
    std::unordered_map<std::string, std::vector<std::string>> byBase_;
};

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

    std::string
        pairSig; /**< The found surface's pair signature, if captured. */

    // Provenance: which search produced this, and under what budget.
    std::string sourceRow;
    int thickenLayers = 1;
    long long maxFaces = 0; /**< 0 means the search was unbounded. */
};

/** Whether `witnesses` already contains an equivalent witness. This is the
 * dedup key that keeps a harvest run from capturing thousands of pair
 * signatures. */
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
std::unordered_map<std::string, Bounds>
propagate(const std::vector<Witness> &witnesses, const NameTable &names);

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
    std::vector<BoundarySide> otherSides;
};

BoundarySplit
splitBoundary(const std::vector<BoundaryComponentNames> &boundaryComponents,
              size_t searchSideBC, const std::string &rowOwnName);

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
    std::unordered_map<size_t, size_t>
        headOf; // tail vertex index -> head vertex index
};

/**
 * Builds `rowEdges`/`rowReversed`'s RowOrientation against `searchSideTri`.
 *
 * \throws regina::InvalidArgument if `rowEdges` is empty, or if its own
 * triangulation is not isomorphic to `searchSideTri` (should not happen
 * when `searchSideTri` is actually built from the ambient boundary
 * component the row's diagram was seeded into).
 */
RowOrientation
buildRowOrientation(const std::vector<const regina::Edge<3> *> &rowEdges,
                    const std::vector<bool> &rowReversed,
                    const regina::Triangulation<3> &searchSideTri);

/**
 * Whether `curves` (one found surface's own induced boundary curves on the
 * row's search-side ambient boundary component, from
 * KnottedSurface::orientedBoundaryLinks(), for one arbitrary choice of the
 * surface's two orientations) matches `row` -- either everywhere or
 * nowhere, across every curve (component). A mixed pattern (some curves
 * match, some don't) means the surface's own relative orientation between
 * components doesn't match this row's PD convention, and is rejected here
 * as a fully-mismatched pattern would be reversed by simply picking the
 * surface's other orientation.
 */
bool matchesRowOrientation(const RowOrientation &row,
                           const std::vector<OrientedCurve> &curves);

} // namespace cobordismgraph

#endif // COBORDISMGRAPH_H

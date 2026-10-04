//
//  solver.h
//
//  The atlas solver: bounds from cobordisms and the literature.
//

#ifndef SURFER_COBOUND_SOLVER_H
#define SURFER_COBOUND_SOLVER_H

#include <climits>
#include <string>
#include <unordered_map>
#include <vector>

#include "cobound/cobordisms/cobordism.h"
#include "cobound/solver/literature.h"

/*! \file utils/surfer/cobound/solver/solver.h
 *  \brief The name/genus resolution graph built up across the --input rows:
 *  given a set of cobordisms between named knots/links, works out what each
 *  row's slice genus must be.
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
 *  \section cg_outgoing Which outgoing links may carry a bound at all
 *
 *  census::nameComplement() names an outgoing link by its COMPLEMENT. Whether that
 *  name may feed either inequality above depends entirely on how many
 *  components the outgoing link has, and the rule is decided by outgoingBearsBound():
 *
 *  - **One component.** Gordon-Luecke: a knot is determined by its
 *    complement up to mirroring, and `g_4` is mirror-invariant. The name is
 *    the knot, and it bounds.
 *  - **An `n`-component unlink.** census::nameComplement() emits that name only from a
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
 *    outgoing link is still RECORDED (it is the honest observation, and the input
 *    a later peripheral resolution needs) but the solver skips it.
 *
 *  Orientation is the lesser half of the same problem: `L4a1{0}` (slice genus
 *  0) and `L4a1{1}` (slice genus 1) are the same manifold, which is why
 *  NameTable::candidates() still widens a base name to its oriented variants.
 *  That widening is correct as far as it goes; it is simply not sufficient
 *  on its own, and so it is never the thing that licenses a bound.
 */

namespace solver {

/** Sentinels for "no bound derived yet" */
constexpr int NO_UPPER_BOUND = INT_MAX;
constexpr int NO_LOWER_BOUND = INT_MIN;

/**
 * Whether `w`'s outgoing link is allowed to supply or receive a slice-genus bound.
 *
 * True exactly when the outgoing link has one observed component (a knot; sound by
 * Gordon-Luecke) or is a structurally proven `"Unknot"` /
 * `"<n>-component unlink"`. False for every other multi-component outgoing link,
 * whatever it is called: a Thistlethwaite name, a census name, a bare
 * isoSig -- none of these determines a link from a complement alone (see
 * \ref cg_outgoing). Gated on the OBSERVED component count rather than the
 * spelling of the name, so an alias that renames a two-component outgoing link to
 * something knot-shaped cannot slip through.
 */
bool outgoingBearsBound(const cobordisms::Cobordism &w);

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

    // Provenance of whichever cobordism last improved `hi`.
    cobordisms::CobordismKind kind = cobordisms::CobordismKind::direct;
    std::string viaName;
    int viaGenus = 0;
    std::string pairSig;
    long long pairSigOffset = -1; /**< The cobordism's Cobordism::fileOffset, for
                                       reading pairSig back when it is not in
                                       memory. */
    bool tubed = false;

    bool haveUpper() const { return hi != NO_UPPER_BOUND; }
    bool haveLower() const { return lo != NO_LOWER_BOUND; }
};

/**
 * Derives every bound the cobordism set supports, by relaxing to a fixpoint.
 *
 * For each cobordism, with `g` its genus, `n_a` the subject's component count
 * and `n_b` the outgoing link's:
 * \f[
 *   hi[a] \leftarrow \min\bigl(hi[a],\ \max_{c \in cand(b)} hi(c) + g + n_b -
 * 1\bigr),
 * \f]
 * \f[
 *   lo[a] \leftarrow \max\bigl(lo[a],\ \min_{c \in cand(b)} lo(c) - g - n_a +
 * 1\bigr),
 * \f]
 * and `hi[a] <- min(hi[a], g)` for a direct cobordism. Every cobordism is relaxed in
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
 * An upper bound proved outside this solver's cobordisms: a certified
 * goal run's certificate (cobound run with a goal, checked by cascade_check.py),
 * read from the atlas's data/cascade_proofs.csv (the certified_bounds key).
 * `support` names the literature values the proof's leaves use; empty for a
 * constructive proof. Seeded like a direct cobordism, so it grounds chains
 * exactly as one does.
 */
struct ExternalProof {
    std::string name;
    int genus = 0;
    std::vector<std::string> support;
    std::string source; ///< e.g. "cascade:2026-09-28_knots8/8_16"
};

std::unordered_map<std::string, Bounds>
propagate(const std::vector<cobordisms::Cobordism> &cobordisms, const NameTable &names,
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
       impossible, so it means a bug -- see driver/fatal.h's fatal-bug
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
 * ultimately rests on (a direct cobordism, the Unknot, or an unlink), joining
 * the path with ';'. No cycles. */
std::string
buildDependsOn(const std::string &via,
               const std::unordered_map<std::string, Bounds> &bounds);

} // namespace solver

#endif // SURFER_COBOUND_SOLVER_H

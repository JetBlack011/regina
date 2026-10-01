//
//  literature.h
//
//  What the tables say about a name, for the atlas solver.
//

#ifndef SURFER_COBOUND_LITERATURE_H
#define SURFER_COBOUND_LITERATURE_H

#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

/*! \file utils/surfer/cobound/solver/literature.h
 *  \brief What the knot and link tables say about a name, independently of
 *  any surface found: its component count, its literature 4-genus bounds,
 *  its oriented variants, a knot's symmetry type, and the composites that
 *  are slice by elementary concordance (the anchors).
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
} // namespace cobordismgraph

#endif // SURFER_COBOUND_LITERATURE_H

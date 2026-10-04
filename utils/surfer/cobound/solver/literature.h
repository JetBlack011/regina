//
//  literature.h
//
//  What the tables say about a name, for the atlas solver.
//

#ifndef SURFER_COBOUND_LITERATURE_H
#define SURFER_COBOUND_LITERATURE_H

#include <filesystem>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

#include "linknaming/tables.h"

/*! \file utils/surfer/cobound/solver/literature.h
 *  \brief What the knot and link tables say about a name, independently of
 *  any surface found: its component count, its literature 4-genus bounds,
 *  its oriented variants and a knot's symmetry type. The tables are read by
 *  linknaming/tables.h, which also holds the symmetry types and the anchors
 *  (linknaming::isElementarySlice()).
 */

namespace solver {

/** One --input row: a name to verify/bound the slice genus of, plus its
 * literature bounds and the PD code to search from. */
struct InputRow {
    std::string name;
    std::string pdNotation;
    int lo = 0, hi = 0; // literature genus bounds; lo == hi except for a
                        // handful of [lo;hi] range rows
    int crossings = 0;
};

/** A prime knot's symmetry type: the one enum, linknaming/tables.h's. */
using linknaming::SymmetryType;

/** What we know about one name independently of any surface we have found. */
struct NameInfo {
    int components = 1;
    bool haveLiterature = false;
    int litLo = 0;
    int litHi = 0;
};

/**
 * Populated from the --input tables. A links run and a knots run pointed at
 * the same table set see the same NameTable, which is what lets a knot's
 * search use an outgoing link named from the link table and vice versa.
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
     * This is a description of a NAME, not a claim about an outgoing link: for a
     * multi-component outgoing link the true candidate set is not enumerable from
     * a complement at all (see \ref cg_outgoing), so callers must never treat
     * this list as licensing a bound. outgoingBearsBound() decides that.
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

    /** Every recorded symmetry type, for linknaming::isElementarySlice(). */
    const linknaming::SymmetryTable &symmetries() const { return symmetry_; }

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
    linknaming::SymmetryTable symmetry_;
    std::unordered_map<std::string, NameInfo> info_;
    std::unordered_map<std::string, std::vector<std::string>> byBase_;
};

} // namespace solver

namespace solver {

/// Registers every row of a literature table in `names` (names and bounds
/// only; PD codes are skipped). Returns the number loaded.
size_t loadNameTable(const std::filesystem::path &path,
                     solver::NameTable &names);

/// The knot and link tables' names and literature bounds (the sign step's
/// candidate sets) and, with `knotSymmetry`, the knots' symmetry types (the
/// slice-composite anchors; their count into `symmetryTypes`): a goal run's
/// NameTable, and `sign`'s (without symmetry types).
solver::NameTable loadTableNames(const std::string &knotTable,
                                         const std::string &linkTable,
                                         const std::string &knotSymmetry,
                                         size_t *symmetryTypes = nullptr);

} // namespace solver

#endif // SURFER_COBOUND_LITERATURE_H

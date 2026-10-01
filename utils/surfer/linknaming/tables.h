//
//  exacttables.h
//
//  The knot and link tables as oriented diagrams, for naming far sides
//  exactly.
//

/*! \file utils/surfer/linknaming/tables.h
 *  \brief The knot and link tables indexed so that a diagram can be named
 *  with its orientation and chirality pinned.
 *
 *  Every table entry (a knot, or one oriented variant `L7n1{1}` of a link) is
 *  held as the diagram its PD code draws, in four versions: as written,
 *  mirrored, with every component reversed, and both. Each version is keyed
 *  by its Regina signature with neither reflection nor reversal allowed
 *  (rotation, an isotopy, allowed), so a diagram whose own such signature is
 *  a key IS that version of that entry, exactly. A second index, by the
 *  signature that allows both, gives the entry's base (the link up to
 *  orientation and mirror) for diagrams reached by a search that may reflect
 *  or reverse (regina::Link::rewrite()).
 *
 *  A slice genus sees neither a mirror nor a global reversal, so a table name
 *  denotes its link up to both; the versions matter only relative to each
 *  other, inside a sum or a split.
 */

#ifndef SURFER_EXACTNAMING_EXACTTABLES_H
#define SURFER_EXACTNAMING_EXACTTABLES_H

#include <filesystem>
#include <string>
#include <unordered_map>
#include <vector>

#include <link/link.h>

namespace exactnaming {

/** One table entry. */
struct TableEntry {
    std::string name;      /**< "3_1", "11n_34", "L7n1{1}" */
    std::string base;      /**< name without its orientation tag */
    size_t components = 0;
    std::string g4;        /**< the literature 4-genus as the table writes it: "1", "[0;1]" */
    regina::Link diagram;  /**< the table's PD code, as oriented */
};

/** A diagram that is exactly one version of a table entry. */
struct VersionMatch {
    const TableEntry *entry = nullptr;
    bool mirror = false;   /**< the entry reflected */
    bool reverse = false;  /**< every component of the entry reversed */
};

/** KnotInfo symmetry types (data/knot_symmetry.csv). */
enum class Symmetry { unknown, chiral, reversible, positiveAmphicheiral,
                      negativeAmphicheiral, fullyAmphicheiral };

class ExactTables {
  public:
    ExactTables() = default;
    ExactTables(ExactTables &&) noexcept = default;
    ExactTables &operator=(ExactTables &&) noexcept = default;
    ExactTables(const ExactTables &) = delete; // the indices point into entries_
    ExactTables &operator=(const ExactTables &) = delete;

    /**
     * \param knotTable, linkTable "Name,PD,Genus-4D" CSVs (the atlas tables).
     * \param knotSymmetry "name,symmetry_type,..." (may be empty: every knot
     *        then counts as chiral, which only ever withholds exactness).
     * \exception regina::InvalidArgument a table cannot be read, or two
     *        different bases share an unoriented diagram.
     */
    static ExactTables load(const std::string &knotTable, const std::string &linkTable,
                            const std::string &knotSymmetry);

    /** The versions whose exact signature (Link::sig<2>(false, false, true)) this is. */
    const std::vector<VersionMatch> *exact(const std::string &sig) const;
    /** The base whose unoriented signature (Link::sig<2>(true, true, true)) this is. */
    const std::string *base(const std::string &sig) const;
    /** Every entry of a base: a knot's one entry, or a link's oriented variants. */
    const std::vector<const TableEntry *> &variants(const std::string &base) const;
    const TableEntry *entry(const std::string &name) const;
    Symmetry symmetry(const std::string &knot) const;
    size_t size() const { return entries_.size(); }
    /** Every entry, in load order. */
    const std::vector<TableEntry> &entries() const { return entries_; }

    /**
     * The name that stands for this entry's link. Entries whose diagrams
     * coincide up to mirror and global reversal (the Hopf link's two
     * orientations, say, are each other's mirrors) are one link up to both,
     * which is all a table name denotes; each such class is named by its
     * least name.
     */
    const std::string &canonical(const std::string &name) const;
    /** Classes whose members' literature 4-genus disagree (a table error:
     *  one link cannot have two). Each as "name=g4 name=g4 ...". */
    const std::vector<std::string> &inconsistentClasses() const { return inconsistent_; }

  private:
    std::vector<TableEntry> entries_;
    std::unordered_map<std::string, size_t> byName_;
    std::unordered_map<std::string, std::string> canonical_;
    std::vector<std::string> inconsistent_;
    std::unordered_map<std::string, std::vector<VersionMatch>> exact_;
    std::unordered_map<std::string, std::string> unoriented_;
    std::unordered_map<std::string, std::vector<const TableEntry *>> byBase_;
    std::unordered_map<std::string, Symmetry> symmetry_;
};

/** A table PD code ("[[1;5;2;4];...]" or "PD[X[4; 1; 3; 2]; ...]") as Regina reads it. */
regina::Link linkFromTablePD(const std::string &pd);

} // namespace exactnaming

namespace witnessstore {

// ---- literature tables (Name,PD Notation,Genus-4D) ----

/// One row of a literature table, as written: its name, its PD code and its
/// 4-genus interval (lo == hi but for a few "[lo;hi]" rows).
struct LiteratureRow {
    std::string name;
    std::string pd;
    int lo = 0, hi = 0;
};

/// Splits one table row into its name, PD and genus field (no quoting
/// appears in these files). False for a row without two commas.
bool splitInputLine(const std::string &line, std::string &name,
                    std::string &pd, std::string &genusField);

/// Parses "N" or "[lo;hi]" into lo/hi (lo == hi in the plain-integer case).
void parseGenusField(const std::string &field, int &lo, int &hi);

/// Every row of a literature table, in file order: the header line, empty
/// lines, rows without two commas and rows whose genus field does not parse
/// are skipped. PD codes are kept as text, never parsed.
/// \throws std::runtime_error if the file cannot be opened.
std::vector<LiteratureRow> readLiteratureRows(const std::filesystem::path &path);

} // namespace witnessstore

#endif

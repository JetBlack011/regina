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
#include <optional>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <link/link.h>

namespace exactnaming {

// ---- the table files ----

/** One row of a knot or link table ("Name,PD,Genus-4D"), as written. */
struct TableRow {
    std::string name; /**< "3_1", "L7n1{1}" */
    std::string pd;   /**< the PD field, as written */
    std::string g4;   /**< the 4-genus field, as written: "1", "[0;1]" */
};

/**
 * Every row of a knot or link table, in file order: the name before the
 * first comma, the PD code between the first and the second, the 4-genus
 * after the second (no quoting appears in these files: a PD code uses ';'
 * internally, never a comma). The header line, empty lines and lines without
 * two commas are skipped; a trailing '\r' is dropped. The one reader of the
 * tables: the solver's literature, the search's input rows, the namer's
 * indices and the cascade's PD lookups all read through it.
 * \exception regina::InvalidArgument the file cannot be opened.
 */
std::vector<TableRow> readTableRows(const std::filesystem::path &path);

/**
 * A table's literature 4-genus, "2" or "[0;1]", as (lo, hi); nullopt if the
 * field is malformed (anything but digits, or a "[lo;hi]" with lo > hi), so
 * that a malformed value never becomes a bound.
 */
std::optional<std::pair<int, int>> parseTableG4(const std::string &field);

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

/** Knot name -> symmetry type. */
using SymmetryTable = std::unordered_map<std::string, SymmetryType>;

/**
 * data/knot_symmetry.csv ("name,symmetry_type,..."; header skipped): every
 * row whose type parses. The one reader of that file.
 * \exception regina::InvalidArgument the file cannot be opened.
 */
SymmetryTable readSymmetryTable(const std::filesystem::path &path);

/**
 * Whether the knot `name` is slice by ELEMENTARY concordance-group reasoning,
 * making it an anchor exactly like the unknot.
 *
 * `K # -K` bounds a ribbon disc for every K, where `-K = m(K^r)` is the
 * concordance inverse, and a sum of slice knots is slice. So a composite knot
 * is elementarily slice when its summands pair off into inverse pairs. The
 * inverse depends on the summand's symmetry: for an invertible K, `-K = mK`;
 * if K is also amphicheiral, `-K = K`. Summands without a certified symmetry
 * (absent from `symmetry`), and every NON-invertible summand, are refused:
 * our names record chirality but not reversal, so for a non-invertible K the
 * name cannot say whether its neighbour is -K or its reverse (8_17 is the
 * trap).
 *
 * The two long-standing anchors "3_1#m3_1" and "4_1#4_1" are accepted even
 * with no symmetry data loaded, so a run without --knot-symmetry loses
 * nothing it had before. Consumed by both the atlas solver and the cascade.
 */
bool isElementarySlice(const std::string &name, const SymmetryTable &symmetry);

// ---- the tables, indexed for naming ----

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
    /** A knot's symmetry type, or nullopt when knot_symmetry.csv has none. */
    std::optional<SymmetryType> symmetry(const std::string &knot) const;
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
    SymmetryTable symmetry_;
};

/**
 * A PD code as text ("[[1;5;2;4];...]", "PD[X[4; 1; 3; 2]; ...]", or any
 * other punctuation: the integers in fours) as Regina reads it, labels as
 * written. The one parser of PD text into a regina::Link: the tables', and a
 * cascade row's. A code holding a 0 is taken as 0-based and shifted up by
 * one, as diagramtriangulation's parsePDCode() reads it.
 * \exception regina::InvalidArgument no labels, or a count not divisible by 4.
 */
regina::Link linkFromTablePD(const std::string &pd);

} // namespace exactnaming

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

} // namespace farside

#endif

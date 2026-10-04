//
//  tables.h
//
//  The knot and link tables as oriented diagrams, for naming outgoing links
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

#ifndef SURFER_LINKNAMING_TABLES_H
#define SURFER_LINKNAMING_TABLES_H

#include <filesystem>
#include <optional>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <link/link.h>

namespace linknaming {

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
 * tables: the solver's literature, a run's targets, the namer's indices
 * and a goal run's PD lookups all read through it.
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
 * is elementarily slice when its summands pair off into inverse pairs, each
 * summand read with its marks (`m` mirrored, `r` reversed, relative to the
 * table's diagram) and its symmetry type, which says which marks change the
 * knot: for a reversible K, `-K = mK` (`3_1#m3_1`, the `r` meaningless); for a
 * fully amphicheiral one, `-K = K`; for a negative amphicheiral one `-K = K`
 * (`8_17#8_17`, while `8_17#m8_17` = `8_17 # 8_17^r` is not slice: the
 * trap); for a positive amphicheiral one `-K = K^r` (`12a_1#r12a_1`); and for
 * a chiral one only the mark `mr` gives the inverse (`9_32#mr9_32`, never
 * `9_32#m9_32`). A summand without a symmetry type (absent from `symmetry`)
 * is refused, never guessed.
 *
 * The two long-standing anchors "3_1#m3_1" and "4_1#4_1" are accepted even
 * with no symmetry data loaded, so a run without `knot_symmetry` still
 * has them. Consumed by both the atlas solver and goal runs.
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

class Tables {
  public:
    Tables() = default;
    Tables(Tables &&) noexcept = default;
    Tables &operator=(Tables &&) noexcept = default;
    Tables(const Tables &) = delete; // the indices point into entries_
    Tables &operator=(const Tables &) = delete;

    /**
     * \param knotTable, linkTable "Name,PD,Genus-4D" CSVs (the atlas tables).
     * \param knotSymmetry "name,symmetry_type,..." (may be empty: every knot
     *        then counts as chiral, which only ever withholds exactness).
     * \exception regina::InvalidArgument a table cannot be read, or two
     *        different bases share an unoriented diagram.
     */
    static Tables load(const std::string &knotTable, const std::string &linkTable,
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
    /** Every entry, in load order: the knot table's rows, then the link table's. */
    const std::vector<TableEntry> &entries() const { return entries_; }
    /** How many of entries() came from the knot table (they come first). */
    size_t knotEntries() const { return knotEntries_; }

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
    size_t knotEntries_ = 0;
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
 * searched link's. A code holding a 0 is taken as 0-based and shifted up by
 * one, as diagramtriangulation's parsePDCode() reads it.
 * \exception regina::InvalidArgument no labels, or a count not divisible by 4.
 */
regina::Link linkFromTablePD(const std::string &pd);

} // namespace linknaming

#endif

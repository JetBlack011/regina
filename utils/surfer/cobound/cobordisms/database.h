//
//  database.h
//
//  The append-only cobordism database (cobordisms.csv), shared by
//  every run and command. Moved from verifyslicegenus.cpp on
//  2026-09-28 (as witnessstore.h) so that every program that finds
//  cobordisms writes them in one format, through one append path.
//

#pragma once

#include <filesystem>
#include <mutex>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "cobound/cobordisms/cobordism.h"

namespace cobordisms {

/// The database's header (13 columns, resolved_vertices last).
inline constexpr const char *COBORDISMS_HEADER =
    "kind,subject,subject_components,other,other_candidates,other_components,"
    "genus,tubed,pairsig,source_row,thicken_layers,max_faces,"
    "resolved_vertices";

// ---- the database ----

/// One cobordism as a line of the database (no newline).
std::string formatCobordism(const cobordisms::Cobordism &w);

/// A cobordism from the first 12 or 13 fields of a cobordism line (split by
/// parseCsvLine()); fields past the 13th are ignored. The pair signature is
/// kept only if `keepPairSig`; `path` is for messages. False if malformed.
bool cobordismFromFields(std::vector<std::string> fields,
                       cobordisms::Cobordism &w, bool keepPairSig,
                       bool wantPairSigKey, const std::filesystem::path &path);

/// The other_candidates field split into names: at every ';' outside a
/// `{...}` orientation tag (a tag lists its entries with ';'), empty pieces
/// dropped.
std::vector<std::string> splitCandidates(const std::string &field);

/// One cobordism line (12 or 13 fields, parseCsvLine()) into `w`, as
/// cobordismFromFields(); false if malformed.
bool parseCobordismLine(const std::string &line, cobordisms::Cobordism &w,
                      bool keepPairSig, bool wantPairSigKey,
                      const std::filesystem::path &path);

/// Every complete cobordism line of `path`, WITHOUT pair signatures: each
/// keeps its line's byte offset (Cobordism::fileOffset). A torn last line is
/// ignored. A missing file is an empty database.
std::vector<cobordisms::Cobordism>
loadCobordisms(const std::filesystem::path &path, bool wantPairSigKeys);

/// As loadCobordisms(), but every cobordism keeps its pair signature: for a
/// small file read whole (a child search's cob.csv).
std::vector<cobordisms::Cobordism> readCobordisms(const std::filesystem::path &path);

/// A cobordism's pair signature read back from its line (Cobordism::fileOffset),
/// memoized: the few places that print one ask for the same few repeatedly.
class PairSigReader {
  public:
    void setPath(std::filesystem::path path);
    /// "" for a negative offset. Throws if there is no cobordism line there.
    std::string at(long long offset);

  private:
    std::mutex mutex_;
    std::filesystem::path path_;
    std::unordered_map<long long, std::string> cache_;
};

/// A stored cobordism as an index reads it back: its cobordism (pair signature
/// kept) and, when the database's `.rows.csv` sidecar records one, the
/// diagram its search ran on (empty: the subject's table PD).
struct StoredCobordism {
    cobordisms::Cobordism cobordism;
    std::string incomingPD;
};

/**
 * A database read only, indexed: one scan records where each subject's lines
 * start and which lines name each outgoing base; lines are read on demand.
 * A torn last line is left out, as loadCobordisms() leaves it out.
 */
class DatabaseIndex {
  public:
    /// A table name's base: orientation tag and a knot's mirror prefix
    /// dropped. Used only to FIND candidate cobordisms; every one found is
    /// redrawn and identified exactly before it means anything.
    static std::string base(std::string name);

    /// \exception std::runtime_error `path` cannot be read.
    explicit DatabaseIndex(const std::string &path);

    bool has(const std::string &subject) const { return offsets_.count(subject) > 0; }
    /// The cobordisms whose subject is `subject`.
    std::vector<StoredCobordism> ofSubject(const std::string &subject) const;
    /// The cobordisms of OTHER subjects whose recorded outgoing link has this
    /// base name (a hint only); at most `cap`.
    std::vector<StoredCobordism> byOutgoing(const std::string &base, size_t cap) const;
    size_t subjects() const { return offsets_.size(); }

  private:
    std::vector<StoredCobordism> read(const std::vector<std::streamoff> &offsets,
                                      size_t cap) const;
    std::string path_;
    std::unordered_map<std::string, std::vector<std::streamoff>> offsets_;
    std::unordered_map<std::string, std::vector<std::streamoff>> byOther_;
    std::unordered_map<std::string, std::string> incomingPD_; ///< cobordism key -> incoming PD
};

/// cobordisms::cobordismIdentity() of every cobordism loadCobordisms() would
/// load from `path` -- the same lines, torn last line and malformed lines
/// skipped alike -- parsed in `threads` byte ranges cut at line starts. A
/// database's dedupe needs only this, and the atlas's master is ~3 GB, so a
/// run's signing step spent 8-12 s of one thread here (2026-09-30).
std::unordered_set<std::string>
cobordismIdentities(const std::filesystem::path &path, unsigned threads,
                  std::streamoff minRangeBytes = 64 << 20);

/// Appends cobordisms[from..] to `path` and fsyncs; creates it with the
/// header if absent, truncates a torn last line first, and refuses a file
/// with any other header. Each appended cobordism gets its fileOffset and
/// drops its pair signature. Throws on failure, leaving them untouched.
void appendCobordisms(const std::filesystem::path &path,
                     std::vector<cobordisms::Cobordism> &cobordisms,
                     size_t from);

} // namespace cobordisms

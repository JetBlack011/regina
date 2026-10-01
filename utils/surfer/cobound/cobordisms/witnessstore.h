//
//  witnessstore.h
//
//  The append-only witness file (cobordisms.csv) and the literature tables'
//  CSV rows, shared by verifyslicegenus and cascadesearch. Moved from
//  verifyslicegenus.cpp on 2026-09-28 so that every program that finds
//  witnesses writes them in one format, through one append path.
//

#pragma once

#include <filesystem>
#include <string>
#include <unordered_set>
#include <vector>

#include "cobound/cobordisms/cobordism.h"
#include "cobound/solver/literature.h"

namespace witnessstore {

/// The witness file's header (13 columns, resolved_vertices last).
inline constexpr const char *COBORDISMS_HEADER =
    "kind,subject,subject_components,other,other_candidates,other_components,"
    "genus,tubed,pairsig,source_row,thicken_layers,max_faces,"
    "resolved_vertices";

// ---- literature tables (Name,PD Notation,Genus-4D) ----

/// Splits one table row into its name, PD and genus field (no quoting
/// appears in these files). False for a row without two commas.
bool splitInputLine(const std::string &line, std::string &name,
                    std::string &pd, std::string &genusField);

/// Parses "N" or "[lo;hi]" into lo/hi (lo == hi in the plain-integer case).
void parseGenusField(const std::string &field, int &lo, int &hi);

/// Registers every row of a literature table in `names` (names and bounds
/// only; PD codes are skipped). Returns the number loaded.
size_t loadNameTable(const std::filesystem::path &path,
                     cobordismgraph::NameTable &names);

// ---- the witness file ----

/// One witness as a line of the witness file (no newline).
std::string formatWitness(const cobordismgraph::Witness &w);

/// A witness from the first 12 or 13 fields of a witness line (split by
/// parseCsvLine()); fields past the 13th are ignored. The pair signature is
/// kept only if `keepPairSig`; `path` is for messages. False if malformed.
bool witnessFromFields(std::vector<std::string> fields,
                       cobordismgraph::Witness &w, bool keepPairSig,
                       bool wantPairSigKey, const std::filesystem::path &path);

/// Every complete witness line of `path`, WITHOUT pair signatures: each
/// keeps its line's byte offset (Witness::fileOffset). A torn last line is
/// ignored. A missing file is an empty store.
std::vector<cobordismgraph::Witness>
loadWitnesses(const std::filesystem::path &path, bool wantPairSigKeys);

/// cobordismgraph::witnessIdentity() of every witness loadWitnesses() would
/// load from `path` -- the same lines, torn last line and malformed lines
/// skipped alike -- parsed in `threads` byte ranges cut at line starts. A
/// store's dedupe needs only this, and the atlas's master is ~3 GB, so a
/// run's store step spent 8-12 s of one thread here (2026-09-30).
std::unordered_set<std::string>
witnessIdentities(const std::filesystem::path &path, unsigned threads,
                  std::streamoff minRangeBytes = 64 << 20);

/// Appends witnesses[from..] to `path` and fsyncs; creates it with the
/// header if absent, truncates a torn last line first, and refuses a file
/// with any other header. Each appended witness gets its fileOffset and
/// drops its pair signature. Throws on failure, leaving them untouched.
void appendWitnesses(const std::filesystem::path &path,
                     std::vector<cobordismgraph::Witness> &witnesses,
                     size_t from);

/// --rewrite-witnesses: the one full rewrite (the 12->13 column
/// migration), verifying that every line round-trips. Returns the count.
size_t rewriteWitnessFile(const std::filesystem::path &path);

} // namespace witnessstore

// pending.h
//
// Every surface a search keeps, into the atlas's database, as a run without
// a goal would record it (README.md, "The database and signing").
//
// An in-process search keeps each surface's faces, not its pair signature: a
// signature's cost is almost all the thickening's own (~50 s for a 10-crossing
// link), and a search takes seconds. So each search appends its kept surfaces to
// <search dir>/kept.csv at once (cheap, durable), and the pair signatures are
// computed after the search: one ambient per incoming diagram, diagrams in
// parallel (pairSigsOf()), only for surfaces whose cobordism is new to the database.

#pragma once

#include <atomic>
#include <filesystem>
#include <functional>
#include <mutex>
#include <string>
#include <unordered_set>
#include <vector>

#include "cobound/cobordisms/cobordism.h"
#include "cobound/solver/literature.h"

namespace cobordisms {

/// A kept surface bound for the database: its cobordism (every column
/// but the pair signature), and what its signature is computed from.
struct PendingCobordism {
  cobordisms::Cobordism cobordism;
  std::string incomingPD;
  int layers = 2;
  std::vector<int> faces;
};

/// One kept.csv line (with its newline): the cobordism's database line (pair
/// signature empty), then its faces (space-separated), its incoming PD and layers.
std::string formatKept(const PendingCobordism &p);

/// Appends `kept` to <searchDir>/kept.csv and fsyncs (formatKept()).
void appendKept(const std::string &searchDir, const std::vector<PendingCobordism> &kept);

/**
 * A search's pending file (kept.csv), written while the search runs: every
 * surface it keeps is added (add(), from drain threads), and what was added
 * reaches the file by appending and fsyncing -- once a minute
 * (checkpoint(false), from the search's and the drain's progress), at once
 * on demand (checkpoint(true): the first constructive find), and at the end
 * (flush()).
 * A kill therefore loses at most a minute of found surfaces, and `sign`
 * (cobound sign) signs what reached the file. Plan divergence 7.
 */
class PendingWriter {
public:
  explicit PendingWriter(std::filesystem::path path);
  PendingWriter(const PendingWriter &) = delete;
  PendingWriter &operator=(const PendingWriter &) = delete;

  /// Queues one kept surface's line. Thread-safe.
  void add(const PendingCobordism &p);
  /// Appends and fsyncs what was added since the last write, unless the
  /// last write was under 60 s ago; `force` writes at once. Returns why a
  /// write failed, or "" -- never throws, since it runs on the search's
  /// worker, drain and judge threads: the caller ends its search as an I/O
  /// error (search.h, SearchResult::ioFailure). What is queued stays queued,
  /// and the next write retries it. Thread-safe.
  std::string checkpoint(bool force);
  /// The end of the search: appends and fsyncs the rest. Throws on failure.
  void flush();
  /// The file's length after the last successful write (its length when
  /// opened, before any): what a frontier taken now may rely on.
  long long syncedBytes() const;
  const std::filesystem::path &path() const { return path_; }

private:
  void write_(); // under mutex_
  std::filesystem::path path_;
  mutable std::mutex mutex_;
  std::string queued_;
  long long synced_ = 0;
  std::atomic<long long> lastWriteTick_{0};
};


struct SignResult {
  size_t kept = 0;     ///< surfaces offered
  size_t fresh = 0;    ///< of those, cobordism identities new to the database
  size_t appended = 0; ///< written
  double dedupeSeconds = 0; ///< reading the databases' identities (dedupe_against too)
  double signSeconds = 0;
};

/// What a database already holds, as a run loaded it: the identities of its
/// first `bytes` bytes. signKept() then reads only what was appended since.
struct LoadedPrefix {
  const std::unordered_set<std::string> *identities = nullptr;
  std::uintmax_t bytes = 0;
};

/// Every <work>/hop_*/kept.csv, in search order, from where `sign` last signed
/// it (signedThrough()) to its last complete line (a torn last line is
/// skipped). `readTo`, if given, gets each file read and the byte it was
/// read to.
std::vector<PendingCobordism>
readKept(const std::string &work,
         std::vector<std::pair<std::string, long long>> *readTo = nullptr);

/// How far `sign` has signed the pending file `path` (its <path>.signed
/// record): 0 if never.
long long signedThrough(const std::string &path);

/// `sign`: signKept() over readKept(work), then each pending file's
/// <path>.signed record set to what was read (plan divergence 1: a frontier
/// whose pending file is signed at least as far as it recorded may be
/// resumed). Returns signKept()'s result.
SignResult signPending(const std::string &work, const std::string &database,
                        const std::vector<std::string> &dedupeAgainst,
                        const solver::NameTable &names, unsigned threads,
                        const std::string &pairSigCache = "",
                        const LoadedPrefix &loaded = {},
                        const std::function<bool(const PendingCobordism &)> &sidecarLine = {});


/**
 * Signs kept surfaces into the database. A surface is fresh when its
 * cobordisms::cobordismIdentity() is in neither `database` nor any of
 * `dedupeAgainst` (read only; the master, say), nor earlier in `pending`:
 * exactly a run without a goal's rule, which records one cobordism per identity. Only
 * fresh surfaces are signed (pairSigsOf(), `threads` diagrams at a time).
 * Their other_candidates come from `names`, as a run's without a goal do. The append
 * holds an exclusive lock on <database>.lock and re-reads the database's
 * identities under it, so runs sharing a database never record one identity
 * twice. With `loaded`, the identities the run already holds stand for the
 * database's first `loaded.bytes` bytes, and only the rest is read.
 * `sidecarLine`, when given, says which appended cobordisms get a
 * `.rows.csv` line (by default every one).
 */
SignResult signKept(std::vector<PendingCobordism> pending, const std::string &database,
                      const std::vector<std::string> &dedupeAgainst,
                      const solver::NameTable &names, unsigned threads,
                      const std::string &pairSigCache = "",
                      const LoadedPrefix &loaded = {},
                      const std::function<bool(const PendingCobordism &)> &sidecarLine = {});

/**
 * The cobordism database as a run loaded it: its cobordisms (the solver's
 * input), their identities (what a search dedupes its finds against: plan
 * divergence 7, one cobordism per identity across the database), and the
 * file's length when read (what signKept() need not read again). A search
 * never writes it: its finds go to its pending file, and `sign` appends them.
 */
class LoadedDatabase {
public:
  /// `loaded`: the cobordisms of the database file at `path` (loadCobordisms()).
  LoadedDatabase(std::filesystem::path path, std::vector<cobordisms::Cobordism> loaded);
  LoadedDatabase(const LoadedDatabase &) = delete;
  LoadedDatabase &operator=(const LoadedDatabase &) = delete;

  /// Every cobordism, as loaded.
  const std::vector<cobordisms::Cobordism> &all() const { return cobordisms_; }
  /// cobordisms::cobordismIdentity() of each.
  const std::unordered_set<std::string> &identities() const { return identities_; }
  /// The database's first bytes() bytes are what was loaded.
  LoadedPrefix loaded() const { return {&identities_, bytes_}; }

private:
  std::filesystem::path path_;
  std::vector<cobordisms::Cobordism> cobordisms_;
  std::unordered_set<std::string> identities_;
  std::uintmax_t bytes_ = 0;
};

} // namespace cobordisms

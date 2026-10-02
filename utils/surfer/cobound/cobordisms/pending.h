// keptstore.h
//
// Every surface a hop keeps, into the atlas's witness store, as the sweep
// would have recorded it (README.md, "Witnesses for the atlas").
//
// An in-process hop keeps each surface's faces, not its pair signature: a
// signature's cost is almost all the row's ambient (~50 s for a 10-crossing
// row), and a hop takes seconds. So each hop appends its kept surfaces to
// <hop dir>/kept.csv at once (cheap, durable), and the pair signatures are
// computed after the search: one ambient per row, rows in parallel
// (pairSigsOf()), only for surfaces whose witness is new to the store.

#pragma once

#include <atomic>
#include <filesystem>
#include <mutex>
#include <string>
#include <unordered_set>
#include <vector>

#include "cobound/cobordisms/cobordism.h"
#include "cobound/solver/literature.h"

namespace cascade {

/// A kept surface bound for the witness store: its witness (every column
/// but the pair signature), and what its signature is computed from.
struct PendingWitness {
  cobordismgraph::Witness witness;
  std::string rowPD;
  int layers = 2;
  std::vector<int> faces;
};

/// Appends `kept` to <hopDir>/kept.csv and fsyncs. Each line is the
/// witness's store line (pair signature empty), then its faces
/// (space-separated), its row PD and layers.
void appendKept(const std::string &hopDir, const std::vector<PendingWitness> &kept);

/// Every <work>/hop_*/kept.csv, in hop order. A torn last line is skipped.
std::vector<PendingWitness> readKept(const std::string &work);

struct StoreResult {
  size_t kept = 0;     ///< surfaces offered
  size_t fresh = 0;    ///< of those, witness identities new to the store
  size_t appended = 0; ///< written
  double dedupeSeconds = 0; ///< reading the stores' identities (--dedupe-against too)
  double signSeconds = 0;
};

/**
 * Signs and stores kept surfaces. A surface is fresh when its
 * cobordismgraph::witnessIdentity() is in neither `store` nor any of
 * `dedupeAgainst` (read only; the master, say), nor earlier in `pending`:
 * exactly the sweep's rule, which records one witness per identity. Only
 * fresh surfaces are signed (pairSigsOf(), `threads` rows at a time).
 * Their other_candidates come from `names`, as the sweep's do. The append
 * holds an exclusive lock on <store>.lock and re-reads the store's
 * identities under it, so runs sharing a store never record one identity
 * twice.
 */
StoreResult storeKept(std::vector<PendingWitness> pending, const std::string &store,
                      const std::vector<std::string> &dedupeAgainst,
                      const cobordismgraph::NameTable &names, unsigned threads,
                      const std::string &pairSigCache = "");

/**
 * The witnesses of a run that signs during its searches (verifyslicegenus,
 * SearchPolicy::Signing::duringSearch): the database's own, as loaded, and
 * every one its searches have found since, one per witness identity across
 * them all.
 *
 * A search claims a new witness's identity when it finds it (claim()), so
 * the dedupe is exact, and publishes it once it is signed (publish(), as
 * WitnessSigner's output). What is published reaches the database file only
 * by appending -- at checkpoints (checkpoint()), and at the end of each
 * search (flush()) -- so nothing on disk is ever rewritten.
 *
 * Checkpointing exists because a search's witnesses would otherwise be
 * written only when it ends, so a search that is interrupted -- by a crash,
 * a shutdown, or an operator stopping a run that looks unproductive --
 * loses everything it found. That is not hypothetical: a 4-hour L9n2{1} row
 * lost 13 hours to a shutdown mid-drain, and an L9a26{1} row was killed five
 * hours after it had already found a constructive genus-0 witness that had
 * never reached disk. Under --harvest the exposure is worst, because a row
 * that has ALREADY resolved deliberately keeps running to bank more edges.
 *
 * The plan's divergence 7 (signing) replaces this with pending surfaces,
 * appended unsigned and signed at the end of the run.
 */
class RecordedWitnesses {
public:
  /// `loaded`: the witnesses of the database file at `path` (loadWitnesses()),
  /// which every later write appends to.
  RecordedWitnesses(std::filesystem::path path, std::vector<cobordismgraph::Witness> loaded);
  RecordedWitnesses(const RecordedWitnesses &) = delete;
  RecordedWitnesses &operator=(const RecordedWitnesses &) = delete;

  /// Claims `w`'s witness identity: false when a witness with it is recorded
  /// (or claimed) already. Thread-safe.
  bool claim(const cobordismgraph::Witness &w);

  /// A claimed witness, signed, into the record. Thread-safe.
  void publish(cobordismgraph::Witness &&w);

  /// Appends whatever was published since the last write, unless the last
  /// checkpoint was under 60 s ago; `force` writes at once regardless. A
  /// failure is reported on stderr, never thrown: a failed checkpoint must
  /// not kill a running search, and the next write retries. Thread-safe.
  void checkpoint(bool force);

  /// The end-of-search and end-of-run write: appends whatever is new and --
  /// unlike a checkpoint -- lets a failure propagate. Writes nothing at all
  /// when nothing is new, so a run that searched nothing never touches the
  /// file.
  void flush();

  /// Restarts the quiescence clock, as a search begins.
  void markActivity();
  /// Milliseconds since the last new witness was claimed (or markActivity()).
  long long millisSinceNew() const;

  /// Every witness, loaded and found. Read only while no search runs.
  const std::vector<cobordismgraph::Witness> &all() const { return witnesses_; }

private:
  std::filesystem::path path_;
  std::mutex mutex_;
  std::vector<cobordismgraph::Witness> witnesses_;
  // witnessIdentity() of every witness in witnesses_, and of every witness
  // claimed and not yet published, under mutex_: the dedupe is a hash
  // lookup, not a scan of every witness ever recorded while holding the one
  // lock.
  std::unordered_set<std::string> identities_;
  // How far the file is current. witnesses_ is append-only -- publish() is
  // its only mutation after construction -- so its size changes if and only
  // if its content does, which makes the count an exact dirty flag rather
  // than a heuristic one. Read and written only under mutex_. Initialised
  // from the loaded set, so a run that resumes does not rewrite what it
  // just read.
  size_t lastCheckpointedCount_;
  std::atomic<long long> lastCheckpointTick_{0};
  std::atomic<long long> lastNewTick_{0};
};

} // namespace cascade

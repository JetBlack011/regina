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

#include <string>
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

} // namespace cascade

// hoprunner.h
//
// One hop, searched in process: the row's own thickening (the one its
// HopAssembler certified and reads from), the campaign's search shape, and
// the gates and accounting verifyslicegenus applies (rowsearch.h), with the
// accepted surfaces read straight off the search instead of written out as
// pair signatures and read back.
//
// What that saves, per hop: a process start and table load, a pair
// signature per kept surface (3.4 s each on a 10-crossing row, 2026-09-28:
// 64 s of a 152 s hop), and decoding them again. A kept surface keeps its
// faces instead, from which its pair signature can be computed later
// (pairSig() of the rebuilt thickening), for the few witnesses a
// certificate needs.

#pragma once

#include <functional>
#include <string>
#include <vector>

#include "exactnaming/exacttables.h"
#include "farsidecurves.h"
#include "farsidenaming.h"
#include "farsideredraw.h"

namespace cascade {

/// The campaign's search shape (atlas tools/orchestrate/hosts.conf
/// [campaign], root budget 840 from c5 on). The layers are the row's.
struct HopShape {
  long long maxFaces = 5;
  unsigned iddfsIterations = 2;
  long long iddfsStart = 4;
  long long iddfsStep = 1;
  long long rootBudgetStart = 840;
  long long rootBudgetGrowth = 2;
  bool resolveUnlinked = true;
  /// hosts.conf's per-host limits, identical on every host. The pending
  /// cap in particular: at the binary's default (500,000) a search pauses
  /// to drain its queue, which a campaign row never does.
  size_t pendingSurfaceCap = 20'000'000;
  size_t petalCacheLimit = 12'000'000;
  size_t boundarySignatureCacheLimit = 1'000'000;
  /// Process-wide (identify::recognitionCacheLimit); set by the driver.
  size_t recognitionCacheLimit = 1'500'000;
};

/// A surface a hop kept.
struct KeptSurface {
  farside::OutgoingLink link; ///< far side on T, oriented, and incoming side
  int genus = 0;              ///< tubed genus
  int resolvedVertices = 0;
  std::string farName;        ///< the namer's name; "" for no far side
  std::vector<int> faces;     ///< triangles of the row's thickening
  std::string key;            ///< its dedupe key (see HopSearcher::run())
};

/// What a hop did.
struct HopRun {
  std::vector<KeptSurface> kept;
  long long accepted = 0;
  std::string accounting;        ///< the `accounting:` body (rowsearch.h)
  std::string accountingFailure; ///< empty iff every surface is accounted for
  std::string outcome;           ///< surface-target, exhausted, timeout, stopped
  double wall = 0;               ///< seconds
  double cpu = 0;                ///< process CPU seconds over the hop
  double setup = 0;              ///< wall before the search: namer, search, seed checks
  double search = 0;             ///< wall of the search itself, drain included
  std::vector<double> rounds;    ///< each IDDFS round's wall (SearchStats::Profile)
  size_t drainTail = 0;          ///< surfaces left for the drain when the search ended
  double drainTailSeconds = 0;   ///< and the wall it took to describe them
};

class HopSearcher {
public:
  /// `signatures` and `exact` name far sides as verifyslicegenus names them
  /// (farside::DiagramNamer), which is what surfaces are deduplicated by.
  /// Both must outlive the searcher.
  HopSearcher(const farside::SignatureTable &signatures,
              const exactnaming::ExactTables *exact, HopShape shape,
              unsigned threads);

  /**
   * Searches `row`'s thickening, seeded with its collar, under `proper`,
   * until `surfaceTarget` surfaces satisfy it (or `seconds` pass; the drain
   * then finishes), and keeps one surface per key: verifyslicegenus's
   * witness identity (far-side name, its component count, genus, whether
   * tubed or resolved) together with which row components and how many
   * far-side curves each surface component carries. The second part is
   * what profiles read (partitions), and the witness identity alone
   * collapses it.
   *
   * `stop`, when given, is polled with each newly kept surface (from drain
   * threads, one at a time) and ends the search when it returns true.
   *
   * \throws std::runtime_error when the row's seed invariant fails.
   */
  HopRun run(const farside::WitnessRedrawer &row, const std::string &rowName,
             long long surfaceTarget, double seconds,
             const std::function<bool(const KeptSurface &)> &stop = {}) const;

private:
  const farside::SignatureTable &signatures_;
  const exactnaming::ExactTables *exact_;
  HopShape shape_;
  unsigned threads_;
};

/// The pair signature of a kept surface, from its faces and the thickening
/// it was found in (rebuilt from the row's PD by rowsearch::buildRow(), or
/// the searched one itself).
std::string pairSigOf(const regina::Triangulation<4> &thickening,
                      const std::vector<int> &faces);

/// A kept surface to sign: the row it was found in, and its faces there.
struct SignRequest {
  std::string rowPD;
  int layers = 2;
  std::vector<int> faces;
};

/// pairSigOf() for many surfaces at once. Almost all of a signature's cost is
/// the ambient's own (~50 s for a 10-crossing row, 2026-09-28), so each
/// distinct row is rebuilt (rowsearch::buildRow(), deterministic, so faces
/// index it as they did the searched one) and its ambient part computed once
/// (PairSigContext), up to `threads` rows at a time. In request order.
std::vector<std::string> pairSigsOf(const std::vector<SignRequest> &requests,
                                    unsigned threads);

} // namespace cascade

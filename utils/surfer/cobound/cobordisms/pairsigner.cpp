//
//  pairsigner.cpp
//

#include "cobound/cobordisms/pairsigner.h"

#include <algorithm>
#include <atomic>
#include <map>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <thread>

#include "cobound/parallelfor.h"
#include "surfer/pairsig/pairsig.h"
#include "cobound/search/incoming.h"

namespace cascade {

std::string pairSigOf(const regina::Triangulation<4> &thickening,
                      const std::vector<int> &faces) {
  return pairSig<4, 2>(thickening, faces);
}

std::unique_ptr<PairSigContext<4, 2>>
pairSigContextFor(const regina::Triangulation<4> &thickening, const std::string &cacheDir,
                  unsigned threads, bool *loaded) {
  if (cacheDir.empty()) {
    if (loaded) *loaded = false;
    return std::make_unique<PairSigContext<4, 2>>(thickening, threads);
  }
  return PairSigContext<4, 2>::cached(thickening, cacheDir, loaded, threads);
}

std::vector<std::string> pairSigsOf(const std::vector<SignRequest> &requests,
                                    unsigned threads, const std::string &cacheDir) {
  std::map<std::pair<std::string, int>, std::vector<size_t>> byRow;
  for (size_t i = 0; i < requests.size(); ++i)
    byRow[{requests[i].rowPD, requests[i].layers}].push_back(i);
  std::vector<const std::pair<const std::pair<std::string, int>, std::vector<size_t>> *> rows;
  for (const auto &entry : byRow) rows.push_back(&entry);

  std::vector<std::string> out(requests.size());
  std::mutex errorMutex;
  std::string error;
  // Rows share the threads: one row per worker, and each row's context (its
  // ambient's isoSig, nearly all of signing) on threads / workers of them,
  // so a run that kept surfaces from one row -- every 12-crossing knot row
  // of a campaign -- signs on all of them rather than one.
  if (rows.empty()) return out;
  const size_t n = std::min<size_t>(std::max(threads, 1u), rows.size());
  const unsigned inner = std::max(1u, static_cast<unsigned>(std::max(threads, 1u) / n));
  parallelFor(rows.size(), static_cast<unsigned>(n), [&](size_t r) {
    try {
      const auto &[row, indices] = *rows[r];
      rowsearch::RowBuild rb;
      rowsearch::buildRow(row.first, row.second, row.second, /*useCone=*/false, rb);
      const std::unique_ptr<PairSigContext<4, 2>> context =
          pairSigContextFor(rb.tri, cacheDir, inner);
      for (size_t i : indices) out[i] = context->sig(requests[i].faces);
    } catch (const std::exception &e) {
      std::lock_guard<std::mutex> lock(errorMutex);
      error = e.what();
    }
  });
  if (!error.empty()) throw std::runtime_error("pairSigsOf: " + error);
  return out;
}

} // namespace cascade

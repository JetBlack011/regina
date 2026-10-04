//
//  pairsigner.h
//
//  Pair signatures for cobordisms found by a search.
//

/*! \file utils/surfer/cobound/cobordisms/pairsigner.h
 *  \brief Pair signatures for the surfaces a search kept, from their faces
 *  and the thickening they were found in: one at a time, or many at once
 *  with each thickening's ambient part (PairSigContext) computed once.
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include <triangulation/dim4.h>

#include "cobound/cobordisms/cobordism.h"
#include "surfer/pairsig/pairsig.h"

namespace cobordisms {

/// The pair signature of a kept surface, from its faces and the thickening
/// it was found in (rebuilt from the incoming PD by search::buildIncoming(), or
/// the searched one itself).
std::string pairSigOf(const regina::Triangulation<4> &thickening,
                      const std::vector<int> &faces);

/// A kept surface to sign: the incoming diagram it was found on, and its faces there.
struct SignRequest {
  std::string incomingPD;
  int layers = 2;
  std::vector<int> faces;
};

/// pairSigOf() for many surfaces at once. Almost all of a signature's cost is
/// the ambient's own (~50 s for a 10-crossing link, 2026-09-28), so each
/// distinct incoming diagram is rebuilt (search::buildIncoming(), deterministic, so faces
/// index it as they did the searched one) and its ambient part computed once
/// (PairSigContext), up to `threads` diagrams at a time. In request order.
/// With a `cacheDir`, each diagram's context is read from there when stored and
/// stored when built (PairSigContext::cached()), so a diagram signed again in a
/// later run -- a link goal runs meet often -- costs no rebuild.
std::vector<std::string> pairSigsOf(const std::vector<SignRequest> &requests,
                                    unsigned threads, const std::string &cacheDir = "");

/// A thickening's pair-signature context (its ambient part, nearly all of a
/// signature's cost) on `threads` threads: built in memory, or with a
/// `cacheDir`, through the one cache format (PairSigContext::cached(): a
/// `.pairsigctx` file per ambient, checked on every use). `loaded` (if
/// given) says whether the cache held it.
std::unique_ptr<PairSigContext<4, 2>>
pairSigContextFor(const regina::Triangulation<4> &thickening, const std::string &cacheDir,
                  unsigned threads = 1, bool *loaded = nullptr);

} // namespace cobordisms

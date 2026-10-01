//
//  sha1.h
//
//  SHA-1, as the digests the search code keys on need it.
//

#ifndef SURFER_PAIRSIG_SHA1_H
#define SURFER_PAIRSIG_SHA1_H

#include <string>

/*! \file utils/surfer/surfer/pairsig/sha1.h
 *  \brief The SHA-1 digest the search code keys on: a thin wrapper over
 *  OpenSSL's SHA1().
 *
 *  Three digests rest on it: a search frontier's fingerprint
 *  (EmbeddingSearch::frontierFingerprint_()), a pair-signature context's
 *  ambient key (PairSigContext::ambientKey()), and cobound's cobordism key,
 *  sha1(pairsig)[:12]. Each is stored -- in frontier files, in context cache
 *  file names, in the atlas's tables -- and the atlas computes the last with
 *  Python's hashlib.sha1, so this must agree with FIPS 180-1 byte for byte.
 *  It replaced a hand-written implementation only after agreeing with it on
 *  every input the gates produce (tests/sha1_test.cpp keeps the known-answer
 *  vectors).
 */

namespace pairsig {

/** The full 40-character lowercase hex SHA-1 digest of \a data. */
std::string sha1Hex(const std::string &data);

} // namespace pairsig

#endif // SURFER_PAIRSIG_SHA1_H

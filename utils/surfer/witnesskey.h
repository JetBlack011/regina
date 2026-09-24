//
//  witnesskey.h
//
//  Created by John Teague on 09/08/2026.
//

#ifndef WITNESSKEY_H

#define WITNESSKEY_H

#include <string>

/*! \file utils/surfer/witnesskey.h
 *  \brief The per-witness key a far-side resolution is stored under.
 *
 *  \section wk_why Why a witness key rather than a name
 *
 *  A far side reaches cobordisms.csv as a NAME produced by identify() from
 *  the complement alone. For a knot that is enough -- Gordon-Luecke makes
 *  the complement determine the knot up to mirroring, and g_4 is
 *  mirror-invariant. For a LINK it is not: a complement does not determine
 *  a link (Rolfsen twisting), and one observed name is genuinely different
 *  links on different witnesses -- "m129" alone stands for a dozen of them
 *  in our data. So a name-keyed table (--name-aliases) cannot express a
 *  link far side's identity without being wrong on most of its rows.
 *
 *  The pair signature DOES determine the far side, so it is the honest key.
 *  It is also long (~11k characters), so what is stored is a digest of it.
 *
 *  \section wk_compat Matching the Python side
 *
 *  cobordism-atlas/tools/frontier.py and identify_far_sides.py key the same
 *  table on `hashlib.sha1(pairsig.encode()).hexdigest()[:12]`, so this must
 *  agree with Python's hashlib byte for byte -- witnessKey() is checked
 *  against it in tests/witnesskey_test.cpp.
 *
 *  A prefix of the signature itself is deliberately NOT used: pair
 *  signatures share long prefixes, and truncating one silently merged
 *  distinct witnesses once.
 */

namespace witnesskey {

/** The full 40-character lowercase hex SHA-1 digest of \a data. */
std::string sha1Hex(const std::string &data);

/**
 * The key a far-side resolution for this witness is stored under: the
 * first 12 hex characters of the pair signature's SHA-1.
 *
 * 12 hex characters is 48 bits. Across the ~19k witnesses recorded so far
 * the birthday probability of any collision is ~6e-10, and a collision
 * would have to additionally survive the component-count check in
 * applyFarSideResolutions() to do any harm.
 */
std::string witnessKey(const std::string &pairSig);

} // namespace witnesskey

#endif

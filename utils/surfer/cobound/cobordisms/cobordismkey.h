//
//  cobordismkey.h
//
//  Created by John Teague on 09/08/2026.
//

#ifndef SURFER_COBOUND_COBORDISMKEY_H

#define SURFER_COBOUND_COBORDISMKEY_H

#include <string>

/*! \file utils/surfer/cobound/cobordisms/cobordismkey.h
 *  \brief The per-cobordism key an outgoing resolution is stored under.
 *
 *  \section wk_why Why a cobordism key rather than a name
 *
 *  An outgoing link reaches cobordisms.csv as a NAME produced by
 *  census::nameComplement() from
 *  the complement alone. For a knot that is enough -- Gordon-Luecke makes
 *  the complement determine the knot up to mirroring, and g_4 is
 *  mirror-invariant. For a LINK it is not: a complement does not determine
 *  a link (Rolfsen twisting), and one observed name is genuinely different
 *  links on different cobordisms -- "m129" alone stands for a dozen of them
 *  in our data. So a name-keyed table (--name-aliases) cannot express a
 *  link's identity without being wrong on most of its cobordisms.
 *
 *  The pair signature DOES determine the outgoing link, so it is the honest key.
 *  It is also long (~11k characters), so what is stored is a digest of it.
 *
 *  \section wk_compat Matching the Python side
 *
 *  cobordism-atlas/tools/frontier.py and identify_far_sides.py key the same
 *  table on `hashlib.sha1(pairsig.encode()).hexdigest()[:12]`, so this must
 *  agree with Python's hashlib byte for byte -- cobordismKey() is checked
 *  against it in tests/cobordismkey_test.cpp, and the SHA-1 itself
 *  (surfer/pairsig/sha1.h) in surfer's tests/sha1_test.cpp.
 *
 *  A prefix of the signature itself is deliberately NOT used: pair
 *  signatures share long prefixes, and truncating one silently merged
 *  distinct cobordisms once.
 */

namespace cobordisms {

/**
 * The key an outgoing resolution for this cobordism is stored under: the
 * first 12 hex characters of the pair signature's SHA-1.
 *
 * 12 hex characters is 48 bits. Across the ~19k cobordisms recorded so far
 * the birthday probability of any collision is ~6e-10, and a collision
 * would have to additionally survive the component-count check in
 * applyOutgoingResolutions() to do any harm.
 */
std::string cobordismKey(const std::string &pairSig);

} // namespace cobordisms

#endif

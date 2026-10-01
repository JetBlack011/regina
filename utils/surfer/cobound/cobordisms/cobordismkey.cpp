//
//  witnesskey.cpp
//
//  Created by John Teague on 09/08/2026.
//

#include "cobound/cobordisms/cobordismkey.h"

#include "surfer/pairsig/sha1.h"

namespace witnesskey {

std::string witnessKey(const std::string &pairSig) {
    return pairsig::sha1Hex(pairSig).substr(0, 12);
}

} // namespace witnesskey

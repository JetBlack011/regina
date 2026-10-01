//
//  sha1.cpp
//

#include "surfer/pairsig/sha1.h"

#include <openssl/sha.h>

namespace pairsig {

std::string sha1Hex(const std::string &data) {
    unsigned char digest[SHA_DIGEST_LENGTH];
    SHA1(reinterpret_cast<const unsigned char *>(data.data()), data.size(),
         digest);
    static constexpr char hex[] = "0123456789abcdef";
    std::string out(2 * SHA_DIGEST_LENGTH, '0');
    for (int i = 0; i < SHA_DIGEST_LENGTH; ++i) {
        out[2 * i] = hex[digest[i] >> 4];
        out[2 * i + 1] = hex[digest[i] & 0x0f];
    }
    return out;
}

} // namespace pairsig

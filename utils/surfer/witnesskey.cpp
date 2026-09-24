//
//  witnesskey.cpp
//
//  Created by John Teague on 09/08/2026.
//

#include "witnesskey.h"

#include <array>
#include <cstdint>
#include <cstdio>

namespace witnesskey {

namespace {

inline uint32_t rotl(uint32_t v, int n) { return (v << n) | (v >> (32 - n)); }

// FIPS 180-1. Written out rather than pulled from a dependency: the only
// hash needed anywhere in surfer is this one, and it must agree exactly
// with Python's hashlib.sha1 (see witnesskey.h).
void processBlock(const unsigned char *p, std::array<uint32_t, 5> &h) {
    uint32_t w[80];
    for (int i = 0; i < 16; ++i)
        w[i] = (uint32_t(p[i * 4]) << 24) | (uint32_t(p[i * 4 + 1]) << 16) |
               (uint32_t(p[i * 4 + 2]) << 8) | uint32_t(p[i * 4 + 3]);
    for (int i = 16; i < 80; ++i)
        w[i] = rotl(w[i - 3] ^ w[i - 8] ^ w[i - 14] ^ w[i - 16], 1);

    uint32_t a = h[0], b = h[1], c = h[2], d = h[3], e = h[4];
    for (int i = 0; i < 80; ++i) {
        uint32_t f, k;
        if (i < 20) {
            f = (b & c) | ((~b) & d);
            k = 0x5A827999u;
        } else if (i < 40) {
            f = b ^ c ^ d;
            k = 0x6ED9EBA1u;
        } else if (i < 60) {
            f = (b & c) | (b & d) | (c & d);
            k = 0x8F1BBCDCu;
        } else {
            f = b ^ c ^ d;
            k = 0xCA62C1D6u;
        }
        const uint32_t tmp = rotl(a, 5) + f + e + k + w[i];
        e = d;
        d = c;
        c = rotl(b, 30);
        b = a;
        a = tmp;
    }
    h[0] += a;
    h[1] += b;
    h[2] += c;
    h[3] += d;
    h[4] += e;
}

} // namespace

std::string sha1Hex(const std::string &data) {
    std::array<uint32_t, 5> h{0x67452301u, 0xEFCDAB89u, 0x98BADCFEu,
                              0x10325476u, 0xC3D2E1F0u};

    const auto *bytes = reinterpret_cast<const unsigned char *>(data.data());
    const size_t n = data.size();

    size_t full = n / 64;
    for (size_t i = 0; i < full; ++i)
        processBlock(bytes + i * 64, h);

    // Tail: the remaining bytes, the 0x80 terminator, zero padding, and the
    // 64-bit big-endian bit length. Two blocks are needed exactly when the
    // remainder leaves no room for the length field.
    unsigned char tail[128] = {0};
    const size_t rem = n - full * 64;
    for (size_t i = 0; i < rem; ++i)
        tail[i] = bytes[full * 64 + i];
    tail[rem] = 0x80;
    const size_t tailBlocks = (rem >= 56) ? 2 : 1;
    const uint64_t bits = uint64_t(n) * 8;
    for (int i = 0; i < 8; ++i)
        tail[tailBlocks * 64 - 1 - i] =
            static_cast<unsigned char>((bits >> (8 * i)) & 0xFF);
    for (size_t i = 0; i < tailBlocks; ++i)
        processBlock(tail + i * 64, h);

    char out[41];
    for (int i = 0; i < 5; ++i)
        std::snprintf(out + i * 8, 9, "%08x", h[i]);
    return std::string(out, 40);
}

std::string witnessKey(const std::string &pairSig) {
    return sha1Hex(pairSig).substr(0, 12);
}

} // namespace witnesskey

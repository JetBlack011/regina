//
//  sha1_test.cpp
//
//  pairsig::sha1Hex() (OpenSSL's SHA1()) against FIPS 180-1's own vectors,
//  the length-boundary cases a hand-written padding step gets wrong, and the
//  values Python's hashlib gives (which the atlas keys cobordisms on). Before
//  it replaced the hand-written SHA-1 (witnesskey.cpp, retired in phase 2 of
//  the refactor), this test also compared the two on every length from 0 to
//  4200 bytes of pseudo-random data and on inputs up to 16 MiB: no input
//  differed.
//

#include <iostream>
#include <string>

#include "surfer/pairsig/sha1.h"

namespace {

int failures = 0;

void check(const std::string &what, const std::string &got,
           const std::string &want) {
    if (got == want) {
        std::cout << "  ok: " << what << "\n";
    } else {
        std::cout << "  FAIL: " << what << "\n    got  " << got
                  << "\n    want " << want << "\n";
        ++failures;
    }
}

} // namespace

int main() {
    std::cout << "sha1\n";

    // FIPS 180-1 appendix vectors.
    check("empty string", pairsig::sha1Hex(""),
          "da39a3ee5e6b4b0d3255bfef95601890afd80709");
    check("\"abc\"", pairsig::sha1Hex("abc"),
          "a9993e364706816aba3e25717850c26c9cd0d89d");
    check("56-byte message (two-block tail)",
          pairsig::sha1Hex(
              "abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq"),
          "84983e441c3bd26ebaae4aa1f95129e5e54670f1");
    check("one million 'a'", pairsig::sha1Hex(std::string(1000000, 'a')),
          "34aa973cd4c4daa4f61eeb2bdbad27316534016f");

    // Length-boundary cases: 55 bytes still fits one padding block, 56 and 64
    // force a second. These are exactly where a hand-written tail goes wrong.
    check("55 bytes", pairsig::sha1Hex(std::string(55, 'a')),
          "c1c8bbdc22796e28c0e15163d20899b65621d65a");
    check("56 bytes", pairsig::sha1Hex(std::string(56, 'a')),
          "c2db330f6083854c99d4b5bfb6e8f29f201be699");
    check("64 bytes", pairsig::sha1Hex(std::string(64, 'a')),
          "0098ba824b5c16427bd7a1122a5a442a25ec644d");
    check("1000 bytes", pairsig::sha1Hex(std::string(1000, 'a')),
          "291e9a6c66994949b57ba5e650361e98fc36b1ba");

    if (failures) {
        std::cout << failures << " failure(s)\n";
        return 1;
    }
    std::cout << "all passed\n";
    return 0;
}

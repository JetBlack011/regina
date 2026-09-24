//
//  witnesskey_test.cpp
//
//  Checks witnesskey::sha1Hex against known FIPS 180-1 vectors and against
//  values produced by Python's hashlib, which is what the cobordism-atlas
//  tooling keys far-side resolutions on. If these two ever disagree, a
//  resolution written by the Python side silently fails to match any
//  witness on the C++ side -- and vice versa -- with no error anywhere.
//

#include <iostream>
#include <string>

#include "witnesskey.h"

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
    std::cout << "witnesskey\n";

    // FIPS 180-1 appendix vectors.
    check("empty string", witnesskey::sha1Hex(""),
          "da39a3ee5e6b4b0d3255bfef95601890afd80709");
    check("\"abc\"", witnesskey::sha1Hex("abc"),
          "a9993e364706816aba3e25717850c26c9cd0d89d");
    check("56-byte message (two-block tail)",
          witnesskey::sha1Hex(
              "abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq"),
          "84983e441c3bd26ebaae4aa1f95129e5e54670f1");

    // Length-boundary cases: 55 bytes still fits one padding block, 56 and 64
    // force a second. These are exactly where a hand-written tail goes wrong.
    check("55 bytes", witnesskey::sha1Hex(std::string(55, 'a')),
          "c1c8bbdc22796e28c0e15163d20899b65621d65a");
    check("56 bytes", witnesskey::sha1Hex(std::string(56, 'a')),
          "c2db330f6083854c99d4b5bfb6e8f29f201be699");
    check("64 bytes", witnesskey::sha1Hex(std::string(64, 'a')),
          "0098ba824b5c16427bd7a1122a5a442a25ec644d");
    check("1000 bytes", witnesskey::sha1Hex(std::string(1000, 'a')),
          "291e9a6c66994949b57ba5e650361e98fc36b1ba");

    // A real pair signature prefix, keyed the way frontier.py keys it.
    check("witnessKey is the 12-char prefix",
          witnesskey::witnessKey("abc"), "a9993e364706");

    if (failures) {
        std::cout << failures << " failure(s)\n";
        return 1;
    }
    std::cout << "all passed\n";
    return 0;
}

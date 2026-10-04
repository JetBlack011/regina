//
//  cobordismkey_test.cpp
//
//  Checks cobordisms::cobordismKey() against values produced by Python's
//  hashlib, which is what the cobordism-atlas tooling keys outgoing
//  link resolutions on. If these two ever disagree, a resolution written by the
//  Python side silently fails to match any cobordism on the C++ side -- and
//  vice versa -- with no error anywhere. The SHA-1 underneath
//  (surfer/pairsig/sha1.h) is checked against FIPS 180-1 in surfer's
//  tests/sha1_test.cpp.
//

#include <iostream>
#include <string>

#include "cobound/cobordisms/cobordismkey.h"

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

    // A real pair signature prefix, keyed the way frontier.py keys it.
    check("witnessKey is the 12-char prefix",
          cobordisms::cobordismKey("abc"), "a9993e364706");
    check("the empty signature's key", cobordisms::cobordismKey(""),
          "da39a3ee5e6b");
    check("a 1000-byte signature's key",
          cobordisms::cobordismKey(std::string(1000, 'a')), "291e9a6c6699");

    if (failures) {
        std::cout << failures << " failure(s)\n";
        return 1;
    }
    std::cout << "all passed\n";
    return 0;
}

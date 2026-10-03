//
//  cobordism.cpp
//

#include "cobound/cobordisms/cobordism.h"

#include <algorithm>

namespace cobordisms {

std::string witnessIdentity(const Witness &w) {
    // Unit separators: no field (names, keys, numbers) can contain one.
    constexpr char SEP = '\x1f';
    std::string key;
    key.reserve(w.subject.size() + w.other.size() + 32);
    key += w.kind == WitnessKind::direct ? 'd' : 'c';
    key += SEP;
    key += w.subject;
    key += SEP;
    key += std::to_string(w.subjectComponents);
    key += SEP;
    key += w.other;
    key += SEP;
    key += std::to_string(w.otherComponents);
    key += SEP;
    key += std::to_string(w.genus);
    key += SEP;
    key += w.tubed ? 't' : 'f';
    key += SEP;
    key += w.resolvedVertices > 0 ? 'r' : 'e';
    return key;
}

bool haveWitness(const std::vector<Witness> &witnesses, const Witness &w) {
    const std::string key = witnessIdentity(w);
    return std::any_of(witnesses.begin(), witnesses.end(),
                       [&key](const Witness &existing) {
                           return witnessIdentity(existing) == key;
                       });
}
} // namespace cobordisms

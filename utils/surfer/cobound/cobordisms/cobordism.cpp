//
//  cobordism.cpp
//

#include "cobound/cobordisms/cobordism.h"

#include <algorithm>

namespace cobordisms {

std::string cobordismIdentity(const Cobordism &w) {
    // Unit separators: no field (names, keys, numbers) can contain one.
    constexpr char SEP = '\x1f';
    std::string key;
    key.reserve(w.subject.size() + w.other.size() + 32);
    key += w.kind == CobordismKind::direct ? 'd' : 'c';
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

bool haveCobordism(const std::vector<Cobordism> &cobordisms, const Cobordism &w) {
    const std::string key = cobordismIdentity(w);
    return std::any_of(cobordisms.begin(), cobordisms.end(),
                       [&key](const Cobordism &existing) {
                           return cobordismIdentity(existing) == key;
                       });
}
} // namespace cobordisms

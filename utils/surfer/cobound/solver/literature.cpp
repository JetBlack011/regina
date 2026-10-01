//
//  literature.cpp
//

#include "cobound/solver/literature.h"

#include <map>
#include <tuple>

#include "linknaming/names.h"

namespace cobordismgraph {

void NameTable::addLiterature(const std::string &name, int lo, int hi) {
    auto [it, inserted] = info_.try_emplace(name);
    if (inserted) {
        it->second.components = componentsFromName(name);
        byBase_[baseName(name)].push_back(name);
    }
    it->second.haveLiterature = true;
    it->second.litLo = lo;
    it->second.litHi = hi;
}

const NameInfo *NameTable::find(const std::string &name) const {
    auto it = info_.find(name);
    return it == info_.end() ? nullptr : &it->second;
}

int NameTable::components(const std::string &name) const {
    auto it = info_.find(name);
    return it == info_.end() ? componentsFromName(name) : it->second.components;
}

std::vector<std::string>
NameTable::candidates(const std::string &name,
                      std::optional<int> observedComponents) const {
    auto it = byBase_.find(baseName(name));
    if (it == byBase_.end() || it->second.empty())
        return {name};

    if (!observedComponents)
        return it->second;

    std::vector<std::string> filtered;
    for (const std::string &candidate : it->second)
        if (components(candidate) == *observedComponents)
            filtered.push_back(candidate);
    // An empty result means every registered variant of this base has the
    // wrong component count, i.e. the identification and the geometry
    // disagree. The far side is none of those variants, so a max/min over
    // them bounds nothing: hand back the name itself, which is unregistered
    // and so bears no bound.
    return filtered.empty() ? std::vector<std::string>{name} : filtered;
}

namespace {

/**
 * Composite knots that bound a smooth disk, and so ground a chain exactly as
 * the unknot does.
 *
 * `K # m(K^r)` -- K summed with the reverse of its mirror -- is the identity
 * of the concordance group and bounds an explicit ribbon disk, for every K.
 * Which *spelling* of that qualifies depends on the symmetry of K, so this is
 * an explicit allowlist and not a pattern:
 *
 *   - "3_1#m3_1": 3_1 is invertible, so 3_1^r = 3_1 and m(3_1^r) = m3_1.
 *   - "4_1#4_1":  4_1 is invertible AND amphichiral, so m(4_1^r) = 4_1 and
 *                 the sum of 4_1 with *itself* is already the ribbon case.
 *
 * DO NOT generalize this to the string pattern "A#mA". That is wrong as soon
 * as A is non-invertible (the first such knot is 8_17), where A # mA and
 * A # m(A^r) are different knots and only the latter is slice. Adding an
 * entry means checking the symmetry of the summand first.
 *
 * Names are matched exactly, in the canonical spelling that
 * tools/identify_by_retriangulation.py emits: summands sorted, and the
 * lexicographically smaller of the name and its overall mirror.
 */
bool isSliceComposite(const std::string &name) {
    return name == "3_1#m3_1" || name == "4_1#4_1";
}
} // namespace

bool isElementarySlice(const std::string &name, const NameTable &names) {
    if (isSliceComposite(name))
        return true;
    std::vector<std::string> parts = knotSummands(name);
    if (parts.empty())
        return false;
    // Each summand reduced to what its symmetry type leaves meaningful of its
    // marks (m: mirrored, r: reversed), as (knot, m, r) with the meaningless
    // marks cleared; its concordance inverse -K = mrK reduced the same way.
    // Slice when every reduced summand pairs off with its inverse (a
    // self-inverse one with another copy of itself).
    using Key = std::tuple<std::string, bool, bool>;
    auto reduce = [](const std::string &knot, SymmetryType type, bool m, bool r) -> Key {
        switch (type) {
        case SymmetryType::fullyAmphicheiral: return {knot, false, false};
        case SymmetryType::reversible: return {knot, m, false};         // K^r = K
        case SymmetryType::negativeAmphicheiral: return {knot, m != r, false}; // K^r = mK
        case SymmetryType::positiveAmphicheiral: return {knot, false, r};      // mK = K
        default: return {knot, m, r};                                   // chiral
        }
    };
    std::map<Key, int> count;
    std::map<Key, Key> inverse;
    for (const std::string &p : parts) {
        if (p == "Unknot" || p == "mUnknot")
            continue;
        size_t k = 0;
        const bool m = k < p.size() && p[k] == 'm';
        if (m) ++k;
        const bool r = k < p.size() && p[k] == 'r';
        if (r) ++k;
        const std::string knot = p.substr(k);
        const SymmetryType *type = names.symmetry(knot);
        if (!type)
            return false;
        const Key key = reduce(knot, *type, m, r);
        ++count[key];
        inverse[key] = reduce(knot, *type, !m, !r);
    }
    for (const auto &[key, n] : count) {
        const Key &inv = inverse.at(key);
        if (inv == key) {
            if (n % 2 != 0)
                return false;
        } else {
            auto it = count.find(inv);
            if (it == count.end() || it->second != n)
                return false;
        }
    }
    return true;
}

std::optional<SymmetryType> parseSymmetryType(const std::string &text) {
    if (text == "chiral")
        return SymmetryType::chiral;
    if (text == "reversible")
        return SymmetryType::reversible;
    if (text == "positive amphicheiral")
        return SymmetryType::positiveAmphicheiral;
    if (text == "negative amphicheiral")
        return SymmetryType::negativeAmphicheiral;
    if (text == "fully amphicheiral")
        return SymmetryType::fullyAmphicheiral;
    return std::nullopt;
}
} // namespace cobordismgraph

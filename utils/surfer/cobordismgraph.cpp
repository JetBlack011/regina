//
//  cobordismgraph.cpp
//
//  Created by John Teague on 08/02/2026.
//

#include "cobordismgraph.h"

#include <algorithm>
#include <cctype>
#include <sstream>
#include <unordered_set>

namespace cobordismgraph {

/* Names of knots/links that appear in the cobordism graph */

// TODO: Probably these shouldn't be here... Might make sense to move them to
// their own file at some point.

const std::string kSplitSeparator = " u ";

std::string baseName(const std::string &name) {
    // A SPLIT name has no base in this sense. Stripping at the first '{'
    // would turn "L2a1{0} u Unknot" into "L2a1", so NameTable::candidates()
    // would hand back L2a1's orientation variants as if the far side were
    // that two-component link rather than a three-component split one --
    // enumerating variants of a summand as variants of the whole. Returning
    // the name unchanged makes the byBase_ lookup miss, which is exactly
    // right: a split name's alternatives live in its FACTORS, and upperOf()/
    // lowerOf() take them from there.
    if (name.find(kSplitSeparator) != std::string::npos)
        return name;
    size_t brace = name.find('{');
    return brace == std::string::npos ? name : name.substr(0, brace);
}

std::string normalizeIdentifiedName(const std::string &name) {
    if (name.size() < 4 || name.back() != ')')
        return name;
    size_t open = name.rfind(" (");
    if (open == std::string::npos || open == 0)
        return name;
    return name.substr(0, open);
}

/* Split (disjoint-union) far sides */

// The separator decompose_far_sides.py writes between the factors of a split
// link, as recorded in results/split_far_sides.csv: "3_1 u Unknot".

// The factors of a split name, or an empty vector if `name` is not one.
//
// A far side is split exactly when its exterior is reducible, which is the
// commonest thing an unnameable multi-component far side turns out to be:
// measured over 70 sampled unidentified link far sides, 43 of them. The
// factors are named separately, each by its own exterior, and the link is
// their disjoint union.
std::vector<std::string> splitFactors(const std::string &name) {
    std::vector<std::string> factors;
    size_t pos = 0;
    for (;;) {
        size_t sep = name.find(kSplitSeparator, pos);
        if (sep == std::string::npos)
            break;
        factors.push_back(name.substr(pos, sep - pos));
        pos = sep + kSplitSeparator.size();
    }
    if (factors.empty())
        return {};
    factors.push_back(name.substr(pos));
    return factors;
}

std::vector<std::string> factorAlternatives(const std::string &factor) {
    std::vector<std::string> alts;
    size_t pos = 0;
    for (;;) {
        size_t bar = factor.find('|', pos);
        if (bar == std::string::npos)
            break;
        alts.push_back(factor.substr(pos, bar - pos));
        pos = bar + 1;
    }
    alts.push_back(factor.substr(pos));
    return alts;
}

std::optional<CompositeName> compositeParts(const std::string &name) {
    // "<knot> #_<c> <link>", knot = m?<digits>_<digits>,
    // link = L<digits><a|n><digits>. Anything else is not composite.
    const std::string sep = " #_";
    const size_t at = name.find(sep);
    if (at == std::string::npos || at == 0)
        return std::nullopt;
    const std::string knot = name.substr(0, at);
    const size_t digits = at + sep.size();
    const size_t space = name.find(' ', digits);
    if (space == std::string::npos || space == digits)
        return std::nullopt;
    const std::string comp = name.substr(digits, space - digits);
    const std::string link = name.substr(space + 1);
    auto allDigits = [](const std::string &t, size_t from, size_t to) {
        if (from >= to)
            return false;
        for (size_t i = from; i < to; ++i)
            if (!std::isdigit(static_cast<unsigned char>(t[i])))
                return false;
        return true;
    };
    const size_t k0 = (!knot.empty() && knot[0] == 'm') ? 1 : 0;
    const size_t us = knot.find('_', k0);
    if (us == std::string::npos)
        return std::nullopt;
    const size_t de =
        (us > k0 && (knot[us - 1] == 'a' || knot[us - 1] == 'n')) ? us - 1 : us;
    if (!allDigits(knot, k0, de) || !allDigits(knot, us + 1, knot.size()))
        return std::nullopt;
    if (!allDigits(comp, 0, comp.size()))
        return std::nullopt;
    size_t i = 1;
    if (link.empty() || link[0] != 'L')
        return std::nullopt;
    while (i < link.size() && std::isdigit(static_cast<unsigned char>(link[i])))
        ++i;
    if (i == 1 || i >= link.size() || (link[i] != 'a' && link[i] != 'n') ||
        !allDigits(link, i + 1, link.size()))
        return std::nullopt;
    return CompositeName{knot, std::stoi(comp), link};
}

std::vector<std::string> knotSummands(const std::string &name) {
    if (name.find(" : ") != std::string::npos ||
        name.find(" #_") != std::string::npos ||
        name.find(kSplitSeparator) != std::string::npos ||
        name.find('{') != std::string::npos)
        return {};
    std::vector<std::string> parts;
    size_t pos = 0;
    for (;;) {
        size_t hash = name.find('#', pos);
        parts.push_back(name.substr(pos, hash == std::string::npos
                                             ? std::string::npos
                                             : hash - pos));
        if (hash == std::string::npos)
            break;
        pos = hash + 1;
    }
    if (parts.size() < 2)
        return {};
    for (const std::string &p : parts) {
        if (p == "Unknot" || p == "mUnknot")
            continue;
        // m? <digits> [a|n]? _ <digits>: "5_2", "m5_2", "11a_367".
        const size_t k0 = (!p.empty() && p[0] == 'm') ? 1 : 0;
        const size_t us = p.find('_', k0);
        if (us == std::string::npos || us == k0 || us + 1 >= p.size())
            return {};
        size_t digitsEnd = us;
        if (p[us - 1] == 'a' || p[us - 1] == 'n')
            --digitsEnd;
        if (digitsEnd == k0)
            return {};
        for (size_t i = k0; i < p.size(); ++i)
            if (i != us && i != digitsEnd &&
                !std::isdigit(static_cast<unsigned char>(p[i])))
                return {};
        if (digitsEnd != us && (p[digitsEnd] != 'a' && p[digitsEnd] != 'n'))
            return {};
    }
    return parts;
}

int componentsFromName(const std::string &name) {
    // A split name's components are its factors' added up. Without this a far
    // side recorded as "3_1 u Unknot" would claim ONE component, and every
    // check matching a name's component count against the curve count actually
    // observed on that boundary would reject it.
    //
    // LIMIT: an UNTAGGED link base name does not state its own count, so
    // "Unknot u L2a1" reads as 2 rather than 3; only "Unknot u L2a1{0}" is
    // counted correctly. Emitters must not write a split name with an
    // untagged link factor. frontier.py's components_from_name() carries the
    // identical limit on purpose, so --check still compares like with like;
    // fixing it properly needs the name TABLE rather than the name.
    if (std::vector<std::string> factors = splitFactors(name);
        !factors.empty()) {
        int total = 0;
        for (const std::string &f : factors)
            total += componentsFromName(factorAlternatives(f).front());
        return total;
    }

    // "<n>-component unlink"
    if (name.ends_with("-component unlink")) {
        size_t dash = name.find('-');
        if (dash != std::string::npos && dash > 0) {
            bool allDigits = std::all_of(
                name.begin(), name.begin() + static_cast<long>(dash),
                [](unsigned char c) { return std::isdigit(c) != 0; });
            if (allDigits) {
                try {
                    return std::stoi(name.substr(0, dash));
                } catch (const std::exception &) {
                    // fall through to the default below
                }
            }
        }
        return 1;
    }

    // An orientation tag lists a choice per component after the first, so
    // the component count is one more than the number of entries. See
    // componentsFromName()'s doc comment.
    size_t brace = name.find('{');
    if (brace == std::string::npos)
        return 1;
    size_t close = name.find('}', brace);
    if (close == std::string::npos)
        return 1;
    int entries = 1;
    for (size_t i = brace + 1; i < close; ++i)
        if (name[i] == ';')
            ++entries;
    return entries + 1;
}

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
    // An empty result would mean every registered variant of this base has
    // the wrong component count, i.e. the identification and the geometry
    // disagree. Hand back the unfiltered set and let the max/min in propagate()
    // stay okay.
    return filtered.empty() ? it->second : filtered;
}

/* Witness cobordisms (cobordisms that verify slice genus somehow) */

bool haveWitness(const std::vector<Witness> &witnesses, const Witness &w) {
    return std::any_of(
        witnesses.begin(), witnesses.end(), [&w](const Witness &existing) {
            return existing.kind == w.kind && existing.subject == w.subject &&
                   existing.other == w.other && existing.genus == w.genus &&
                   existing.tubed == w.tubed &&
                   existing.otherComponents == w.otherComponents &&
                   existing.subjectComponents == w.subjectComponents;
        });
}

bool farSideBearsBound(const Witness &w) {
    // The observed count, never the name: componentsFromName() is an
    // inference from a string, and an alias could make a two-curve far side
    // read like a knot.
    // A far side proved per witness -- meridians carried, or an all-knot
    // split -- is the one case where a multi-component name is a complete
    // candidate set; see Witness::farSideProved.
    return w.otherComponents == 1 || w.farSideProved ||
           identify::isOrientationSafeName(w.other);
}

namespace {

/** Saturating add that keeps NO_UPPER_BOUND */
int addUpper(int bound, int delta) {
    if (bound == NO_UPPER_BOUND)
        return NO_UPPER_BOUND;
    return bound + delta;
}

/** One candidate's contribution to an upper bound, together with each name
 * whose literature value that contribution rests on. */
struct UpperContribution {
    int value = NO_UPPER_BOUND;
    std::vector<std::string> support;
};

/** Sorted union: cheaper than keeping a set for a small list of names. */
void mergeSupport(std::vector<std::string> &into,
                  const std::vector<std::string> &from) {
    into.insert(into.end(), from.begin(), from.end());
    std::sort(into.begin(), into.end());
    into.erase(std::unique(into.begin(), into.end()), into.end());
}

bool supportContains(const std::vector<std::string> &support,
                     const std::string &name) {
    return std::binary_search(support.begin(), support.end(), name);
}

namespace {

// A prime knot's name, possibly mirrored: "3_1", "m3_1", "11n_34".
bool isPrimeKnotName(const std::string &alt) {
    const std::string s = (alt.size() > 1 && alt[0] == 'm') ? alt.substr(1) : alt;
    size_t i = 0;
    while (i < s.size() && std::isdigit(static_cast<unsigned char>(s[i])))
        ++i;
    if (i == 0)
        return false;
    if (i < s.size() && (s[i] == 'a' || s[i] == 'n'))
        ++i;
    if (i >= s.size() || s[i] != '_')
        return false;
    size_t j = ++i;
    while (j < s.size() && std::isdigit(static_cast<unsigned char>(s[j])))
        ++j;
    return j > i && j == s.size();
}

// Whether one alternative of a split factor is a KNOT: the Unknot, a table
// knot name possibly mirrored ("3_1", "m3_1", "11n_34"), or a composite knot
// ("3_1#5_1"). Only knot factors take part in the split rule's lower bound.
bool isKnotFactor(const std::string &alt) {
    return alt == "Unknot" || !knotSummands(alt).empty() ||
           isPrimeKnotName(alt);
}

// A prime knot's mirror image has the same g_4, and the name table holds only
// the unmirrored name, so a bound for "m3_1" is 3_1's. A COMPOSITE is left
// alone: stripping the leading m of "m3_1#3_1" (the square knot, slice)
// would give "3_1#3_1" (the granny, g_4 = 2), a different knot. Its summands
// are unmirrored where the composite rule reads them.
std::string unmirrored(const std::string &alt) {
    return (alt.size() > 1 && alt[0] == 'm' && isPrimeKnotName(alt))
               ? alt.substr(1)
               : alt;
}

} // namespace

UpperContribution upperOf(const std::string &name,
                          const std::unordered_map<std::string, Bounds> &bounds,
                          const NameTable &names) {
    UpperContribution best;
    // Both sources are valid upper bounds, so take whichever is better.
    auto it = bounds.find(name);
    if (it != bounds.end() && it->second.haveUpper())
        best = {.value = it->second.hi, .support = it->second.support};
    if (const NameInfo *info = names.find(name); info && info->haveLiterature)
        if (best.value == NO_UPPER_BOUND || info->litHi < best.value)
            best = {.value = info->litHi, .support = {name}};

    // A split far side is not in any table, so its UPPER bound comes from its
    // factors: g_4(A u B) <= g_4(A) + g_4(B), tubing the factors' minimal
    // surfaces together in B^4 (a tube leaves the genus the sum) --
    // constructive. There is no matching lower bound by addition: see
    // lowerOf(), where K u -K bounding an annulus is the reason.
    if (std::vector<std::string> factors = splitFactors(name);
        !factors.empty()) {
        int total = 0;
        std::vector<std::string> support;
        bool haveAll = true;
        for (const std::string &f : factors) {
            // Worst case over a factor's alternatives, the same convention
            // candidates() uses: the bound has to hold whichever it is.
            int worstAlt = NO_UPPER_BOUND;
            std::vector<std::string> altSupport;
            for (const std::string &alt : factorAlternatives(f)) {
                UpperContribution u = upperOf(
                    isKnotFactor(alt) ? unmirrored(alt) : alt, bounds, names);
                if (u.value == NO_UPPER_BOUND) {
                    worstAlt = NO_UPPER_BOUND;
                    break;
                }
                worstAlt = worstAlt == NO_UPPER_BOUND
                               ? u.value
                               : std::max(worstAlt, u.value);
                mergeSupport(altSupport, u.support);
            }
            if (worstAlt == NO_UPPER_BOUND) {
                haveAll = false;
                break;
            }
            total += worstAlt;
            mergeSupport(support, altSupport);
        }
        if (haveAll && (best.value == NO_UPPER_BOUND || total < best.value))
            best = {.value = total, .support = std::move(support)};
    }
    // A composite K #_c L: g_4(K # L) <= g_4(K) + g_4(L), boundary-connect-
    // summing the two minimal surfaces in B^4 # B^4 = B^4 -- constructive,
    // and true for every component c. (The matching lower bound,
    // g_4(L) - g_4(K), is in lowerOf.) L is a base name proved up to
    // orientation, so take the WORST of its registered orientation variants;
    // if it has none, no bound.
    if (std::optional<CompositeName> cp = compositeParts(name)) {
        const std::string knot =
            cp->knot[0] == 'm' ? cp->knot.substr(1) : cp->knot;
        UpperContribution k = upperOf(knot, bounds, names);
        const std::vector<std::string> variants = names.candidates(cp->link);
        const bool registered = !variants.empty() && variants.front() != cp->link;
        if (k.value != NO_UPPER_BOUND && registered) {
            int worst = NO_UPPER_BOUND;
            std::vector<std::string> support = k.support;
            bool haveAll = true;
            for (const std::string &v : variants) {
                UpperContribution u = upperOf(v, bounds, names);
                if (u.value == NO_UPPER_BOUND) {
                    haveAll = false;
                    break;
                }
                worst = worst == NO_UPPER_BOUND ? u.value : std::max(worst, u.value);
                mergeSupport(support, u.support);
            }
            if (haveAll && worst != NO_UPPER_BOUND &&
                (best.value == NO_UPPER_BOUND || k.value + worst < best.value))
                best = {.value = k.value + worst, .support = std::move(support)};
        }
    }
    // A composite KNOT A#B#...: g_4 is subadditive under connected sum, so
    // g_4 <= sum of the summands' g_4 (mirrors look up as their knot, g_4
    // being mirror-invariant).
    if (std::vector<std::string> parts = knotSummands(name); !parts.empty()) {
        int total = 0;
        std::vector<std::string> support;
        bool haveAll = true;
        for (const std::string &p : parts) {
            if (p == "Unknot" || p == "mUnknot")
                continue;
            UpperContribution u =
                upperOf(p[0] == 'm' ? p.substr(1) : p, bounds, names);
            if (u.value == NO_UPPER_BOUND) {
                haveAll = false;
                break;
            }
            total += u.value;
            mergeSupport(support, u.support);
        }
        if (haveAll && (best.value == NO_UPPER_BOUND || total < best.value))
            best = {.value = total, .support = std::move(support)};
    }
    return best;
}

/** As UpperContribution, for lower bounds. */
struct LowerContribution {
    int value = NO_LOWER_BOUND;
    std::vector<std::string> support;
};

LowerContribution lowerOf(const std::string &name,
                          const std::unordered_map<std::string, Bounds> &bounds,
                          const NameTable &names) {
    LowerContribution best;
    if (const NameInfo *info = names.find(name); info && info->haveLiterature)
        best = {.value = info->litLo, .support = {name}};
    auto it = bounds.find(name);
    if (it != bounds.end() && it->second.haveLower() &&
        (best.value == NO_LOWER_BOUND || it->second.lo > best.value))
        best = {.value = it->second.lo, .support = it->second.lowerSupport};

    // The split rule in the LOWER direction. It is NOT additive: for every
    // knot K the split link K u -K (-K the mirror reverse) bounds an annulus
    // in B^4, so g_4(4_1 u 4_1) = 0 although each factor has g_4 = 1. (The
    // additive version, adopted as an assumption on 2026-09-14, produced
    // "lower bound 1, ABOVE the literature upper bound 0" for L10a91{0}, whose
    // genus-0 witnesses end at 4_1 u 4_1.) What does hold, by one band either
    // way -- a band joining two components of a connected surface raises its
    // genus by one, a band splitting one component leaves it unchanged -- is
    //     g_4(#factors) - (f - 1) <= g_4(u factors) <= g_4(#factors),
    // after dropping split UNKNOT factors, which change nothing:
    // g_4(L u U) = g_4(L), since a split unknot is capped off by a disc in a
    // collar (or tubed on). g_4 of the sum is bounded below by the
    // composite-knot rule, g_4(K_i) - sum_{j != i} g_4(K_j), so this applies
    // only when every remaining factor is a knot. A factor's alternatives
    // ("3_1|m3_1") must all satisfy it: the MIN of their lower bounds and
    // the MAX of their upper bounds.
    if (std::vector<std::string> factors = splitFactors(name);
        !factors.empty()) {
        std::vector<std::vector<std::string>> knots;   // nontrivial factors
        bool allKnots = true;
        for (const std::string &f : factors) {
            std::vector<std::string> alts;
            bool trivial = true;
            for (const std::string &alt : factorAlternatives(f)) {
                if (!isKnotFactor(alt)) {
                    allKnots = false;
                    break;
                }
                alts.push_back(unmirrored(alt));
                if (alt != "Unknot")
                    trivial = false;
            }
            if (!allKnots)
                break;
            if (!trivial)
                knots.push_back(std::move(alts));
        }
        if (allKnots && !knots.empty()) {
            const int f = static_cast<int>(knots.size());
            for (int i = 0; i < f; ++i) {
                int value = NO_LOWER_BOUND;
                std::vector<std::string> support;
                bool ok = true;
                for (const std::string &alt : knots[i]) {
                    LowerContribution l = lowerOf(alt, bounds, names);
                    if (l.value == NO_LOWER_BOUND) {
                        ok = false;
                        break;
                    }
                    value = value == NO_LOWER_BOUND ? l.value
                                                    : std::min(value, l.value);
                    mergeSupport(support, l.support);
                }
                for (int j = 0; ok && j < f; ++j) {
                    if (j == i)
                        continue;
                    int worst = NO_UPPER_BOUND;
                    for (const std::string &alt : knots[j]) {
                        UpperContribution u = upperOf(alt, bounds, names);
                        if (u.value == NO_UPPER_BOUND) {
                            ok = false;
                            break;
                        }
                        worst = worst == NO_UPPER_BOUND ? u.value
                                                        : std::max(worst, u.value);
                        mergeSupport(support, u.support);
                    }
                    if (ok)
                        value -= worst;
                }
                if (!ok)
                    continue;
                value -= f - 1;
                if (best.value == NO_LOWER_BOUND || value > best.value)
                    best = {.value = value, .support = std::move(support)};
            }
        }
    }
    // A composite KNOT: K_i is concordant to (A # -rest), so
    //     g_4(K_i) <= g_4(A) + sum_{j != i} g_4(K_j),
    // i.e. g_4(A) >= g_4(K_i) - sum_{j != i} g_4(K_j), for every i. This is
    // the only lower bound a sum admits in terms of its summands: K # -K is
    // slice, so nothing like g_4(K) + g_4(J) - c holds.
    if (std::vector<std::string> parts = knotSummands(name); !parts.empty()) {
        std::vector<std::string> knots;
        for (const std::string &p : parts)
            if (p != "Unknot" && p != "mUnknot")
                knots.push_back(p[0] == 'm' ? p.substr(1) : p);
        for (size_t i = 0; i < knots.size(); ++i) {
            LowerContribution l = lowerOf(knots[i], bounds, names);
            if (l.value == NO_LOWER_BOUND)
                continue;
            int value = l.value;
            std::vector<std::string> support = l.support;
            bool haveAll = true;
            for (size_t j = 0; j < knots.size(); ++j) {
                if (j == i)
                    continue;
                UpperContribution u = upperOf(knots[j], bounds, names);
                if (u.value == NO_UPPER_BOUND) {
                    haveAll = false;
                    break;
                }
                value -= u.value;
                mergeSupport(support, u.support);
            }
            if (haveAll && (best.value == NO_LOWER_BOUND || value > best.value))
                best = {.value = value, .support = std::move(support)};
        }
    }

    // A composite K #_c L: summing -K into the same component undoes K up to
    // concordance (K # -K is slice, and summing a slice knot into a component
    // gives a concordant link), so g_4(L) <= g_4(K #_c L) + g_4(K), i.e.
    //     g_4(K #_c L) >= g_4(L) - g_4(K).
    // L is proved up to orientation, so the MIN over its variants.
    if (std::optional<CompositeName> cp = compositeParts(name)) {
        const std::string knot =
            cp->knot[0] == 'm' ? cp->knot.substr(1) : cp->knot;
        UpperContribution k = upperOf(knot, bounds, names);
        const std::vector<std::string> variants = names.candidates(cp->link);
        const bool registered = !variants.empty() && variants.front() != cp->link;
        if (k.value != NO_UPPER_BOUND && registered) {
            int least = NO_LOWER_BOUND;
            std::vector<std::string> support = k.support;
            bool haveAll = true;
            for (const std::string &v : variants) {
                LowerContribution l = lowerOf(v, bounds, names);
                if (l.value == NO_LOWER_BOUND) {
                    haveAll = false;
                    break;
                }
                least = least == NO_LOWER_BOUND ? l.value : std::min(least, l.value);
                mergeSupport(support, l.support);
            }
            if (haveAll && least != NO_LOWER_BOUND &&
                (best.value == NO_LOWER_BOUND || least - k.value > best.value))
                best = {.value = least - k.value, .support = std::move(support)};
        }
    }
    return best;
}

/** Applies `hi[name] <- min(hi[name], value)`, returning whether it changed. */
bool relaxUpper(std::unordered_map<std::string, Bounds> &bounds,
                const NameTable &names, const std::string &name, int value,
                std::vector<std::string> support, const Witness &w,
                const std::string &via) {
    if (value == NO_UPPER_BOUND)
        return false;
    // No circular dependencies
    if (supportContains(support, name))
        return false;
    const Basis basis =
        support.empty() ? Basis::constructive : Basis::literatureAssisted;
    value = std::max(value, 0); // a genus is never negative
    if (const NameInfo *info = names.find(name);
        info && info->haveLiterature && value > info->litHi)
        return false;
    Bounds &b = bounds[name];
    bool better = value < b.hi;
    bool firmer = value == b.hi && b.basis == Basis::literatureAssisted &&
                  basis == Basis::constructive;
    if (!better && !firmer)
        return false;
    b.hi = value;
    b.basis = basis;
    b.support = std::move(support);
    b.kind = w.kind;
    b.viaName = via;
    b.viaGenus = w.genus;
    b.pairSig = w.pairSig;
    b.tubed = w.tubed;
    return true;
}

/** Applies `lo[name] <- max(lo[name], value)`, returning whether it changed. */
bool relaxLower(std::unordered_map<std::string, Bounds> &bounds,
                const NameTable &names, const std::string &name, int value,
                std::vector<std::string> support) {
    if (value == NO_LOWER_BOUND)
        return false;
    // No circular dependencies
    if (supportContains(support, name))
        return false;
    if (value <= 0) // 0 isn't a new lower bound
        return false;
    if (const NameInfo *info = names.find(name);
        info && info->haveLiterature && value < info->litLo)
        return false;
    Bounds &b = bounds[name];
    if (b.haveLower() && value <= b.lo)
        return false;
    b.lo = value;
    b.lowerSupport = std::move(support);
    return true;
}

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

/** Seeds the bounds that need neither a search nor the literature. */
void seedAxioms(std::unordered_map<std::string, Bounds> &bounds,
                const std::vector<Witness> &witnesses,
                const NameTable &names) {
    auto axiom = [&bounds](const std::string &name) {
        Bounds &b = bounds[name];
        b.hi = 0;
        b.lo = 0;
        b.basis = Basis::constructive;
        b.kind = WitnessKind::direct;
    };
    axiom("Unknot");
    // Every "<n>-component unlink" actually mentioned anywhere: it bounds
    // n disks, which tube into a connected planar surface of genus 0.
    // Likewise every slice composite mentioned anywhere: it bounds a disk.
    for (const Witness &w : witnesses)
        for (const std::string &side : {w.other, w.subject})
            if (side.ends_with("-component unlink") ||
                isElementarySlice(side, names))
                axiom(side);
}

} // namespace

bool isElementarySlice(const std::string &name, const NameTable &names) {
    if (isSliceComposite(name))
        return true;
    std::vector<std::string> parts = knotSummands(name);
    if (parts.empty())
        return false;
    // Tally each summand's concordance class against its inverse.
    std::unordered_map<std::string, int> balance; // reversible: +1 K, -1 mK
    std::unordered_map<std::string, int> copies;  // self-inverse classes
    for (const std::string &p : parts) {
        if (p == "Unknot" || p == "mUnknot")
            continue;
        const bool mirrored = p[0] == 'm';
        const std::string knot = mirrored ? p.substr(1) : p;
        const SymmetryType *type = names.symmetry(knot);
        if (!type)
            return false;
        switch (*type) {
        case SymmetryType::reversible: // -K = mK
            balance[knot] += mirrored ? -1 : 1;
            break;
        case SymmetryType::fullyAmphicheiral: // -K = K = mK
            ++copies[knot];
            break;
        case SymmetryType::negativeAmphicheiral: // -K = K, -(mK) = mK
            ++copies[p];
            break;
        default: // -K needs a reversal marker the name does not carry
            return false;
        }
    }
    for (const auto &[knot, b] : balance)
        if (b != 0)
            return false;
    for (const auto &[cls, n] : copies)
        if (n % 2 != 0)
            return false;
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

std::unordered_map<std::string, Bounds>
propagate(const std::vector<Witness> &witnesses, const NameTable &names) {
    std::unordered_map<std::string, Bounds> bounds;
    seedAxioms(bounds, witnesses, names);

    // Direct witnesses are the constructive base case (surfaces that straight
    // up bound the link)
    for (const Witness &w : witnesses)
        if (w.kind == WitnessKind::direct)
            relaxUpper(bounds, names, w.subject, w.genus, /*support=*/{}, w,
                       "");

    // Relax every cobordism in both directions until nothing moves (this will
    // halt eventually: see propagate()'s doc comment).
    bool changed = true;
    while (changed) {
        changed = false;
        for (const Witness &w : witnesses) {
            if (w.kind != WitnessKind::cobordism)
                continue;
            // A multi-component far side that is not a proven unlink bounds
            // nothing in either direction -- its complement does not
            // determine which link it is (\ref cg_farside). The witness
            // stays in the file; it is only the solver that declines it.
            if (!farSideBearsBound(w))
                continue;

            // Both endpoints, with the component count that belongs to each.
            // `far` is the side supplying the bound; `near` is the side
            // receiving it.
            struct Direction {
                std::string near;
                int nearComponents;
                std::vector<std::string> far;
                int farComponents;
                std::string viaLabel; // how to name the far side
            };
            std::vector<Direction> directions;
            directions.push_back({w.subject, w.subjectComponents,
                                  w.otherCandidates, w.otherComponents,
                                  w.other});
            // The reverse direction bounds the far side FROM the subject, so
            // it needs the far side's identity, not just a bound over a set:
            // only a single-component far side (a knot, by Gordon-Luecke)
            // qualifies. An unlink passes farSideBearsBound() but is an
            // axiom already, and a bound onto it would be meaningless.
            if (w.otherComponents == 1 && w.otherCandidates.size() == 1)
                directions.push_back({w.otherCandidates.front(),
                                      w.otherComponents,
                                      {w.subject},
                                      w.subjectComponents,
                                      w.subject});

            for (const Direction &d : directions) {
                if (d.far.empty())
                    continue;

                // Never bound a name using itself (e.g. the product surface is
                // not interesting).
                if (std::find(d.far.begin(), d.far.end(), d.near) !=
                    d.far.end())
                    continue;

                int worst = 0;
                std::vector<std::string> support;
                bool haveAll = true;
                for (const std::string &c : d.far) {
                    UpperContribution u = upperOf(c, bounds, names);
                    if (u.value == NO_UPPER_BOUND) {
                        haveAll = false;
                        break;
                    }
                    worst = std::max(worst, u.value);
                    mergeSupport(support, u.support);
                }
                const bool unlinkFar =
                    !d.far.empty() &&
                    std::all_of(d.far.begin(), d.far.end(),
                                [](const std::string &c) {
                                    return c.ends_with("-component unlink");
                                });
                const int penalty = unlinkFar ? 0 : d.farComponents - 1;
                if (haveAll && relaxUpper(bounds, names, d.near,
                                          addUpper(worst, w.genus + penalty),
                                          std::move(support), w, d.viaLabel))
                    changed = true;

                // Lower: sound only if it holds for every candidate too,
                // hence the min.
                int best = NO_LOWER_BOUND;
                std::vector<std::string> lowerSupport;
                bool haveAllLower = true;
                for (const std::string &c : d.far) {
                    LowerContribution l = lowerOf(c, bounds, names);
                    if (l.value == NO_LOWER_BOUND) {
                        haveAllLower = false;
                        break;
                    }
                    best = best == NO_LOWER_BOUND ? l.value
                                                  : std::min(best, l.value);
                    mergeSupport(lowerSupport, l.support);
                }
                if (haveAllLower && best != NO_LOWER_BOUND &&
                    relaxLower(bounds, names, d.near,
                               best - w.genus - d.nearComponents + 1,
                               std::move(lowerSupport)))
                    changed = true;
            }
        }
    }

    return bounds;
}

/* Reporting output */

Verdict judge(const std::string &name, const Bounds &bounds,
              const NameTable &names) {
    Verdict v;
    v.bounds = bounds;
    if (const NameInfo *info = names.find(name); info && info->haveLiterature) {
        v.haveLiterature = true;
        v.litLo = info->litLo;
        v.litHi = info->litHi;
    }

    if (v.haveLiterature) {
        if (bounds.haveUpper() && bounds.hi < v.litLo) {
            v.status = Status::contradiction;
            std::ostringstream msg;
            msg << name << ": constructed a surface of genus " << bounds.hi
                << ", BELOW the literature lower bound " << v.litLo << ".";
            v.reason = msg.str();
            return v;
        }
        if (bounds.haveLower() && bounds.lo > v.litHi) {
            v.status = Status::contradiction;
            std::ostringstream msg;
            msg << name << ": derived a lower bound of " << bounds.lo
                << ", ABOVE the literature upper bound " << v.litHi << ".";
            v.reason = msg.str();
            return v;
        }
        if (bounds.haveUpper() && bounds.hi == v.litLo) {
            // Only an independently-constructed bound counts as a
            // verification
            v.status = bounds.basis == Basis::constructive
                           ? Status::verified
                           : Status::verifiedAssisted;
            v.value = v.litLo;
            return v;
        }
        if (bounds.haveUpper() && bounds.hi < v.litHi) {
            v.status = Status::improved;
            v.value = bounds.hi;
            return v;
        }
    }

    if (bounds.haveUpper() && bounds.haveLower() && bounds.hi == bounds.lo) {
        v.status = Status::pinned;
        v.value = bounds.hi;
        return v;
    }
    if (bounds.haveUpper() || bounds.haveLower()) {
        v.status = Status::bounded;
        v.value = bounds.haveUpper() ? bounds.hi : bounds.lo;
        return v;
    }
    v.status = Status::unresolved;
    return v;
}

std::string
buildDependsOn(const std::string &via,
               const std::unordered_map<std::string, Bounds> &bounds) {
    std::vector<std::string> chain;
    std::unordered_set<std::string> seen;
    std::string cur = via;
    while (!cur.empty() && seen.insert(cur).second) {
        chain.push_back(cur);
        if (cur == "Unknot" || cur.ends_with("-component unlink"))
            break;
        auto it = bounds.find(cur);
        if (it == bounds.end() || it->second.viaName.empty())
            break;
        cur = it->second.viaName;
    }
    std::ostringstream out;
    for (size_t i = 0; i < chain.size(); ++i) {
        if (i)
            out << ';';
        out << chain[i];
    }
    return out.str();
}

/* Boundary classification */

BoundarySplit
splitBoundary(const std::vector<BoundaryComponentNames> &boundaryComponents,
              size_t searchSideBC, const std::string &rowOwnName) {
    BoundarySplit result;
    for (const auto &info : boundaryComponents) {
        std::optional<std::string> name =
            info.curveNames.size() == 1
                ? std::optional<std::string>(info.curveNames.front())
                : info.linkName;

        if (info.component == searchSideBC && name && *name == rowOwnName) {
            result.searchCurveCount = info.curveNames.size();
            continue;
        }

        if (!name)
            continue; // describeBoundary_() never leaves a multi-curve
                      // component without a linkName

        result.otherSides.push_back(
            {*name, static_cast<int>(info.curveNames.size())});
    }
    return result;
}

/* Orientation matching */

namespace {
// Maps ambient vertex `v` to its corresponding vertex in `dest`, via `iso`
// (which must map `v`'s own triangulation to `dest`).
size_t mapVertexIndex(const regina::Vertex<3> *v,
                      const regina::Triangulation<3> &dest,
                      const regina::Isomorphism<3> &iso) {
    auto emb = v->front();
    size_t destTet = iso.simpImage(emb.tetrahedron()->index());
    int destLocal = iso.facetPerm(emb.tetrahedron()->index())[emb.vertex()];
    return dest.tetrahedron(destTet)->vertex(destLocal)->index();
}
} // namespace

RowOrientation
buildRowOrientation(const std::vector<const regina::Edge<3> *> &rowEdges,
                    const std::vector<bool> &rowReversed,
                    const regina::Triangulation<3> &searchSideTri) {
    if (rowEdges.empty())
        throw regina::InvalidArgument(
            "buildRowOrientation(): rowEdges must not be empty");

    const regina::Triangulation<3> &rowTri = rowEdges.front()->triangulation();
    std::optional<regina::Isomorphism<3>> iso =
        rowTri.isIsomorphicTo(searchSideTri);
    if (!iso)
        throw regina::InvalidArgument(
            "buildRowOrientation(): the row's own triangulation is not "
            "isomorphic to searchSideTri");

    RowOrientation result;
    for (size_t i = 0; i < rowEdges.size(); ++i) {
        const regina::Vertex<3> *tail =
            rowReversed[i] ? rowEdges[i]->vertex(1) : rowEdges[i]->vertex(0);
        const regina::Vertex<3> *head =
            rowReversed[i] ? rowEdges[i]->vertex(0) : rowEdges[i]->vertex(1);
        result.headOf[mapVertexIndex(tail, searchSideTri, *iso)] =
            mapVertexIndex(head, searchSideTri, *iso);
    }
    return result;
}

bool matchesRowOrientation(const RowOrientation &row,
                           const std::vector<OrientedCurve> &curves) {
    std::optional<bool> overallMatch;
    for (const OrientedCurve &curve : curves) {
        if (curve.empty())
            continue;

        std::optional<bool> curveMatch;
        for (const OrientedEdge &oe : curve) {
            const regina::Vertex<3> *tail =
                oe.reversed ? oe.edge->vertex(1) : oe.edge->vertex(0);
            const regina::Vertex<3> *head =
                oe.reversed ? oe.edge->vertex(0) : oe.edge->vertex(1);

            bool edgeMatches;
            auto it = row.headOf.find(tail->index());
            if (it != row.headOf.end() && it->second == head->index()) {
                edgeMatches = true;
            } else {
                auto it2 = row.headOf.find(head->index());
                if (it2 != row.headOf.end() && it2->second == tail->index()) {
                    edgeMatches = false;
                } else {
                    return false; // not one of the row's own tagged edges
                }
            }

            if (!curveMatch)
                curveMatch = edgeMatches;
            else if (*curveMatch != edgeMatches)
                return false; // shouldn't happen
        }

        if (!curveMatch)
            continue;
        if (!overallMatch)
            overallMatch = curveMatch;
        else if (*overallMatch != *curveMatch)
            return false; // mixed pattern across components, reject
    }

    return overallMatch.has_value();
}

} // namespace cobordismgraph

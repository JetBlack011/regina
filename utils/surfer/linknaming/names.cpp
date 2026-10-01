//
//  names.cpp
//

#include "linknaming/names.h"

#include <algorithm>
#include <cctype>

namespace cobordismgraph {

/* Names of knots/links that appear in the cobordism graph */

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
    return stripOrientationTag(name);
}

std::string stripOrientationTag(const std::string &name) {
    return name.substr(0, name.find('{'));
}

std::string stripCensusSuffix(const std::string &name) {
    return name.substr(0, name.find(" : "));
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
    size_t k0 = 0;
    if (k0 < knot.size() && knot[k0] == 'm') ++k0;
    if (k0 < knot.size() && knot[k0] == 'r') ++k0;
    const size_t us = knot.find('_', k0);
    if (us == std::string::npos)
        return std::nullopt;
    const size_t de =
        (us > k0 && (knot[us - 1] == 'a' || knot[us - 1] == 'n')) ? us - 1 : us;
    if (!allDigits(knot, k0, de) || !allDigits(knot, us + 1, knot.size()))
        return std::nullopt;
    if (comp != "?" && !allDigits(comp, 0, comp.size()))
        return std::nullopt;
    size_t i = 1;
    if (link.empty() || link[0] != 'L')
        return std::nullopt;
    while (i < link.size() && std::isdigit(static_cast<unsigned char>(link[i])))
        ++i;
    // An exact name writes L with its orientation tag, "L7n1{1}".
    size_t end = link.size();
    const size_t brace = link.find('{');
    if (brace != std::string::npos) {
        if (link.back() != '}' || brace + 2 > link.size())
            return std::nullopt;
        for (size_t k = brace + 1; k + 1 < link.size(); ++k)
            if (!std::isdigit(static_cast<unsigned char>(link[k])) && link[k] != ';')
                return std::nullopt;
        end = brace;
    }
    if (i == 1 || i >= end || (link[i] != 'a' && link[i] != 'n') ||
        !allDigits(link, i + 1, end))
        return std::nullopt;
    return CompositeName{knot, comp == "?" ? -1 : std::stoi(comp), link};
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
        // m? r? <digits> [a|n]? _ <digits>: "5_2", "m5_2", "mr8_17", "11a_367".
        size_t k0 = 0;
        if (k0 < p.size() && p[k0] == 'm') ++k0;
        if (k0 < p.size() && p[k0] == 'r') ++k0;
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

std::string stripKnotMarks(const std::string &name) {
    size_t k = 0;
    if (k < name.size() && name[k] == 'm') ++k;
    if (k < name.size() && name[k] == 'r') ++k;
    if (k == 0 || k == name.size() || !std::isdigit(static_cast<unsigned char>(name[k])))
        return name;
    const std::string rest = name.substr(k);
    // Only a knot name: digits, optional a/n, '_', digits.
    size_t i = 0;
    while (i < rest.size() && std::isdigit(static_cast<unsigned char>(rest[i]))) ++i;
    if (i < rest.size() && (rest[i] == 'a' || rest[i] == 'n')) ++i;
    if (i >= rest.size() || rest[i] != '_') return name;
    size_t j = ++i;
    while (j < rest.size() && std::isdigit(static_cast<unsigned char>(rest[j]))) ++j;
    return (j > i && j == rest.size()) ? rest : name;
}

namespace {

// A table knot's name, unmarked: digits, an optional a/n, '_', digits.
bool isKnotNameShape(const std::string &s) {
    size_t i = 0;
    while (i < s.size() && std::isdigit(static_cast<unsigned char>(s[i]))) ++i;
    if (i == 0) return false;
    if (i < s.size() && (s[i] == 'a' || s[i] == 'n')) ++i;
    if (i >= s.size() || s[i] != '_') return false;
    size_t j = ++i;
    while (j < s.size() && std::isdigit(static_cast<unsigned char>(s[j]))) ++j;
    return j > i && j == s.size();
}

std::vector<std::string> splitOn(const std::string &s, const std::string &sep) {
    std::vector<std::string> out;
    size_t pos = 0;
    for (;;) {
        size_t at = s.find(sep, pos);
        out.push_back(s.substr(pos, at == std::string::npos ? std::string::npos : at - pos));
        if (at == std::string::npos) return out;
        pos = at + sep.size();
    }
}

} // namespace

// A tagged table link, "L7n1{1}": its component count is in the tag.
bool isTaggedLinkName(const std::string &s) {
    if (s.size() < 5 || s[0] != 'L' || s.back() != '}') return false;
    const size_t brace = s.find('{');
    if (brace == std::string::npos) return false;
    size_t i = 1;
    while (i < brace && std::isdigit(static_cast<unsigned char>(s[i]))) ++i;
    if (i == 1 || i >= brace || (s[i] != 'a' && s[i] != 'n')) return false;
    for (size_t k = i + 1; k < brace; ++k)
        if (!std::isdigit(static_cast<unsigned char>(s[k]))) return false;
    for (size_t k = brace + 1; k + 1 < s.size(); ++k)
        if (!std::isdigit(static_cast<unsigned char>(s[k])) && s[k] != ';') return false;
    return brace + 1 < s.size() - 1;
}

std::optional<std::vector<std::pair<std::string, int>>> sumPieces(const std::string &name) {
    std::vector<std::string> pieces;
    if (name.size() > 3 && name.starts_with("#{") && name.back() == '}') {
        for (const std::string &site : splitOn(name.substr(2, name.size() - 3), " ; "))
            for (std::string tok : splitOn(site, " # ")) {
                // Drop the component index, "[?]" or "[2]".
                if (!tok.empty() && tok.back() == ']') {
                    const size_t open = tok.rfind('[');
                    if (open == std::string::npos) return std::nullopt;
                    tok = tok.substr(0, open);
                }
                pieces.push_back(tok);
            }
    } else if (name.find(" #_") != std::string::npos &&
               name.find(kSplitSeparator) == std::string::npos) {
        std::vector<std::string> terms = splitOn(name, " #_");
        pieces.push_back(terms[0]);
        for (size_t t = 1; t < terms.size(); ++t) {
            const size_t space = terms[t].find(' ');
            if (space == std::string::npos) return std::nullopt;
            pieces.push_back(terms[t].substr(space + 1)); // after "? " or "<c> "
        }
    } else {
        return std::nullopt;
    }
    std::vector<std::pair<std::string, int>> out;
    for (const std::string &p : pieces) {
        if (p.find('|') != std::string::npos) return std::nullopt;
        if (p == "Unknot") {
            out.emplace_back(p, 1);
        } else if (isKnotNameShape(stripKnotMarks(p))) {
            out.emplace_back(stripKnotMarks(p), 1);
        } else if (!knotSummands(p).empty()) {
            out.emplace_back(p, 1);
        } else if (isTaggedLinkName(p)) {
            out.emplace_back(p, componentsFromName(p));
        } else {
            return std::nullopt;
        }
    }
    return out;
}

std::vector<std::string> exactCandidates(const std::string &name) {
    if (name.find(kSplitSeparator) != std::string::npos ||
        name.find(" #_") != std::string::npos || name.find("#{") != std::string::npos)
        return {name};
    return factorAlternatives(name);
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

} // namespace cobordismgraph

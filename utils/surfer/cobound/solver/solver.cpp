//
//  solver.cpp
//

#include "cobound/solver/solver.h"

#include <algorithm>
#include <cctype>
#include <optional>
#include <sstream>
#include <unordered_set>

#include "linknaming/complement/unlinknaming.h"
#include "linknaming/names.h"

namespace solver {

bool farSideBearsBound(const cobordisms::Cobordism &w) {
    // The observed count, never the name: componentsFromName() is an
    // inference from a string, and an alias could make a two-curve far side
    // read like a knot.
    // A far side proved per witness -- meridians carried, or an all-knot
    // split -- is the one case where a multi-component name is a complete
    // candidate set; see Witness::farSideProved.
    return w.otherComponents == 1 || w.farSideProved ||
           complement::isUnlinkName(w.other);
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

// A prime knot's name, possibly marked: "3_1", "m3_1", "mr8_17", "11n_34".
bool isPrimeKnotName(const std::string &alt) {
    size_t k0 = 0;
    if (k0 < alt.size() && alt[k0] == 'm') ++k0;
    if (k0 < alt.size() && alt[k0] == 'r') ++k0;
    const std::string s = alt.substr(k0);
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
    return alt == "Unknot" || !linknaming::knotSummands(alt).empty() ||
           isPrimeKnotName(alt);
}

// A prime knot's mirror image has the same g_4, and the name table holds only
// the unmirrored name, so a bound for "m3_1" is 3_1's. A COMPOSITE is left
// alone: stripping the leading m of "m3_1#3_1" (the square knot, slice)
// would give "3_1#3_1" (the granny, g_4 = 2), a different knot. Its summands
// are unmirrored where the composite rule reads them.
std::string unmirrored(const std::string &alt) {
    return isPrimeKnotName(alt) ? linknaming::stripKnotMarks(alt) : alt;
}

// A split factor's component count, when its name states it (every
// alternative agreeing): 1 for a knot, the tag's for a table link, the link's
// for "K #_? L". nullopt otherwise -- "#{...}" does not say how many of its
// written occurrences are one piece.
std::optional<int> factorComponents(const std::string &factor) {
    std::optional<int> count;
    for (const std::string &alt : linknaming::factorAlternatives(factor)) {
        int n;
        if (isKnotFactor(alt)) {
            n = 1;
        } else if (linknaming::isTaggedLinkName(alt)) {
            n = linknaming::componentsFromName(alt);
        } else if (auto pieces = linknaming::sumPieces(alt);
                   pieces && alt.find("#{") == std::string::npos) {
            n = pieces->back().second;
        } else {
            return std::nullopt;
        }
        if (count && *count != n) return std::nullopt;
        count = n;
    }
    return count;
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
    if (std::vector<std::string> factors = linknaming::splitFactors(name);
        !factors.empty()) {
        int total = 0;
        std::vector<std::string> support;
        bool haveAll = true;
        for (const std::string &f : factors) {
            // Worst case over a factor's alternatives, the same convention
            // candidates() uses: the bound has to hold whichever it is.
            int worstAlt = NO_UPPER_BOUND;
            std::vector<std::string> altSupport;
            for (const std::string &alt : linknaming::factorAlternatives(f)) {
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
    if (std::optional<linknaming::CompositeName> cp = linknaming::compositeParts(name)) {
        const std::string knot = linknaming::stripKnotMarks(cp->knot);
        UpperContribution k = upperOf(knot, bounds, names);
        // A TAGGED L comes from an exact name: it is its own only variant.
        const bool tagged = cp->link.find('{') != std::string::npos;
        const std::vector<std::string> variants =
            tagged ? std::vector<std::string>{cp->link} : names.candidates(cp->link);
        const bool registered =
            tagged || (!variants.empty() && variants.front() != cp->link);
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
    // A sum along components (--sum-rules): g_4(A # B) <= g_4(A) + g_4(B),
    // boundary-connect-summing the pieces' minimal connected surfaces along
    // the summed components -- constructive, whichever the components are.
    if (names.sumRules())
        if (auto pieces = linknaming::sumPieces(name)) {
            int total = 0;
            std::vector<std::string> support;
            bool haveAll = true;
            for (const auto &[piece, n] : *pieces) {
                UpperContribution u = upperOf(piece, bounds, names);
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
    // A composite KNOT A#B#...: g_4 is subadditive under connected sum, so
    // g_4 <= sum of the summands' g_4 (mirrors look up as their knot, g_4
    // being mirror-invariant).
    if (std::vector<std::string> parts = linknaming::knotSummands(name); !parts.empty()) {
        int total = 0;
        std::vector<std::string> support;
        bool haveAll = true;
        for (const std::string &p : parts) {
            if (p == "Unknot" || p == "mUnknot")
                continue;
            UpperContribution u = upperOf(linknaming::stripKnotMarks(p), bounds, names);
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
    if (std::vector<std::string> factors = linknaming::splitFactors(name);
        !factors.empty()) {
        std::vector<std::vector<std::string>> knots;   // nontrivial factors
        bool allKnots = true;
        for (const std::string &f : factors) {
            std::vector<std::string> alts;
            bool trivial = true;
            for (const std::string &alt : linknaming::factorAlternatives(f)) {
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
    // --sum-rules, a split with ANY factors: cap every other factor off in a
    // collar with its minimal connected surface. Gluing a connected surface
    // on along k circles adds k - 1 to the genus, so
    //     g_4(F_i) <= g_4(u F) + sum_{j != i} (g_4(F_j) + n(F_j) - 1).
    // For knot factors that is the rule above, less its -(f - 1).
    if (names.sumRules())
        if (std::vector<std::string> factors = linknaming::splitFactors(name); !factors.empty())
            for (size_t i = 0; i < factors.size(); ++i) {
                int value = NO_LOWER_BOUND;
                std::vector<std::string> support;
                bool ok = true;
                for (const std::string &alt : linknaming::factorAlternatives(factors[i])) {
                    LowerContribution l = lowerOf(unmirrored(alt), bounds, names);
                    if (l.value == NO_LOWER_BOUND) {
                        ok = false;
                        break;
                    }
                    value = value == NO_LOWER_BOUND ? l.value : std::min(value, l.value);
                    mergeSupport(support, l.support);
                }
                for (size_t j = 0; ok && j < factors.size(); ++j) {
                    if (j == i)
                        continue;
                    const std::optional<int> nj = factorComponents(factors[j]);
                    if (!nj) {
                        ok = false;
                        break;
                    }
                    int worst = NO_UPPER_BOUND;
                    for (const std::string &alt : linknaming::factorAlternatives(factors[j])) {
                        UpperContribution u = upperOf(unmirrored(alt), bounds, names);
                        if (u.value == NO_UPPER_BOUND) {
                            ok = false;
                            break;
                        }
                        worst = worst == NO_UPPER_BOUND ? u.value : std::max(worst, u.value);
                        mergeSupport(support, u.support);
                    }
                    if (ok)
                        value -= worst + *nj - 1;
                }
                if (ok && (best.value == NO_LOWER_BOUND || value > best.value))
                    best = {.value = value, .support = std::move(support)};
            }
    // --sum-rules, a sum along components: undoing a piece B (summing -B into
    // the same component; B # -B bounds (B^3, T_B) x I, a disc and n(B) - 1
    // annuli, each annulus a handle once glued on) costs g_4(B) + n(B) - 1,
    //     g_4(P_i) <= g_4(sum) + sum_{j != i} (g_4(P_j) + n(P_j) - 1).
    if (names.sumRules())
        if (auto pieces = linknaming::sumPieces(name))
            for (size_t i = 0; i < pieces->size(); ++i) {
                LowerContribution l = lowerOf((*pieces)[i].first, bounds, names);
                if (l.value == NO_LOWER_BOUND)
                    continue;
                int value = l.value;
                std::vector<std::string> support = l.support;
                bool haveAll = true;
                for (size_t j = 0; j < pieces->size(); ++j) {
                    if (j == i)
                        continue;
                    UpperContribution u = upperOf((*pieces)[j].first, bounds, names);
                    if (u.value == NO_UPPER_BOUND) {
                        haveAll = false;
                        break;
                    }
                    value -= u.value + (*pieces)[j].second - 1;
                    mergeSupport(support, u.support);
                }
                if (haveAll && (best.value == NO_LOWER_BOUND || value > best.value))
                    best = {.value = value, .support = std::move(support)};
            }
    // A composite KNOT: K_i is concordant to (A # -rest), so
    //     g_4(K_i) <= g_4(A) + sum_{j != i} g_4(K_j),
    // i.e. g_4(A) >= g_4(K_i) - sum_{j != i} g_4(K_j), for every i. This is
    // the only lower bound a sum admits in terms of its summands: K # -K is
    // slice, so nothing like g_4(K) + g_4(J) - c holds.
    if (std::vector<std::string> parts = linknaming::knotSummands(name); !parts.empty()) {
        std::vector<std::string> knots;
        for (const std::string &p : parts)
            if (p != "Unknot" && p != "mUnknot")
                knots.push_back(linknaming::stripKnotMarks(p));
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
    if (std::optional<linknaming::CompositeName> cp = linknaming::compositeParts(name)) {
        const std::string knot = linknaming::stripKnotMarks(cp->knot);
        UpperContribution k = upperOf(knot, bounds, names);
        const bool tagged = cp->link.find('{') != std::string::npos;
        const std::vector<std::string> variants =
            tagged ? std::vector<std::string>{cp->link} : names.candidates(cp->link);
        const bool registered =
            tagged || (!variants.empty() && variants.front() != cp->link);
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
                std::vector<std::string> support, const cobordisms::Cobordism &w,
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
    b.pairSigOffset = w.fileOffset;
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

/** Seeds the bounds that need neither a search nor the literature. */
void seedAxioms(std::unordered_map<std::string, Bounds> &bounds,
                const std::vector<cobordisms::Cobordism> &cobordisms,
                const NameTable &names) {
    auto axiom = [&bounds](const std::string &name) {
        Bounds &b = bounds[name];
        b.hi = 0;
        b.lo = 0;
        b.basis = Basis::constructive;
        b.kind = cobordisms::CobordismKind::direct;
    };
    axiom("Unknot");
    // Every "<n>-component unlink" actually mentioned anywhere: it bounds
    // n disks, which tube into a connected planar surface of genus 0.
    // Likewise every slice composite mentioned anywhere: it bounds a disk.
    for (const cobordisms::Cobordism &w : cobordisms)
        for (const std::string &side : {w.other, w.subject})
            if (complement::isMultiComponentUnlinkName(side) ||
                linknaming::isElementarySlice(side, names.symmetries()))
                axiom(side);
}

} // namespace

std::unordered_map<std::string, Bounds>
propagate(const std::vector<cobordisms::Cobordism> &cobordisms, const NameTable &names,
          const std::vector<ExternalProof> &external) {
    std::unordered_map<std::string, Bounds> bounds;
    seedAxioms(bounds, cobordisms, names);

    // Direct witnesses are the constructive base case (surfaces that straight
    // up bound the link)
    for (const cobordisms::Cobordism &w : cobordisms)
        if (w.kind == cobordisms::CobordismKind::direct)
            relaxUpper(bounds, names, w.subject, w.genus, /*support=*/{}, w,
                       "");

    // Certified proofs from outside (--cascade-proofs): a base case too, as
    // constructive as their leaves. Recorded as direct, via the proof's
    // source, so a report names where the bound came from.
    for (const ExternalProof &p : external) {
        cobordisms::Cobordism w;
        w.kind = cobordisms::CobordismKind::direct;
        w.subject = p.name;
        w.genus = p.genus;
        relaxUpper(bounds, names, p.name, p.genus, p.support, w, p.source);
    }

    // Relax every cobordism in both directions until nothing moves (this will
    // halt eventually: see propagate()'s doc comment).
    bool changed = true;
    while (changed) {
        changed = false;
        for (const cobordisms::Cobordism &w : cobordisms) {
            if (w.kind != cobordisms::CobordismKind::cobordism)
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
            // A far side named EXACTLY is an identity whatever its component
            // count, so it may receive a bound too (never an unlink, which is
            // an axiom).
            if ((w.otherComponents == 1 || w.farSideExact) &&
                w.otherCandidates.size() == 1 &&
                !complement::isMultiComponentUnlinkName(w.otherCandidates.front()))
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
                                complement::isMultiComponentUnlinkName);
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

    if (bounds.haveUpper() && bounds.haveLower() && bounds.lo > bounds.hi) {
        v.status = Status::contradiction;
        std::ostringstream msg;
        msg << name << ": derived a lower bound of " << bounds.lo
            << ", ABOVE the derived upper bound " << bounds.hi << ".";
        v.reason = msg.str();
        return v;
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
        if (complement::isUnlinkName(cur))
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

} // namespace solver

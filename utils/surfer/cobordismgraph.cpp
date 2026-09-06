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

std::string baseName(const std::string &name) {
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

int componentsFromName(const std::string &name) {
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

bool NameTable::hasVariants(const std::string &name) const {
    auto it = byBase_.find(baseName(name));
    return it != byBase_.end() && !it->second.empty();
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
                const std::vector<Witness> &witnesses) {
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
            if (side.ends_with("-component unlink") || isSliceComposite(side))
                axiom(side);
}

} // namespace

std::unordered_map<std::string, Bounds>
propagate(const std::vector<Witness> &witnesses, const NameTable &names) {
    std::unordered_map<std::string, Bounds> bounds;
    seedAxioms(bounds, witnesses);

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
            // The reverse direction is only available when we know which link
            // the far side actually is.
            const bool farSideOrientationKnown =
                w.otherComponents == 1 || names.hasVariants(w.other);

            if (w.otherCandidates.size() == 1 && farSideOrientationKnown)
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

//
//  cobordismgraph.cpp
//
//  Created by John Teague on 08/02/2026.
//

#include "cobordismgraph.h"

#include <iostream>
#include <sstream>
#include <unordered_set>

namespace cobordismgraph {

void addEdge(Graph &graph, const std::string &a, const std::string &b,
            int genus, const std::string &pairSig) {
    graph[a].push_back({b, genus, pairSig});
    graph[b].push_back({a, genus, pairSig});
}

bool hasEdgeWithGenus(const Graph &graph, const std::string &a,
                      const std::string &b, int genus) {
    auto it = graph.find(a);
    if (it == graph.end())
        return false;
    for (const auto &e : it->second)
        if (e.other == b && e.genus == genus)
            return true;
    return false;
}

std::optional<ResolveInfo>
resolvable(const std::string &name, int target, const Graph &graph,
          const std::unordered_map<std::string, int> &knownGenus) {
    auto it = graph.find(name);
    if (it == graph.end())
        return std::nullopt;
    for (const auto &e : it->second) {
        auto hIt = knownGenus.find(e.other);
        if (hIt == knownGenus.end())
            continue;
        int h = hIt->second;
        if (target == e.genus + h || target == h - e.genus)
            return ResolveInfo{e.other, e.genus, e.pairSig};
    }
    return std::nullopt;
}

std::string buildDependsOn(
    const std::string &viaKnot,
    const std::unordered_map<std::string, OutputRow> &outputRows) {
    std::vector<std::string> chain;
    std::unordered_set<std::string> seen;
    std::string cur = viaKnot;
    while (!cur.empty() && seen.insert(cur).second) {
        chain.push_back(cur);
        if (cur == "Unknot")
            break;
        auto it = outputRows.find(cur);
        if (it == outputRows.end() || it->second.viaKnot.empty())
            break;
        cur = it->second.viaKnot;
    }
    std::ostringstream out;
    for (size_t i = 0; i < chain.size(); ++i) {
        if (i)
            out << ';';
        out << chain[i];
    }
    return out.str();
}

std::vector<std::string>
propagateGraph(const std::vector<InputRow> &pending, const Graph &graph,
              std::unordered_map<std::string, int> &knownGenus,
              std::unordered_map<std::string, OutputRow> &outputRows) {
    std::vector<std::string> newlyResolved;

    auto resolve = [&](const std::string &name, int lo, int hi) -> bool {
        auto info = resolvable(name, hi, graph, knownGenus);
        if (!info)
            return false;
        knownGenus[name] = hi;
        OutputRow out;
        out.knot = name;
        out.resolvedGenus = hi;
        out.status = (lo == hi) ? "resolved" : "range";
        out.witnessKind = "propagated";
        out.viaKnot = info->viaKnot;
        out.viaEdgeGenus = info->viaEdgeGenus;
        out.dependsOn = buildDependsOn(info->viaKnot, outputRows);
        out.literatureLo = lo;
        out.literatureHi = hi;
        outputRows[name] = std::move(out);
        newlyResolved.push_back(name);
        return true;
    };

    bool changed = true;
    while (changed) {
        changed = false;

        for (const auto &row : pending) {
            if (knownGenus.contains(row.name))
                continue;
            if (resolve(row.name, row.lo, row.hi))
                changed = true;
        }

        // Snapshot names first: resolve() mutates outputRows, and iterating
        // a map while inserting into it is undefined behavior.
        std::vector<std::string> others;
        others.reserve(outputRows.size());
        for (const auto &[name, row] : outputRows)
            if (!knownGenus.contains(name))
                others.push_back(name);
        for (const auto &name : others) {
            if (knownGenus.contains(name))
                continue; // may have resolved earlier in this same pass
            const OutputRow &existing = outputRows.at(name);
            if (resolve(name, existing.literatureLo, existing.literatureHi))
                changed = true;
        }
    }

    return newlyResolved;
}

namespace {
std::string capture(const std::function<std::string()> &capturePairSig) {
    return capturePairSig ? capturePairSig() : std::string{};
}
} // namespace

WitnessOutcome recordWitness(
    const std::string &subject, int subjectLo, int subjectHi, int target,
    const std::string &viaName, int witnessGenus,
    const std::function<std::string()> &capturePairSig,
    const char *witnessKind, Graph &graph,
    std::unordered_map<std::string, int> &knownGenus,
    std::unordered_map<std::string, OutputRow> &outputRows) {
    WitnessOutcome outcome;

    if (viaName.empty()) {
        // Direct witness: subject's own boundary alone, nothing else.
        if (witnessGenus == target) {
            knownGenus[subject] = target;
            OutputRow out;
            out.knot = subject;
            out.resolvedGenus = target;
            out.status = (subjectLo == subjectHi) ? "resolved" : "range";
            out.witnessKind = witnessKind;
            out.witnessPairSig = capture(capturePairSig);
            out.literatureLo = subjectLo;
            out.literatureHi = subjectHi;
            outputRows[subject] = std::move(out);
            outcome.resolved = true;
        } else if (witnessGenus < subjectLo) {
            std::ostringstream msg;
            msg << subject << ": found a DIRECT genus-" << witnessGenus
                << " surface, BELOW the literature lower bound " << subjectLo
                << ".";
            outcome.fatalBug = true;
            outcome.fatalBugMessage = msg.str();
        } else if (witnessGenus < target) {
            // witnessGenus >= subjectLo is guaranteed here (the branch
            // above would have fired first otherwise), so this is a
            // genuine improvement: a smaller witness than the previously
            // known upper bound (target == subjectHi), still consistent
            // with the proven lower bound -- not resolved, just surfaced
            // so it doesn't go unnoticed.
            std::cout << "\x1b[1;33m[!!] " << subject
                      << ": found a DIRECT genus-" << witnessGenus
                      << " surface, IMPROVING on the previously known upper "
                         "bound "
                      << subjectHi << " (literature range [" << subjectLo
                      << ", " << subjectHi << "])"
                      << (witnessGenus == subjectLo
                              ? " -- matches the literature lower bound "
                                "exactly"
                              : "")
                      << "\x1b[0m\n";
        }
        return outcome;
    }

    // Cobordism witness to viaName.

    // Skip the (expensive -- see capturePairSig's own doc comment) pairSig
    // capture entirely when it wouldn't teach the graph anything new:
    // resolvable()'s check is an exact equality on genus, so two edges to
    // the same neighbor at the SAME genus are fully interchangeable
    // (whichever witness got recorded first is just as good), and with
    // thousands of near-duplicate surfaces per pair, this is the
    // overwhelmingly common case.
    bool haveThisGenusEdge =
        hasEdgeWithGenus(graph, subject, viaName, witnessGenus);

    auto hIt = knownGenus.find(viaName);

    // A genus-h cobordism to a name whose genus is already firmly known
    // caps subject's own genus at hIt->second + h (attach this cobordism
    // to viaName's own established witness). If that cap is BELOW
    // subject's proven literature lower bound, that's a mathematical
    // impossibility, not just "doesn't resolve this time".
    if (hIt != knownGenus.end() &&
            hIt->second + witnessGenus < subjectLo) {
        std::ostringstream msg;
        msg << subject << ": found a COBORDISM genus-" << witnessGenus
            << " surface to " << viaName << " (known genus " << hIt->second
            << "), implying an upper bound of "
            << (hIt->second + witnessGenus) << " for " << subject
            << ", BELOW the literature lower bound " << subjectLo << ".";
        outcome.fatalBug = true;
        outcome.fatalBugMessage = msg.str();
        return outcome;
    }

    bool resolvesNow =
        hIt != knownGenus.end() &&
        (target == witnessGenus + hIt->second ||
         target == hIt->second - witnessGenus);

    if (haveThisGenusEdge && !resolvesNow)
        return outcome;

    // As the direct case's analogous branch: attaching this cobordism to
    // viaName's own already-established witness constructs an explicit
    // genus-(hIt->second + witnessGenus) surface for subject. If that's a
    // genuine improvement over the previously known upper bound (target
    // == subjectHi), and it isn't what resolvable() is already about to
    // accept below, surface it -- otherwise it'd go unnoticed.
    if (hIt != knownGenus.end() && !resolvesNow) {
        int implied = hIt->second + witnessGenus;
        if (implied >= subjectLo && implied < target)
            std::cout
                << "\x1b[1;33m[!!] " << subject << ": found a COBORDISM genus-"
                << witnessGenus << " surface to " << viaName
                << " (known genus " << hIt->second
                << "), implying an upper bound of " << implied << " for "
                << subject << ", IMPROVING on the previously known upper "
                             "bound "
                << subjectHi << " (literature range [" << subjectLo << ", "
                << subjectHi << "])"
                << (implied == subjectLo
                        ? " -- matches the literature lower bound exactly"
                        : "")
                << "\x1b[0m\n";
    }

    std::string pairSig = capture(capturePairSig);
    if (!haveThisGenusEdge)
        addEdge(graph, subject, viaName, witnessGenus, pairSig);

    if (resolvesNow) {
        knownGenus[subject] = target;
        OutputRow out;
        out.knot = subject;
        out.resolvedGenus = target;
        out.status = (subjectLo == subjectHi) ? "resolved" : "range";
        out.witnessKind = witnessKind;
        out.witnessPairSig = pairSig;
        out.viaKnot = viaName;
        out.viaEdgeGenus = witnessGenus;
        out.dependsOn = buildDependsOn(viaName, outputRows);
        out.literatureLo = subjectLo;
        out.literatureHi = subjectHi;
        outputRows[subject] = std::move(out);
        outcome.resolved = true;
    }
    return outcome;
}

BoundarySplit splitBoundary(
    const std::vector<BoundaryComponentNames> &boundaryComponents,
    size_t knotSideBC) {
    BoundarySplit result;
    for (const auto &info : boundaryComponents) {
        if (info.component == knotSideBC) {
            result.mineCurveCount = info.curveNames.size();
        } else if (info.curveNames.size() == 1) {
            result.otherSides.push_back({info.curveNames.front(), true});
        } else if (info.linkName) {
            result.otherSides.push_back(
                {*info.linkName,
                 identify::isOrientationSafeName(*info.linkName)});
        }
        // else: describeBoundary_() never leaves a multi-curve component
        // without a linkName -- nothing to record if it somehow did.
    }
    return result;
}

} // namespace cobordismgraph

//
//  farsidenaming.cpp
//

#include "farsidenaming.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <fstream>

#include <link/link.h>

#include "cobordismgraph.h"
#include "identifycomplement.h"
#include "knotbuilder/knotbuilder.h"

namespace farside {

namespace {

// Rows of a "Name,PD,..." table: the name, and the PD field (everything
// between the first comma and the last).
template <typename Fn>
void eachRow(const std::string &path, Fn &&fn) {
    if (path.empty()) return;
    std::ifstream in(path);
    if (!in) throw regina::InvalidArgument("cannot open " + path);
    std::string line;
    std::getline(in, line); // header
    while (std::getline(in, line)) {
        size_t a = line.find(','), b = line.rfind(',');
        if (a == std::string::npos || b <= a) continue;
        fn(line.substr(0, a), line.substr(a + 1, b - a - 1));
    }
}

// A table PD code as Regina reads it: the integers in fours, labels as
// written (1..2n). Not knotbuilder::parsePDCode(), which renumbers from 0.
regina::Link linkOf(const std::string &pd) {
    std::vector<std::array<long, 4>> code;
    std::vector<long> n;
    long cur = -1;
    for (char ch : pd + ' ') {
        if (ch >= '0' && ch <= '9') {
            cur = (cur < 0 ? 0 : cur * 10) + (ch - '0');
        } else if (cur >= 0) {
            n.push_back(cur);
            cur = -1;
        }
    }
    if (n.empty() || n.size() % 4 != 0)
        throw regina::InvalidArgument("not a PD code: " + pd);
    for (size_t i = 0; i < n.size(); i += 4)
        code.push_back({n[i], n[i + 1], n[i + 2], n[i + 3]});
    return regina::Link::fromPD(code.begin(), code.end());
}

long long microsSince(std::chrono::steady_clock::time_point t) {
    return std::chrono::duration_cast<std::chrono::microseconds>(
               std::chrono::steady_clock::now() - t)
        .count();
}

} // namespace

SignatureTable SignatureTable::fromTables(const std::string &knotTable,
                                          const std::string &linkTable) {
    // A row that does not parse is an error, not a silent gap: an empty or
    // partial table would quietly send every far side back to the
    // complement route (it did, once, when the PD codes were read with
    // 0-based labels).
    SignatureTable t;
    eachRow(knotTable, [&](const std::string &name, const std::string &pd) {
        t.knots_.try_emplace(linkOf(pd).knotSig(true, true), name);
        t.knotNames_.insert(name);
    });
    eachRow(linkTable, [&](const std::string &name, const std::string &pd) {
        t.links_.try_emplace(linkOf(pd).sig<2>(true, true, true),
                             name.substr(0, name.find('{')));
    });
    if ((!knotTable.empty() && t.knots_.empty()) ||
        (!linkTable.empty() && t.links_.empty()))
        throw regina::InvalidArgument("a table yielded no signatures");
    return t;
}

const std::string *SignatureTable::knot(const std::string &sig) const {
    auto it = knots_.find(sig);
    return it == knots_.end() ? nullptr : &it->second;
}

const std::string *SignatureTable::link(const std::string &sig) const {
    auto it = links_.find(sig);
    return it == links_.end() ? nullptr : &it->second;
}

DiagramNamer::DiagramNamer(const regina::Triangulation<3> &knotT, size_t crossings,
                           const CobordismBuilder<3> &cob, const SignatureTable &table)
    : map_(knotT, cob), drawer_(knotT, crossings), table_(table) {}

std::string DiagramNamer::name(const Link &curves) const {
    ++stats_.calls;
    return identify::perturbedForTesting(nameOnce(curves));
}

std::string DiagramNamer::nameOnce(const Link &curves) const {
    const auto start = std::chrono::steady_clock::now();
    const size_t n = curves.comps_.size();
    std::string key; // what a fallback's answer is remembered under
    try {
        std::vector<knotbuilder::EdgeCycle> cycles;
        cycles.reserve(n);
        for (const Knot &k : curves.comps_) cycles.push_back(map_.carryCycle(k.edges()));
        knotbuilder::Diagram d = drawer_.draw(cycles);

        bool someLinking = false;
        for (size_t i = 0; i < n && !someLinking; ++i)
            for (size_t j = i + 1; j < n && !someLinking; ++j)
                someLinking = d.linkingNumber(i, j) != 0;

        regina::Link simplified = d.link();
        simplified.simplify();
        if (simplified.size() == 0) {
            stats_.microsDiagram += microsSince(start);
            if (n == 1) { ++stats_.unknots; return "Unknot"; }
            ++stats_.unlinks;
            return std::to_string(n) + "-component unlink";
        }
        if (n == 1) {
            std::string sig = simplified.knotSig(true, true);
            if (const std::string *hit = table_.knot(sig)) {
                ++stats_.tableKnots;
                stats_.microsDiagram += microsSince(start);
                return *hit;
            }
            key = "K" + sig;
        } else {
            std::string sig = simplified.sig<2>(true, true, true);
            if (const std::string *hit = table_.link(sig)) {
                ++stats_.tableLinks;
                stats_.microsDiagram += microsSince(start);
                return *hit;
            }
            key = "L" + sig;
            if (someLinking) {
                ++stats_.diagramLinks;
                stats_.microsDiagram += microsSince(start);
                return "diagram:" + sig;
            }
        }
        {
            std::lock_guard<std::mutex> lock(learnedMutex_);
            if (auto it = learned_.find(key); it != learned_.end()) {
                ++(n == 1 ? stats_.learnedKnots : stats_.learnedLinks);
                stats_.microsDiagram += microsSince(start);
                return it->second;
            }
        }
        if (n > 1) {
            // Every linking number is zero. A Jones polynomial other than the
            // unlink's still proves this is not an unlink.
            regina::Laurent<regina::Integer> unlink;
            {
                std::lock_guard<std::mutex> lock(learnedMutex_);
                auto it = unlinkJones_.find(n);
                if (it == unlinkJones_.end())
                    it = unlinkJones_.emplace(n, regina::Link(n).jones()).first;
                unlink = it->second;
            }
            if (simplified.jones() != unlink) {
                ++stats_.jonesLinks;
                stats_.microsDiagram += microsSince(start);
                return "diagram:" + key.substr(1);
            }
        }
    } catch (const knotbuilder::NonPlanar &) {
        ++stats_.nonPlanar;
        key.clear(); // a drawer defect: never name from it, and never learn
    } catch (const knotbuilder::Degenerate &) {
        key.clear(); // fall through to the complement route
    } catch (const regina::InvalidArgument &) {
        key.clear();
    }
    stats_.microsDiagram += microsSince(start);

    // The complement route, as before this namer existed; its answer is
    // remembered against the diagram.
    const auto fb = std::chrono::steady_clock::now();
    ++stats_.fallbacks;
    std::string name = identify::identify(curves);
    stats_.microsFallback += microsSince(fb);
    if (!key.empty()) {
        std::lock_guard<std::mutex> lock(learnedMutex_);
        if (learned_.try_emplace(key, name).second) ++stats_.learned;
    }
    return name;
}

void DiagramNamer::enableExactNames(const exactnaming::ExactTables &tables) {
    exactnaming::NamerLimits fast;
    fast.simplifyTries = 2;
    fast.exhaustiveHeight = 0;
    fast.searchHeight = -1; // no Reidemeister search in the search
    fast.deepHeight = -1;
    exact_ = std::make_unique<exactnaming::ExactNamer>(tables, fast);
}

std::optional<std::string> DiagramNamer::orientedName(
    const std::vector<OrientedCurve> &outgoing,
    const std::map<const regina::Edge<3> *, size_t> &surfaceOf,
    const std::map<size_t, int> &flips) const {
    if (!exact_) return std::nullopt;
    try {
        std::vector<knotbuilder::EdgeCycle> cycles;
        for (const OrientedCurve &curve : outgoing) {
            if (curve.empty()) continue;
            auto comp = surfaceOf.find(curve.front().edge);
            if (comp == surfaceOf.end()) return std::nullopt;
            auto flip = flips.find(comp->second);
            if (flip == flips.end()) return std::nullopt;
            knotbuilder::EdgeCycle cyc = map_.carry(curve);
            if (flip->second < 0) {
                std::reverse(cyc.begin(), cyc.end());
                for (auto &de : cyc) de.reversed = !de.reversed;
            }
            cycles.push_back(std::move(cyc));
        }
        const regina::Link drawn = drawer_.draw(cycles).link();
        const std::string key = drawn.sig<2>(false, false, true);
        {
            std::lock_guard<std::mutex> lock(exactMutex_);
            if (auto it = exactCache_.find(key); it != exactCache_.end()) {
                ++stats_.exactCacheHits;
                return it->second;
            }
        }
        std::string name = exact_->name(drawn).name;
        ++stats_.exactNamed;
        std::lock_guard<std::mutex> lock(exactMutex_);
        return exactCache_.try_emplace(key, std::move(name)).first->second;
    } catch (const std::exception &) {
        ++stats_.exactFailed;
        return std::nullopt;
    }
}

} // namespace farside

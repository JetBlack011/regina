//
//  outgoingnamer.cpp
//

#include "cobound/outgoing/outgoingnamer.h"

#include <algorithm>
#include <chrono>

#include <link/link.h>

#include "cobound/driver/timers.h"
#include "diagramtriangulation/fromdiagram.h"
#include "linknaming/census/censusnaming.h"

namespace outgoing {

// The complement route, as SurfaceSearch called it before a namer was
// required: a lone curve on its own (identify(const EdgeComplement&)),
// several curves together (identify(const Link&)).
ComplementNamer::ComplementNamer()
    : ComplementBoundaryNamer(static_cast<KnotRoute>(&census::identify),
                              static_cast<LinkRoute>(&census::identify)) {}

DiagramNamer::DiagramNamer(const regina::Triangulation<3> &knotT, size_t crossings,
                           const CobordismBuilder<3> &cob, const linknaming::SignatureTable &table)
    : map_(knotT, cob), drawer_(knotT, crossings), namer_(table) {}

std::string DiagramNamer::name(const Link &curves) const {
    return namer_.name(curves, [this, &curves] { return draw(curves); });
}

std::string DiagramNamer::nameLink(size_t bc, const Link &curves) const {
    return handles(bc) ? name(curves) : complement_.nameLink(bc, curves);
}

std::string DiagramNamer::nameCurve(size_t bc, const Knot &curve) const {
    return complement_.nameCurve(bc, curve);
}

linknaming::DrawnCurves DiagramNamer::draw(const Link &curves) const {
    linknaming::DrawnCurves out;
    const size_t n = curves.comps_.size();
    try {
        std::vector<knotbuilder::EdgeCycle> cycles;
        cycles.reserve(n);
        for (const Knot &k : curves.comps_) cycles.push_back(map_.carryCycle(k.edges()));
        knotbuilder::Diagram d = drawer_.draw(cycles);

        bool someLinking = false;
        for (size_t i = 0; i < n && !someLinking; ++i)
            for (size_t j = i + 1; j < n && !someLinking; ++j)
                someLinking = d.linkingNumber(i, j) != 0;

        out.diagram = d.link();
        out.someLinking = someLinking;
        out.outcome = linknaming::DrawnCurves::Outcome::drawn;
    } catch (const knotbuilder::NonPlanar &) {
        out.outcome = linknaming::DrawnCurves::Outcome::nonPlanar;
    } catch (const knotbuilder::Degenerate &) {
        out.outcome = linknaming::DrawnCurves::Outcome::failed;
    } catch (const regina::InvalidArgument &) {
        out.outcome = linknaming::DrawnCurves::Outcome::failed;
    }
    return out;
}

void DiagramNamer::enableExactNames(const linknaming::ExactTables &tables,
                                    std::shared_ptr<linknaming::TableCaches> caches) {
    linknaming::NamerLimits fast;
    fast.simplifyTries = 2;
    fast.exhaustiveHeight = 0;
    fast.searchHeight = -1; // no Reidemeister search in the search
    fast.deepHeight = -1;
    // Nor from the table's side: rewrite() outward from every HOMFLY
    // candidate's diagram, up to 3M diagrams per candidate and height. On a
    // far side that is none of them it runs to the end, a minute or more per
    // name (a 50k-surface hop spent 57,000 thread-seconds there, 2026-09-29).
    // farsidename refines such names offline from the pair signature.
    fast.tableSideHeight = -1;
    exact_ = std::make_unique<linknaming::ExactNamer>(tables, fast, std::move(caches));
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
            knotbuilder::EdgeCycle cyc = map_.carry(outgoingCurve(curve));
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
                ++namer_.stats().exactCacheHits;
                return it->second;
            }
        }
        const auto start = std::chrono::steady_clock::now();
        std::string name = exact_->name(drawn).name;
        const long long micros = timers::microsSince(start);
        namer_.stats().microsExact += micros;
        namer_.stats().noteDuration(micros, "exact", name);
        ++namer_.stats().exactNamed;
        std::lock_guard<std::mutex> lock(exactMutex_);
        return exactCache_.try_emplace(key, std::move(name)).first->second;
    } catch (const std::exception &) {
        ++namer_.stats().exactFailed;
        return std::nullopt;
    }
}

} // namespace outgoing

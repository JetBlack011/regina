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
// required: a lone curve on its own (census::nameComplement(const EdgeComplement&)),
// several curves together (census::nameComplement(const Link&)).
ComplementNamer::ComplementNamer()
    : ComplementBoundaryNamer(static_cast<KnotRoute>(&census::nameComplement),
                              static_cast<LinkRoute>(&census::nameComplement)) {}

OutgoingNamer::OutgoingNamer(const regina::Triangulation<3> &knotT, size_t crossings,
                           const CobordismBuilder<3> &cob, const linknaming::SignatureTable &table)
    : map_(knotT, cob), drawer_(knotT, crossings), diagramNamer_(table) {}

std::string OutgoingNamer::name(const Link &curves) const {
    return diagramNamer_.name(curves, [this, &curves] { return draw(curves); });
}

std::string OutgoingNamer::nameLink(size_t bc, const Link &curves) const {
    return handles(bc) ? name(curves) : complement_.nameLink(bc, curves);
}

std::string OutgoingNamer::nameCurve(size_t bc, const Knot &curve) const {
    return complement_.nameCurve(bc, curve);
}

linknaming::DrawnCurves OutgoingNamer::draw(const Link &curves) const {
    linknaming::DrawnCurves out;
    const size_t n = curves.comps_.size();
    try {
        std::vector<diagramtriangulation::EdgeCycle> cycles;
        cycles.reserve(n);
        for (const Knot &k : curves.comps_) cycles.push_back(map_.carryCycle(k.edges()));
        diagramtriangulation::Diagram d = drawer_.draw(cycles);

        bool someLinking = false;
        for (size_t i = 0; i < n && !someLinking; ++i)
            for (size_t j = i + 1; j < n && !someLinking; ++j)
                someLinking = d.linkingNumber(i, j) != 0;

        out.diagram = d.link();
        out.someLinking = someLinking;
        out.outcome = linknaming::DrawnCurves::Outcome::drawn;
    } catch (const diagramtriangulation::NonPlanar &) {
        out.outcome = linknaming::DrawnCurves::Outcome::nonPlanar;
    } catch (const diagramtriangulation::Degenerate &) {
        out.outcome = linknaming::DrawnCurves::Outcome::failed;
    } catch (const regina::InvalidArgument &) {
        out.outcome = linknaming::DrawnCurves::Outcome::failed;
    }
    return out;
}

void OutgoingNamer::enableOrientedNames(const linknaming::Tables &tables,
                                    std::shared_ptr<linknaming::TableCaches> caches) {
    linknaming::NamerLimits fast;
    fast.simplifyTries = 2;
    fast.exhaustiveHeight = 0;
    fast.searchHeight = -1; // no Reidemeister search in the search
    fast.deepHeight = -1;
    // Nor from the table's side: rewrite() outward from every HOMFLY
    // candidate's diagram, up to 3M diagrams per candidate and height. On a
    // outgoing link that is none of them it runs to the end, a minute or more per
    // name (a 50k-surface search spent 57,000 thread-seconds there, 2026-09-29).
    // `cobound name` refines such names offline from the pair signature.
    fast.tableSideHeight = -1;
    linkNamer_ = std::make_unique<linknaming::LinkNamer>(tables, fast, std::move(caches));
}

std::optional<std::string> OutgoingNamer::orientedName(
    const std::vector<OrientedCurve> &outgoing,
    const std::map<const regina::Edge<3> *, size_t> &surfaceOf,
    const std::map<size_t, int> &flips) const {
    if (!linkNamer_) return std::nullopt;
    try {
        std::vector<diagramtriangulation::EdgeCycle> cycles;
        for (const OrientedCurve &curve : outgoing) {
            if (curve.empty()) continue;
            auto comp = surfaceOf.find(curve.front().edge);
            if (comp == surfaceOf.end()) return std::nullopt;
            auto flip = flips.find(comp->second);
            if (flip == flips.end()) return std::nullopt;
            diagramtriangulation::EdgeCycle cyc = map_.carry(outgoingCurve(curve));
            if (flip->second < 0) {
                std::reverse(cyc.begin(), cyc.end());
                for (auto &de : cyc) de.reversed = !de.reversed;
            }
            cycles.push_back(std::move(cyc));
        }
        const regina::Link drawn = drawer_.draw(cycles).link();
        const std::string key = drawn.sig<2>(false, false, true);
        {
            std::lock_guard<std::mutex> lock(orientedMutex_);
            if (auto it = orientedCache_.find(key); it != orientedCache_.end()) {
                ++diagramNamer_.stats().orientedCacheHits;
                return it->second;
            }
        }
        const auto start = std::chrono::steady_clock::now();
        std::string name = linkNamer_->name(drawn).name;
        const long long micros = timers::microsSince(start);
        diagramNamer_.stats().microsOriented += micros;
        diagramNamer_.stats().noteDuration(micros, "exact", name);
        ++diagramNamer_.stats().orientedNamed;
        std::lock_guard<std::mutex> lock(orientedMutex_);
        return orientedCache_.try_emplace(key, std::move(name)).first->second;
    } catch (const std::exception &) {
        ++diagramNamer_.stats().orientedFailed;
        return std::nullopt;
    }
}

} // namespace outgoing

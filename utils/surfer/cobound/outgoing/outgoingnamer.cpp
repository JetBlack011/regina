//
//  outgoingnamer.cpp
//

#include "cobound/outgoing/outgoingnamer.h"

#include <algorithm>
#include <chrono>
#include <stdexcept>

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

linknaming::NamerLimits OutgoingNamer::limits() {
    linknaming::NamerLimits fast;
    // The attempts each route made before they were merged (phase 7.1): a
    // link's three (one, and two more), a knot's one. More for a knot is a
    // proposal measured for John, not adopted.
    fast.simplifyTries = 2;
    fast.knotSimplifyTries = 0;
    fast.exhaustiveHeight = 0;
    fast.searchHeight = -1; // no Reidemeister search in the search
    fast.deepHeight = -1;
    // Nor from the table's side: rewrite() outward from every HOMFLY
    // candidate's diagram, up to 3M diagrams per candidate and height. On an
    // outgoing link that is none of them it runs to the end, a minute or more
    // per name (a 50k-surface search spent 57,000 thread-seconds there,
    // 2026-09-29). `cobound name` refines such names offline from the pair
    // signature.
    fast.tableSideHeight = -1;
    return fast;
}

OutgoingNamer::OutgoingNamer(const regina::Triangulation<3> &knotT, size_t crossings,
                             const CobordismBuilder<3> &cob, const linknaming::Tables &tables,
                             std::shared_ptr<linknaming::TableCaches> caches)
    : map_(knotT, cob), drawer_(knotT, crossings), namer_(tables, limits(), std::move(caches)) {}

namespace {

const char *routeOf(const linknaming::RoutedName &r, bool oriented) {
    return r.complementCalled ? "complement" : oriented ? "oriented" : "diagram";
}

} // namespace

std::string OutgoingNamer::name(const Link &curves) const {
    const auto start = std::chrono::steady_clock::now();
    const size_t n = curves.comps_.size();
    linknaming::DrawnCurves drawing;
    try {
        std::vector<diagramtriangulation::EdgeCycle> cycles;
        cycles.reserve(n);
        for (const Knot &k : curves.comps_) cycles.push_back(map_.carryCycle(k.edges()));
        drawing = draw(cycles);
    } catch (const regina::InvalidArgument &) {
        drawing.outcome = linknaming::DrawnCurves::Outcome::failed;
    }
    const auto complement = [&curves] { return census::nameComplement(curves); };
    linknaming::RoutedName r;
    try {
        r = namer_.nameDrawing(drawing, n, complement);
    } catch (const std::exception &) {
        // A drawing the namer cannot take (a broken invariant, say): the
        // complement route, as for a drawing that failed.
        drawing.outcome = linknaming::DrawnCurves::Outcome::failed;
        r = namer_.nameDrawing(drawing, n, complement);
    }
    stats_.count(r, /*oriented=*/false, n);
    stats_.noteDuration(timers::microsSince(start), routeOf(r, false), r.name);
    return r.name;
}

std::string OutgoingNamer::nameLink(size_t bc, const Link &curves) const {
    return handles(bc) ? name(curves) : complement_.nameLink(bc, curves);
}

std::string OutgoingNamer::nameCurve(size_t bc, const Knot &curve) const {
    return complement_.nameCurve(bc, curve);
}

linknaming::DrawnCurves
OutgoingNamer::draw(const std::vector<diagramtriangulation::EdgeCycle> &cycles) const {
    linknaming::DrawnCurves out;
    try {
        out.diagram = drawer_.draw(cycles).link();
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

std::string OutgoingNamer::orientedName(
    const std::vector<OrientedCurve> &outgoing,
    const std::map<const regina::Edge<3> *, size_t> &surfaceOf,
    const std::map<size_t, int> &flips) const {
    const auto start = std::chrono::steady_clock::now();
    // The curves' edges, for the complement route.
    std::vector<const regina::Edge<3> *> edges;
    size_t n = 0;
    for (const OrientedCurve &curve : outgoing) {
        if (curve.empty()) continue;
        ++n;
        for (const OrientedEdge &e : curve) edges.push_back(e.edge);
    }
    if (edges.empty())
        throw std::logic_error("OutgoingNamer::orientedName(): no outgoing curve");
    const regina::Triangulation<3> &boundary = edges.front()->triangulation();
    const auto complement = [&] { return census::nameComplement(Link(boundary, edges)); };

    linknaming::DrawnCurves drawing;
    try {
        std::vector<diagramtriangulation::EdgeCycle> cycles;
        for (const OrientedCurve &curve : outgoing) {
            if (curve.empty()) continue;
            auto comp = surfaceOf.find(curve.front().edge);
            if (comp == surfaceOf.end())
                throw std::invalid_argument("an outgoing curve on no surface component");
            auto flip = flips.find(comp->second);
            if (flip == flips.end())
                throw std::invalid_argument("a surface component with no orientation");
            diagramtriangulation::EdgeCycle cyc = map_.carry(outgoingCurve(curve));
            if (flip->second < 0) {
                std::reverse(cyc.begin(), cyc.end());
                for (auto &de : cyc) de.reversed = !de.reversed;
            }
            cycles.push_back(std::move(cyc));
        }
        drawing = draw(cycles);
    } catch (const std::exception &) {
        // Not drawn as this surface orients it: the complement route.
        drawing.outcome = linknaming::DrawnCurves::Outcome::failed;
    }
    linknaming::RoutedName r;
    try {
        r = namer_.nameDrawing(drawing, n, complement);
    } catch (const std::exception &) {
        drawing.outcome = linknaming::DrawnCurves::Outcome::failed;
        r = namer_.nameDrawing(drawing, n, complement);
    }
    stats_.count(r, /*oriented=*/true, n);
    stats_.noteDuration(timers::microsSince(start), routeOf(r, true), r.name);
    return r.name;
}

} // namespace outgoing

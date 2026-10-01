//
//  profile_boundarylinks.cpp
//
//  Where does the boundary-identification drain actually spend its time?
//
//  Timing the identification stages in isolation (build complement, simplify,
//  isoSig, group/unlink test) on real boundary complements accounted for
//  ~1.7 ms per surface. The observed drain rate on a large queue was ~123
//  surfaces/sec across 12 threads, i.e. ~98 ms per surface per thread -- a
//  factor of ~57 unaccounted for. Since those identification stages are also
//  short-circuited by BoundarySignatureCache on repeat boundaries, they cannot
//  be the bulk of it.
//
//  The suspect is KnottedSurface::boundaryLinks() (embeddedsubmanifold.cpp),
//  which for every boundary facet of every surface does
//
//      for (int k = 0; k < ambientBC->countEdges(); ++k)
//          if (ambientBC->edge(k) == ambientFacet) { ...; break; }
//
//  -- a linear scan over the ambient boundary component's entire edge list.
//  For a 2-layer cobordism over a 9-crossing diagram that list has thousands
//  of entries, and the scan runs per facet, per surface.
//
//  This program measures that directly rather than by subtraction: it builds
//  the same cobordism the search uses, reconstructs surfaces from the seed,
//  and times boundaryLinks() against the identification that follows it.
//
//  Usage:  profile_boundarylinks <pd-code> [layers] [samples]
//
#include <chrono>
#include <iostream>
#include <string>
#include <vector>

#include <triangulation/dim3.h>
#include <triangulation/dim4.h>

#include "diagramtriangulation/thickening/thickening.h"
#include "surfer/submanifold/submanifold.h"
#include "linknaming/complement/linkcomplement.h"
#include "diagramtriangulation/fromdiagram.h"
#include "surfer/submanifold/skeleton.h"

namespace {

using Clock = std::chrono::steady_clock;

double msSince(Clock::time_point t0) {
    return std::chrono::duration<double, std::milli>(Clock::now() - t0).count();
}

} // namespace

int main(int argc, char *argv[]) {
    if (argc < 2) {
        std::cerr << "usage: " << argv[0] << " <pd-code> [layers] [samples]\n";
        return 1;
    }
    const std::string pd = argv[1];
    const int layers = argc > 2 ? std::stoi(argv[2]) : 2;
    const int samples = argc > 3 ? std::stoi(argv[3]) : 2000;

    ThickenedLink row; // as the search builds it: the collar through every layer
    buildAmbient(pd, layers, layers, /*useCone=*/false, row);
    const size_t searchSideBC = row.searchSideBC;
    const regina::Triangulation<4> &tri = row.tri;
    const std::vector<int> &seed = row.seedFaces;

    std::cout << "cobordism: " << tri.size() << " pentachora, "
              << tri.countTriangles() << " triangles\n"
              << "seed: " << seed.size() << " faces\n";
    const auto *bc = tri.boundaryComponent(searchSideBC);
    std::cout << "search-side boundary component: " << bc->countEdges()
              << " edges  <-- length of the scan in boundaryLinks()\n\n";

    Skeleton<4, 2> skeleton(tri);
    PetalCache petals;

    // One surface, rebuilt from the seed exactly as the drain does.
    auto t0 = Clock::now();
    KnottedSurface surface(skeleton, petals, seed);
    const double buildMs = msSince(t0);

    // Stage 1: boundaryLinks() -- the suspected hot path.
    t0 = Clock::now();
    for (int i = 0; i < samples; ++i) {
        auto links = surface.boundaryLinks();
        asm volatile("" ::"r"(&links) : "memory"); // keep it
    }
    const double linksMs = msSince(t0) / samples;

    // Stage 2: what identification costs on the SAME boundary, for scale.
    auto links = surface.boundaryLinks();
    size_t curves = 0;
    for (auto &[c, l] : links)
        curves += l.comps_.size();
    t0 = Clock::now();
    const int idSamples = std::min(samples, 200);
    for (int i = 0; i < idSamples; ++i)
        for (auto &[c, l] : links) {
            auto complement = l.buildComplement();
            volatile auto sig = complement.isoSig();
            (void)sig;
        }
    const double idMs = msSince(t0) / idSamples;

    std::cout << "reconstruct surface from seed : " << buildMs << " ms (once)\n"
              << "boundaryLinks()               : " << linksMs << " ms/surface\n"
              << "  (" << curves << " boundary curves found)\n"
              << "buildComplement()+isoSig()    : " << idMs << " ms/surface\n"
              << "\nratio boundaryLinks : identification = "
              << (idMs > 0 ? linksMs / idMs : 0.0) << "\n";
    return 0;
}

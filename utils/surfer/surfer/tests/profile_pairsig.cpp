// Scratch benchmark: how much of pairSig<4,2>() depends only on the AMBIENT
// triangulation (constant for a whole search) versus on the surface?
#include <chrono>
#include <iostream>
#include <string>
#include <vector>
#include <triangulation/dim3.h>
#include <triangulation/dim4.h>
#include "diagramtriangulation/thickening/thickening.h"
#include "surfer/submanifold/submanifold.h"
#include "diagramtriangulation/fromdiagram.h"
#include "surfer/pairsig/pairsig.h"
#include "surfer/submanifold/skeleton.h"

using Clock = std::chrono::steady_clock;
static double msSince(Clock::time_point t0) {
    return std::chrono::duration<double, std::milli>(Clock::now() - t0).count();
}

int main(int argc, char *argv[]) {
    if (argc < 2) { std::cerr << "usage: " << argv[0] << " <pd> [layers] [samples]\n"; return 1; }
    const std::string pd = argv[1];
    const int layers  = argc > 2 ? std::stoi(argv[2]) : 2;
    const int samples = argc > 3 ? std::stoi(argv[3]) : 5;

    ThickenedLink thickened; // as the search builds it: the collar through every layer
    buildAmbient(pd, layers, layers, thickened);
    const regina::Triangulation<4> &tri = thickened.tri;
    const std::vector<int> &seed = thickened.seedFaces;

    std::cout << "cobordism: " << tri.size() << " pentachora, "
              << tri.countTriangles() << " triangles; seed " << seed.size() << " faces\n\n";

    Skeleton<4,2> skeleton(tri);
    PetalCache petals;
    KnottedSurface surface(skeleton, petals, seed);

    // (a) isoSigDetail() on the ambient -- AMBIENT ONLY
    auto t0 = Clock::now();
    std::string sig; 
    for (int i = 0; i < samples; ++i) { auto d = tri.isoSigDetail(); sig = d.first; }
    double isoMs = msSince(t0) / samples;

    // (b) fromSig() rebuild -- AMBIENT ONLY
    t0 = Clock::now();
    for (int i = 0; i < samples; ++i) { auto c = regina::Triangulation<4>::fromSig(sig); (void)c.size(); }
    double fromSigMs = msSince(t0) / samples;

    // (c) findAllIsomorphisms(canon,canon) -- AMBIENT ONLY
    auto canon = regina::Triangulation<4>::fromSig(sig);
    long autos = 0;
    t0 = Clock::now();
    for (int i = 0; i < samples; ++i) {
        autos = 0;
        canon.findAllIsomorphisms(canon, [&](const regina::Isomorphism<4>&){ ++autos; return false; });
    }
    double autoMs = msSince(t0) / samples;

    // (d) the whole pairSig, as the drain calls it
    t0 = Clock::now();
    for (int i = 0; i < samples; ++i) {
        auto s = pairSig<4,2>(skeleton.triangulation(), surface.markedFaces());
        asm volatile("" ::"r"(&s) : "memory");
    }
    double pairMs = msSince(t0) / samples;

    double ambientOnly = isoMs + fromSigMs + autoMs;
    std::cout << "isoSigDetail(ambient)      : " << isoMs     << " ms   [AMBIENT-ONLY]\n"
              << "fromSig(sig)               : " << fromSigMs << " ms   [AMBIENT-ONLY]\n"
              << "findAllIsomorphisms(canon) : " << autoMs    << " ms   [AMBIENT-ONLY]  (" << autos << " automorphisms)\n"
              << "                             ------\n"
              << "ambient-only subtotal      : " << ambientOnly << " ms\n"
              << "full pairSig<4,2>()        : " << pairMs << " ms\n\n"
              << "hoistable fraction         : " << (pairMs>0 ? 100.0*ambientOnly/pairMs : 0.0) << " %\n"
              << "per-witness cost after hoist: ~" << (pairMs - ambientOnly) << " ms\n";
    return 0;
}

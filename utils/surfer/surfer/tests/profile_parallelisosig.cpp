// profile_parallelisosig.cpp
//
// How parallelIsoSigDetail() (../parallelisosig.h) scales: a link's
// campaign-shaped thickening (two layers, collared, no cone), its isoSig at
// each thread count given, each checked equal to the one-thread answer.
// Not a CTest test.
//
// Usage: profile_parallelisosig '<PD code>' [threads ...]
//        (default thread counts: 1 2 4 6 8 12)

#include <chrono>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "surfer/pairsig/parallelisosig.h"
#include "diagramtriangulation/thickening/thickening.h"

int main(int argc, char **argv) {
    if (argc < 2) {
        std::cerr << "usage: profile_parallelisosig '<PD code>' [threads ...]\n";
        return 2;
    }
    std::vector<unsigned> counts;
    for (int i = 2; i < argc; ++i) counts.push_back(std::strtoul(argv[i], nullptr, 10));
    if (counts.empty()) counts = {1, 2, 4, 6, 8, 12};

    ThickenedLink rb;
    buildAmbient(argv[1], 2, 2, rb);
    std::cout << "thickening: " << rb.tri.size() << " pentachora\n";

    std::string reference;
    regina::Isomorphism<4> referenceIso(rb.tri.size());
    double one = 0;
    for (unsigned t : counts) {
        const auto t0 = std::chrono::steady_clock::now();
        auto [sig, iso] = parallelIsoSigDetail(rb.tri, t);
        const double s =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        if (reference.empty()) {
            reference = sig;
            referenceIso = iso;
        }
        if (t == 1) one = s;
        const bool same = sig == reference && iso == referenceIso;
        std::cout << t << " threads: " << s << " s"
                  << (one > 0 ? " (" + std::to_string(one / s) + "x)" : std::string())
                  << (same ? "" : "  MISMATCH") << "\n";
        if (!same) return 1;
    }
    return 0;
}

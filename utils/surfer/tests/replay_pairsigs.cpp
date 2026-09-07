//
//  replay_pairsigs.cpp
//
//  Regression check for PairSigContext against the RECORDED corpus.
//
//  The pair signatures in cobordism-atlas/results/cobordisms.csv were all
//  produced by the pre-context code path. This replays them through the new
//  one: decode each recorded signature, build ONE PairSigContext from the
//  decoded ambient (every witness of a row shares that ambient -- which is
//  the premise the whole optimization rests on), re-encode every witness
//  through it, and require the result to be byte-identical to what was
//  recorded.
//
//  That is a stronger statement than the unit tests can make. They check the
//  context against the free function on hand-built triangulations; this
//  checks it against signatures a DIFFERENT BUILD actually wrote to the
//  primary result file, on real search output.
//
//  Usage:  replay_pairsigs <file-of-pairsigs> [...]
//          (one signature per line; see tests' extraction step)
//
//  Reads only. Writes nothing.
//
#include <chrono>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <triangulation/dim4.h>

#include "../pairsig.h"

namespace {

using Clock = std::chrono::steady_clock;

double secsSince(Clock::time_point t0) {
    return std::chrono::duration<double>(Clock::now() - t0).count();
}

// Returns true if every signature in `path` round-trips byte-identically.
bool replayGroup(const std::string &path) {
    std::vector<std::string> sigs;
    {
        std::ifstream in(path);
        if (!in) {
            std::cerr << "[!] cannot open " << path << "\n";
            return false;
        }
        std::string line;
        while (std::getline(in, line))
            if (!line.empty())
                sigs.push_back(line);
    }
    if (sigs.empty()) {
        std::cerr << "[!] " << path << " is empty\n";
        return false;
    }

    // Decode every recorded signature first. Each yields its own ambient,
    // but all are fromSig() of the same isoSig, hence identically numbered --
    // which is what lets a single context serve the whole group, and is
    // exactly the property being relied on in the search.
    auto t0 = Clock::now();
    std::vector<std::vector<int>> marked;
    marked.reserve(sigs.size());
    std::unique_ptr<regina::Triangulation<4>> ambient;
    for (const std::string &s : sigs) {
        auto decoded = fromKnottedSurfaceSig(s);
        if (!ambient)
            ambient = std::move(decoded.ambient);
        marked.push_back(decoded.surface->markedFaces());
    }
    const double decodeSecs = secsSince(t0);

    t0 = Clock::now();
    PairSigContext<4, 2> ctx(*ambient);
    const double ctxSecs = secsSince(t0);

    t0 = Clock::now();
    size_t mismatches = 0;
    for (size_t i = 0; i < sigs.size(); ++i)
        if (ctx.sig(marked[i]) != sigs[i])
            ++mismatches;
    const double reencodeSecs = secsSince(t0);

    std::cout << (mismatches == 0 ? "  PASS  " : "  FAIL  ") << path << ": "
              << sigs.size() << " recorded signatures, " << mismatches
              << " mismatches\n"
              << "          decode " << decodeSecs << "s, context "
              << ctxSecs << "s, re-encode " << reencodeSecs << "s ("
              << (sigs.empty() ? 0.0
                               : reencodeSecs * 1000.0 /
                                     static_cast<double>(sigs.size()))
              << " ms/signature)\n"
              << "          automorphisms=" << ctx.automorphismCount()
              << "\n";
    return mismatches == 0;
}

} // namespace

int main(int argc, char *argv[]) {
    if (argc < 2) {
        std::cerr << "usage: " << argv[0] << " <file-of-pairsigs> [...]\n";
        return 1;
    }
    bool ok = true;
    for (int i = 1; i < argc; ++i)
        ok = replayGroup(argv[i]) && ok;
    std::cout << (ok ? "\nALL GROUPS REPRODUCED EXACTLY\n"
                     : "\nMISMATCHES PRESENT\n");
    return ok ? 0 : 1;
}

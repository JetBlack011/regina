//
//  bogocheck.cpp
//
//  Created by John Teague on 07/18/2024.
//
//  Scratch tool: retriangulates a hardcoded 4-manifold triangulation via
//  Pachner moves (see process()) and reports whether a strictly smaller
//  triangulation exists. The configuration globals below are set directly
//  in source; argc/argv are currently unused.
//

#include <link/link.h>

std::mutex mutex;

int argHeight = 2; // Pachner move search depth, passed to retriangulate()
long argThreads = 10; // worker thread count, passed to retriangulate()
enum {
    FLAVOUR_NONE = 0,
    FLAVOUR_DIM3 = 3,
    FLAVOUR_DIM4 = 4,
    FLAVOUR_KNOT = 100
} flavour = FLAVOUR_NONE; // currently unused
bool internalSig = false; // if true, print retriangulate()'s own signature instead of recomputing the classic isoSig

// Searches (up to argHeight Pachner moves, using argThreads worker
// threads) for a triangulation smaller than `tri` with the same topology,
// printing every isomorphism signature visited along the way. Reports
// whether a smaller triangulation was found, and if so, its signature.
template <int dim>
void process(const regina::Triangulation<dim>& tri) {
    unsigned long nSolns = 0;
    bool nonMinimal = false;
    std::string simpler;

    tri.retriangulate(
        argHeight, argThreads, nullptr /* tracker */,
        [&nSolns, &nonMinimal, &simpler, &tri](
            const std::string& sig, const regina::Triangulation<dim>& t) {
            if (t.size() > tri.size()) return false;

            if (internalSig) {
                std::lock_guard<std::mutex> lock(mutex);
                std::cout << sig << std::endl;

                if (t.size() < tri.size()) {
                    nonMinimal = true;
                    simpler = sig;
                    return true;
                }

                ++nSolns;
                return false;
            } else {
                // Recompute the signature using the default type IsoSigClassic.
                std::string classic = t.isoSig();

                std::lock_guard<std::mutex> lock(mutex);
                std::cout << classic << std::endl;

                if (t.size() < tri.size()) {
                    nonMinimal = true;
                    simpler = std::move(classic);
                    return true;
                }

                ++nSolns;
                return false;
            }
        });

    if (nonMinimal) {
        std::cerr << "Triangulation is non-minimal!" << std::endl;
        std::cerr << "Smaller triangulation: " << simpler << std::endl;
    } else
        std::cerr << "Found " << nSolns << " triangulation(s)." << std::endl;
}

int main(int argc, char* argv[]) {
    // t1/t2 are currently unused -- kept here for quick swapping into the
    // process() call below when trying a different triangulation.
    regina::Triangulation<4> t1(
        "6LvLLwvzzzQwvwzAPPMMLwMQPQMMAPLAvQMQQwQzQQLQQQQQQcbfigmqtmmlqzBBtpwypx"
        "pCxKALLKEwCMOMDDPJJPKISUUEYKROOMSIIYMYU1XX32SX1155WVZVZ3X22053354444aa"
        "aaaaaaaaaaaabababaaaaaaaaaaa2aObuagauabayaqbaayaaaaaJbbaIbyakbTakbyaya"
        "aaaaaaaaaakbaaaaaaNaaaJbJbaaaa3bEbTaTaJboaNaoaaaPbPbTaMayaPbaaaaaakbNa"
        "qbPbqbPbrbNaJbJbJbJbaaaaoaoaoaoaoa");
    regina::Triangulation<4> t2(
        "-cqbLLLLAvPPLPLLMPzLzvPLwLQPwzLvwzMQMwPLzPQQQQPQQQLLMwAAQLQQQQAQQPAAQQ"
        "kbaeafafaiaialalaoamaqajaratavaualaxanayayaCaGaGaJaIaNaNaIayaJaOaOaNaU"
        "aXa0aTa0aUa2a2aVaWaWa6aZa1a9aZaZaZa8a8a8aUa2a5a5aWaRa7a7a6a6aTa8aTaTaU"
        "a0a3a3a0acbebdbabab4ahbhbhbdbibibibgbjbjbababcbcb+"
        "a9aebbbfbfbebmbeblbgbgblbnbnbjbmbobobobobpbnbpbpbnbpbaabaaaaaaaaababaa"
        "afaaabaaaaaaaaabaaacaaabaaaaaaaaaaaaaaaaabaaaaaaaeaaaaaaabaaabaaaaabab"
        "abaaaaaaaaaaaaaaaaaaaaabaaaaaaababaaaaaaaaabaaabababaaaaaaaaaaaaaaaeae"
        "afababababaaaaaaabaaaaafafababafabacabacacabaaabaaababaaaaaaaeaeaaabab"
        "abaaabaaaaabaaa");

    process(regina::Triangulation<4>::fromIsoSig(" YLLLvLLvvALzQzMLAwvzLAQzAzQQLQQAQQQQQAQQQQkccflmrrwnywomxxoCyCJMKNLzLNNGRAMARAKFJIIEEMUUSJVRVSSTWWLOONPPRUOXXXTPWWSTTVXaayaaaaaaaaaaaaavaaaMaSaSaaaaabaEbaaEbaaaaaaaaMaqbMaaaaaSaaaqbkbqbaaqbaavaMakbkblblbkbaaaakbaaaaaaaaTaTaTaEbEbparbrbpapapapaoaoaMaMaaayayaqbqbyayayaqbMa"));
}

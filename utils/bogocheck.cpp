#include <gmpxx.h>
#include <link/link.h>
#include <triangulation/dim2.h>
#include <triangulation/dim3.h>
#include <triangulation/dim4.h>
#include <unistd.h>
#include "triangulation/generic/triangulation.h"

std::mutex mutex;

int argHeight = 2;
long argThreads = 10;
enum {
    FLAVOUR_NONE = 0,
    FLAVOUR_DIM3 = 3,
    FLAVOUR_DIM4 = 4,
    FLAVOUR_KNOT = 100
} flavour = FLAVOUR_NONE;
bool internalSig = false;

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
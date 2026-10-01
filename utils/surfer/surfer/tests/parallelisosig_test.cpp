// parallelisosig_test.cpp
//
// parallelIsoSigDetail() (../parallelisosig.h) must be isoSigDetail() byte
// for byte: the same signature AND the same isomorphism onto the canonical
// form (a pair-signature context stores that isomorphism, pairsig.h), at
// every thread count, including counts that do not divide the number of
// starts and counts above it.
//
// Inputs: campaign-shaped row thickenings built from real PD codes (the
// ambients pair signatures are taken over), each also with its labelling
// randomised (a different serial winner, and ties between starts whenever
// the triangulation has automorphisms), and a few of Regina's 3- and
// 4-dimensional examples.

#include <iostream>
#include <string>
#include <vector>

#include <triangulation/dim3.h>
#include <triangulation/dim4.h>
#include <triangulation/example3.h>
#include <triangulation/example4.h>

#include "surfer/pairsig/parallelisosig.h"
#include "diagramtriangulation/thickening/thickening.h"

static int passed = 0, failed_count = 0;

namespace {

const std::vector<unsigned> allCounts = {1, 2, 3, 5, 8, 16};

template <int dim>
void check(const std::string &name, const regina::Triangulation<dim> &tri,
           const std::vector<unsigned> &counts = allCounts) {
    const auto serial = tri.isoSigDetail();
    for (unsigned threads : counts) {
        const auto par = parallelIsoSigDetail(tri, threads);
        if (par.first != serial.first) {
            std::cout << "  FAIL: " << name << " at " << threads
                      << " threads: signature differs\n";
            ++failed_count;
        } else if (!(par.second == serial.second)) {
            std::cout << "  FAIL: " << name << " at " << threads
                      << " threads: isomorphism differs\n";
            ++failed_count;
        } else {
            ++passed;
        }
    }
}

template <int dim>
void checkWithRelabellings(const std::string &name,
                           const regina::Triangulation<dim> &tri,
                           const std::vector<unsigned> &counts = allCounts,
                           int relabellings = 2) {
    check(name, tri, counts);
    for (int k = 0; k < relabellings; ++k) {
        regina::Triangulation<dim> copy(tri);
        copy.randomiseLabelling(false);
        check(name + " relabelled " + std::to_string(k), copy, counts);
    }
}

} // namespace

int main() {
    std::cout << "parallelIsoSigDetail() against isoSigDetail()\n";

    const std::vector<std::pair<std::string, std::string>> rows = {
        {"L2a1{0}", "PD[X[4; 1; 3; 2]; X[2; 3; 1; 4]]"},
        {"L4a1{0}", "PD[X[6; 1; 7; 2]; X[8; 3; 5; 4]; X[2; 5; 3; 6]; "
                    "X[4; 7; 1; 8]]"},
        {"3_1", "[[1;5;2;4];[3;1;4;6];[5;3;6;2]]"},
        {"4_1", "[[4;2;5;1];[8;6;1;5];[6;3;7;4];[2;7;3;8]]"},
    };
    // Every row's T at every count; one row's thickening (the ambient a
    // pair signature is taken over, 384 pentachora) at a few, since each
    // serial isoSigDetail() of one takes seconds and ctest allows 60.
    for (const auto &[name, pd] : rows) {
        ThickenedLink rb;
        buildAmbient(pd, 2, 2, false, rb);
        checkWithRelabellings(name + " T", rb.link.tri);
        if (name == "L2a1{0}")
            checkWithRelabellings(name + " thickening", rb.tri, {3, 8}, 1);
    }

    checkWithRelabellings("S4", regina::Example<4>::sphere());
    checkWithRelabellings("CP2", regina::Example<4>::cp2());
    checkWithRelabellings("Weber-Seifert", regina::Example<3>::weberSeifert());
    checkWithRelabellings("Poincare", regina::Example<3>::poincare());
    checkWithRelabellings("simplicial S4", regina::Example<4>::simplicialSphere());

    std::cout << passed << " passed, " << failed_count << " failed\n";
    return failed_count ? 1 : 0;
}

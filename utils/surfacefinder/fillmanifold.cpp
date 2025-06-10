//
//  fillmanifold.cpp
//
//  Created by John Teague on 04/12/2025.

#include <gmpxx.h>
#include <link/link.h>
#include <triangulation/dim2.h>
#include <triangulation/dim3.h>
#include <triangulation/dim4.h>
#include <unistd.h>

#include <iostream>
#include <ostream>
#include <string>

#include "triangulation/example3.h"
#include "triangulation/forward.h"
#include "triangulation/generic/boundarycomponent.h"
#include "triangulation/generic/isomorphism.h"
#include "triangulation/generic/triangulation.h"

#include "knotbuilder.h"

namespace {
void usage(const char *progName, const std::string &error = std::string()) {
    if (!error.empty())
        std::cerr << error << "\n\n";

    std::cerr << "Usage:\n";
    std::cerr << "    " << progName << " {dim} <isosig>\n"
              << "    dim  : Specify the dimension of the given isomorphism "
                 "signature (optional, default is 3)\n"
              << "    " << progName << " [ -v, --version | -h, --help ]\n\n";
    std::cerr << "    -v, --version  : Show which version of Regina is "
                 "being used\n";
    std::cerr << "    -h, --help     : Display this help\n";
    exit(1);
}
} // namespace

int main(int argc, char *argv[]) {
    // Check for standard arguments:
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-h" || arg == "--help")
            usage(argv[0]);
        if (arg == "-v" || arg == "--version") {
            if (argc != 2)
                usage(argv[0], "Option --version cannot be used with "
                               "any other arguments.");
            std::cout << PACKAGE_BUILD_STRING << "\n";
            exit(0);
        }
    }

    std::string isoSig;
    int dim = 3;

    if (argc == 2) {
        isoSig = argv[1];
    } else if (argc == 3) {
        isoSig = argv[2];
        dim = std::stoi(argv[1]);
    } else {
        usage(argv[0],
              "Please provide an isomorphism signature and/or dimension.");
    }

    knotbuilder::PDCode pdcode;

    regina::Triangulation<3> threeMfld = regina::Example<3>::sphere600();
    std::cout << "Original = " << threeMfld.isoSig() << "\n";
    regina::Triangulation<4> tri1 = knotbuilder::thicken(threeMfld).value();
    regina::Triangulation<4> tri2 = knotbuilder::cone(threeMfld);
    // regina::Triangulation tri = regina::Example<2>::sphere();
    // tri.subdivide();
    //regina::Triangulation<2> sphere = regina::Example<2>::sphere();
    //regina::Triangulation<3> tri1 = thicken(sphere).value();
    //tri1.subdivide();
    //regina::Triangulation<3> tri2(tri1);

    //std::cout << "Tri1 = " << tri1.isoSig() << "\n";

    // Find all isomorphisms between the two triangulations
    std::vector<regina::Isomorphism<3>> isos;
    tri1.boundaryComponent(0)->build().findAllIsomorphisms(
        tri2.boundaryComponent(0)->build(),
        [&isos](const regina::Isomorphism<3> &iso) {
            isos.push_back(iso);
            return true;
        });
    std::cout << "[+] Found " << isos.size() << " isomorphisms\n";
    for (const auto &iso : isos) {
        std::cout << iso << "\n";
    }
    std::cout << "\n";
    regina::Triangulation<4> tri = knotbuilder::glue(tri1, 0, tri2, 0, isos[0]);

    //std::cout << "Tri = " << tri.isoSig() << "\n";
    //if (tri.isValid()) {
    //    std::cout << "[+] Triangulation is valid\n";
    //} else {
    //    std::cout << "[!] Triangulation is NOT valid\n";
    //}

    //isos.clear();
    //tri.boundaryComponent(0)->build().findAllIsomorphisms(
    //    tri2.boundaryComponent(0)->build(),
    //    [&isos](const regina::Isomorphism<3> &iso) {
    //        isos.push_back(iso);
    //        return true;
    //    });
    //std::cout << "[+] Found " << isos.size() << " isomorphisms\n";
    //for (const auto &iso : isos) {
    //    std::cout << iso << "\n";
    //}
    //std::cout << "\n";
    //regina::Triangulation<4> thickerTri = glue(tri, 0, tri2, 0, isos[0]);

    //std::cout << "Tri = " << thickerTri.isoSig() << "\n";
    //if (tri.isValid()) {
    //    std::cout << "[+] Triangulation is valid\n";
    //} else {
    //    std::cout << "[!] Triangulation is NOT valid\n";
    //}

    // auto conedTri = cone(tri);

    // std::cout << "[+] Number of simplices = " << tri.simplices().size() <<
    // "\n"; std::cout << "[+] Number of cone simplices = "
    //           << conedTri.simplices().size() << "\n";

    // std::vector<regina::Isomorphism<2>> isos;

    // auto res = thicken(tri);
    // if (!res.has_value()) {
    //     std::cout << "[!] Given triangulation could not be thickened! Make "
    //                  "sure it's possible to order it (e.g. make sure the link
    //                  " "diagram you gave to knotbuilder is alternating).\n";
    //     return -1;
    // }
    // const auto &thickTri = res.value();
    // regina::Triangulation evenThickerTri(thickTri);

    // evenThickerTri.boundaryComponent(0)->build().findAllIsomorphisms(
    //     thickTri.boundaryComponent(1)->build(),
    //     [&isos](const regina::Isomorphism<2> &iso) {
    //         isos.push_back(iso);
    //         return true;
    //     });

    // for (int i = 0; i < 2; ++i) {
    //     std::cout << isos[0] << "\n";
    //     evenThickerTri = glue(evenThickerTri, 0, thickTri, 1, isos[0]);
    // }

    // isos.clear();
    // evenThickerTri.boundaryComponent(0)->build().findAllIsomorphisms(
    //     conedTri.boundaryComponent(0)->build(),
    //     [&isos](const regina::Isomorphism<2> &iso) {
    //         isos.push_back(iso);
    //         return false;
    //     });

    //// regina::Triangulation<4> thickenedTri =
    //// regina::Example<4>::iBundle(tri);

    // std::cout << "[+] Number of thickened simplices = "
    //           << thickTri.simplices().size() << "\n";
    // std::cout << "[+] Number of even more thickened simplices = "
    //           << evenThickerTri.simplices().size() << "\n";
    // std::cout << "[+] Total simplices after gluing = "
    //           << conedTri.simplices().size() + thickTri.simplices().size()
    //           << "\n\n";
    // if (thickTri.isValid()) {
    //     std::cout << "[+] Thickened triangulation is valid\n";
    // } else {
    //     std::cout << "[!] Thickened triangulation is NOT valid\n";
    // }
    // if (evenThickerTri.isValid()) {
    //     std::cout << "[+] Even thicker triangulation is valid\n";
    // } else {
    //     std::cout << "[!] Even thicker triangulation is NOT valid\n";
    // }

    // std::cout << "[+] Cone isosig = " << conedTri.isoSig() << "\n\n";
    // std::cout << "[+] Thickened isosig = " << thickTri.isoSig() << "\n\n";
    // std::cout << "[+] Even thicker isosig = " << evenThickerTri.isoSig()
    //           << "\n\n";

    // for (const auto &iso : isos) {
    //     auto gluedTri = glue(evenThickerTri, 0, conedTri, 0, iso);

    //    if (gluedTri.isValid()) {
    //        std::cout << "[+] Thickened + coned triangulation is valid\n";
    //    } else {
    //        std::cout << "[!] Thickened + coned triangulation is NOT valid\n";
    //    }
    //    std::cout << "ISO SIG = " << gluedTri.isoSig() << "\n";
    //}
    return 0;
}

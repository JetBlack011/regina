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

    regina::Triangulation<3> bdryTri(isoSig);
    knotbuilder::CobordismBuilder<3> cob(bdryTri);
    regina::Triangulation<4> tri = cob.thicken(2);

    std::cout << "[+] Thickened isosig = " << tri.isoSig()
              << "\n";

    //regina::Triangulation<3> threeMfld = regina::Example<3>::sphere600();
    //std::cout << "Original = " << threeMfld.isoSig() << "\n";
    //regina::Triangulation<4> tri1 = knotbuilder::thicken(threeMfld).value();
    //regina::Triangulation<4> tri2 = knotbuilder::cone(threeMfld);

    return 0;
}

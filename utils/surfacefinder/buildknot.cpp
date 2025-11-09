#include <gmpxx.h>
#include <link/link.h>
#include <triangulation/dim2.h>
#include <triangulation/dim3.h>
#include <triangulation/dim4.h>
#include <unistd.h>

#include <cassert>
#include <iostream>
#include <string>

#include "triangulation/example2.h"
#include "triangulation/forward.h"
#include "triangulation/generic/triangulation.h"

#include "knotbuilder.h"
#include "knottedsurfaces.h"

namespace {
void usage(const char *progName, const std::string &error = std::string()) {
    if (!error.empty())
        std::cerr << error << "\n\n";

    std::cerr << "Usage:\n";
    std::cerr << "    " << progName << " <PD Code>\n    " << progName
              << " [ -v, --version | -h, --help ]\n\n";
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

    if (argc != 2) {
        usage(argv[0], "Please provide a valid PD Code.");
    }

    //regina::Triangulation<4> tri;
    //auto p1 = SimplicialPrism(tri);
    //auto p2 = SimplicialPrism(tri);
    //auto p3 = SimplicialPrism(tri);
    //p1.glue(0, p2, 0);
    //p2.glue(1, p3, 0);
    //std::cout << "Triangulation = " << tri.isoSig() << "\n";

    regina::Triangulation<3> tri;
    std::vector<const regina::Edge<3> *> edges;
    knotbuilder::buildLink(
        tri, knotbuilder::parsePDCode("	[[4,2,5,1],[8,6,1,5],[6,3,7,4],[2,7,3,8]]"),
        edges);
    knotbuilder::CobordismBuilder cob(tri);
    regina::Triangulation<3> copyTri = tri;
    copyTri.finiteToIdeal();
    Link l(copyTri, edges);
    std::cout << "Triangulation of the link = " << tri.isoSig() << "\n";
    std::cout << "Triangulation of the complement = "
              << l.buildComplement().isoSig() << "\n";
    regina::Triangulation<4> &thickenedTri = cob.thicken();
    std::cout << "Thickened triangulation = " << thickenedTri.isoSig() << "\n";

    auto edgesTimesI = cob.facesTimesI<1>(edges);

    int counter = 0;
    std::cout << "Edges times I:\n";
    for (int i = 0; i < edges.size(); ++i) {

        std::cout << "Edge " << i << " = " << edges[i]->front().vertices()
                  << "\n";
        for (int j = 2 * i; j < 2 * (i + 1); ++j) {
            std::cout << "Num embeddings = "
                      << edgesTimesI[j]->embeddings().size() << "\n";
            const auto &emb = edgesTimesI[j]->front();
            std::cout << "Simplex = " << emb.simplex()->index()
                      << ", face = " << emb.face()
                      << ", Vertices = " << emb.vertices();
            std::cout << "\n";
        }
        std::cout << "\n";
    }

    // regina::Triangulation<3> tri;
    // std::vector<const regina::Edge<3> *> edges;

    // std::string pdcode_str = argv[1];
    // knotbuilder::PDCode pdcode = knotbuilder::parsePDCode(pdcode_str);

    // std::cout << "[*] Building link with " << pdcode.size()
    //           << " crossings.\n";

    // if (pdcode_str.length() <= 4 || !knotbuilder::buildLink(tri, pdcode,
    // edges)) {
    //     usage(argv[0], "Please provide a valid PD Code.");
    // }

    // std::cout << "[+] Number of edges in the link = " << edges.size() <<
    // "\n"; std::cout << "\nTriangulation of S^3 = " << tri.isoSig() << "\n\n";

    // Link l(tri, edges);
    ////std::cout << "Triangulation of the complement = "
    ////          << l.buildComplement().isoSig() << "\n";

    // knotbuilder::CobordismBuilder<3> cob(tri);
    // regina::Triangulation<4> thickenedTri = cob.thicken();

    // std::cout << "[+] Thickened triangulation = " << thickenedTri.isoSig()
    //           << "\n";

    // cob.edgesTimesI(edges);

    return 0;
}

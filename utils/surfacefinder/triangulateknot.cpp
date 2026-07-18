#include "cobordismbuilder.h"
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

    std::string pdcode_str = argv[1];
    knotbuilder::PDCode pdcode = knotbuilder::parsePDCode(pdcode_str);

    std::cout << "[*] Building link with " << pdcode.size() << " crossings.\n";

    if (pdcode_str.length() <= 4) {
        usage(argv[0], "Please provide a valid PD Code.");
    }

    knotbuilder::TriangulationWithLink result;
    try {
        result = knotbuilder::buildLink(pdcode);
    } catch (const regina::InvalidArgument &e) {
        usage(argv[0], "Please provide a valid PD Code.");
    }
    auto &[tri, edges] = result;

    std::cout << "[+] Number of edges in the link = " << edges.size() << "\n";
    std::cout << "\nTriangulation of S^3 = " << tri.isoSig() << "\n\n";

    Link l(tri, edges);
    std::cout << "Triangulation of the complement = "
              << l.buildComplement().isoSig() << "\n";

    CobordismBuilder<3> cob(tri);
    regina::Triangulation<4> thickenedTri = cob.thicken();

    std::cout << "[+] Thickened triangulation = " << thickenedTri.isoSig()
              << "\n";

    return 0;
}

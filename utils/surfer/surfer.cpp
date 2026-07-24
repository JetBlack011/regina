//
//  surfer.cpp
//
//  Created by John Teague on 06/19/2024.
//

#include <thread>

#include <triangulation/dim4.h>

#include "cobordismbuilder.h"
#include "collar.h"
#include "embeddingsearch.h"
#include "knotbuilder.h"

namespace {
void usage(const char *progName, const std::string &error = std::string()) {
  if (!error.empty())
    std::cerr << error << "\n\n";

  std::cerr << "Usage:\n";
  std::cerr
      << "    " << progName
      << " [ -a, --all | -c, --closed | -p, --proper | --connected ]\n"
         "    [ --threads N ] <isosig>\n\n"
      << "    " << progName
      << " [ -a, --all | -c, --closed | -p, --proper | --connected ]\n"
         "    [ --threads N ] [ --layers N ] [ --cone | --no-cone ]\n"
         "    [ --collar | --no-collar ] --pd <pdcode>\n\n"
      << "    " << progName << " [ -v, --version | -h, --help ]\n\n";
  std::cerr
      << "    -a, --all      : Find all embedded submanifolds, regardless of "
         "boundary\n"
         "                     conditions (default)\n";
  std::cerr << "    -c, --closed   : Find only closed embedded submanifolds\n";
  std::cerr << "    -p, --proper   : Find embedded submanifolds whose "
               "boundary is contained\n"
               "                     entirely in the boundary of the ambient "
               "triangulation\n";
  std::cerr << "    --connected    : Same as --proper, but additionally "
               "require at most one\n"
               "                     boundary component of the embedded "
               "submanifold per\n"
               "                     boundary component of the ambient "
               "triangulation\n\n";
  std::cerr << "    --threads N    : Number of worker threads to search "
               "with (default: the\n"
               "                     number of hardware threads available)\n"
               "\n";
  std::cerr
      << "    <isosig>       : Isomorphism signature of a 4-manifold\n"
         "                     triangulation to search directly (default "
         "input mode)\n\n";
  std::cerr
      << "    --pd <pdcode>  : Planar diagram code of a knot/link. Builds a "
         "3-sphere\n"
         "                     triangulation containing it, thickens it "
         "into a\n"
         "                     4-dimensional cobordism, and searches that "
         "instead of\n"
         "                     a bare isosig. Mutually exclusive with "
         "<isosig>.\n";
  std::cerr << "    --layers N     : Number of times to thicken() the "
               "cobordism (default: 1;\n"
               "                     only valid with --pd)\n";
  std::cerr
      << "    --cone, --no-cone      : Whether to cap the cobordism with "
         "cone()\n"
         "                     (default: --cone; only valid with --pd)\n";
  std::cerr
      << "    --collar, --no-collar  : Whether to seed the search with a "
         "collar\n"
         "                     traced through the link's edges (default: "
         "--collar;\n"
         "                     only valid with --pd)\n\n";
  std::cerr
      << "    -v, --version  : Show which version of Regina is being used\n";
  std::cerr << "    -h, --help     : Display this help\n";
  exit(1);
}

// Runs a SurfaceSearch over `tri`, seeded with `seedFaces` if non-empty
// (an empty seed is not a valid seeded-search argument, so this falls back
// to an unseeded search instead), and reports the same way regardless of
// which input mode built `tri`.
void runSearch(const regina::Triangulation<4> &tri,
              const std::vector<int> &seedFaces, BoundaryCondition cond,
              unsigned numThreads) {
  std::cerr << "[+] Running with " << numThreads
            << " threads, condition = " << boundaryConditionName(cond)
            << "\n\n";

  if (seedFaces.empty()) {
    SurfaceSearch e(tri);
    e.search(numThreads, cond);
  } else {
    SurfaceSearch e(tri, seedFaces);
    e.search(numThreads, cond);
  }
}
} // namespace

int main(int argc, char *argv[]) {
  BoundaryCondition cond = BoundaryCondition::all;
  unsigned numThreads = std::thread::hardware_concurrency();
  if (numThreads == 0)
    numThreads = 1;

  std::string isoSig;
  bool haveIsoSig = false;

  std::string pdCode;
  bool havePD = false;

  int layers = 1;
  bool useCone = true;
  bool useCollar = true;
  bool sawLayers = false, sawCone = false, sawCollar = false;

  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "-h" || arg == "--help") {
      usage(argv[0]);
    } else if (arg == "-v" || arg == "--version") {
      if (argc != 2)
        usage(argv[0], "Option --version cannot be used with "
                       "any other arguments.");
      std::cout << PACKAGE_BUILD_STRING << "\n";
      return 0;
    } else if (arg == "-a" || arg == "--all") {
      cond = BoundaryCondition::all;
    } else if (arg == "-c" || arg == "--closed") {
      cond = BoundaryCondition::closed;
    } else if (arg == "-p" || arg == "--proper") {
      cond = BoundaryCondition::proper;
    } else if (arg == "--connected") {
      cond = BoundaryCondition::connected;
    } else if (arg == "--threads") {
      if (i + 1 >= argc)
        usage(argv[0], "--threads requires a value.");
      try {
        numThreads = static_cast<unsigned>(std::stoul(argv[++i]));
      } catch (const std::exception &) {
        usage(argv[0], "--threads requires an integer value.");
      }
    } else if (arg == "--pd") {
      if (i + 1 >= argc)
        usage(argv[0], "--pd requires a value.");
      pdCode = argv[++i];
      havePD = true;
    } else if (arg == "--layers") {
      if (i + 1 >= argc)
        usage(argv[0], "--layers requires a value.");
      try {
        layers = std::stoi(argv[++i]);
      } catch (const std::exception &) {
        usage(argv[0], "--layers requires an integer value.");
      }
      sawLayers = true;
    } else if (arg == "--cone") {
      useCone = true;
      sawCone = true;
    } else if (arg == "--no-cone") {
      useCone = false;
      sawCone = true;
    } else if (arg == "--collar") {
      useCollar = true;
      sawCollar = true;
    } else if (arg == "--no-collar") {
      useCollar = false;
      sawCollar = true;
    } else if (!arg.empty() && arg[0] == '-') {
      usage(argv[0], "Unknown option: " + arg);
    } else if (haveIsoSig) {
      usage(argv[0], "Only one isosig may be specified.");
    } else {
      isoSig = arg;
      haveIsoSig = true;
    }
  }

  if (havePD && haveIsoSig)
    usage(argv[0], "--pd and an isosig are mutually exclusive.");
  if (!havePD && !haveIsoSig)
    usage(argv[0], "Please specify an isosig or --pd <pdcode>.");
  if (!havePD && (sawLayers || sawCone || sawCollar))
    usage(argv[0], "--layers, --cone/--no-cone, and --collar/--no-collar "
                   "are only valid with --pd.");

  std::cout << "------ SurFer (Surface Finder) \U0001F30A ------\n\n";

  if (havePD) {
    knotbuilder::PDCode pdcode = knotbuilder::parsePDCode(pdCode);

    knotbuilder::TriangulationWithLink link;
    try {
      link = knotbuilder::buildLink(pdcode);
    } catch (const regina::InvalidArgument &e) {
      usage(argv[0], std::string("Invalid PD code: ") + e.what());
    }
    auto &[t2, edges2] = link;

    // Indices are preserved across CobordismBuilder's internal copy/reorder
    // of t2 (see CobordismBuilder::baseTriangulation()'s doc comment), so
    // edgeIndices (taken from t2) still identify the same edges in
    // cob.baseTriangulation().
    std::vector<int> edgeIndices;
    edgeIndices.reserve(edges2.size());
    for (const regina::Edge<3> *e : edges2)
      edgeIndices.push_back(static_cast<int>(e->index()));

    CobordismBuilder<3> cob(t2);
    CollarBuilder collarBuilder(edgeIndices);
    for (int i = 0; i < layers; ++i) {
      cob.thicken();
      // Must run every layer, not just the first: CollarBuilder::addLayer
      // only captures the most-recently-built thickening layer's prisms, so
      // tracing the collar all the way from the link's original position up
      // through every layer requires calling it once per thicken().
      if (useCollar)
        collarBuilder.addLayer(cob);
    }
    if (useCone)
      cob.cone();
    regina::Triangulation<4> tri = cob.getCobordism();

    std::vector<int> seedFaces;
    if (useCollar) {
      for (regina::Triangle<4> *t : collarBuilder.resolve())
        seedFaces.push_back(static_cast<int>(t->index()));
      std::cerr << "[+] Collar seed faces = " << seedFaces.size() << "\n";
    }

    runSearch(tri, seedFaces, cond, numThreads);
  } else {
    regina::Triangulation<4> tri;
    try {
      tri = regina::Triangulation<4>(isoSig);
    } catch (const regina::InvalidArgument &e) {
      usage(argv[0], std::string("Invalid isosig: ") + e.what());
    }

    runSearch(tri, {}, cond, numThreads);
  }

  return 0;
}

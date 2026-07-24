//
//  surfer.cpp
//
//  Created by John Teague on 06/19/2024.
//

#include <maths/perm.h>
#include <triangulation/dim4.h>
#include <triangulation/example2.h>
#include <triangulation/example3.h>
#include <triangulation/example4.h>

#include "cobordismbuilder.h"
#include "collar.h"
#include "embeddingsearch.h"
#include "knotbuilder.h"
#include "skeleton.h"

namespace {
void usage(const char *progName, const std::string &error = std::string()) {
    if (!error.empty())
        std::cerr << error << "\n\n";

    std::cerr << "Usage:\n";
    std::cerr << "    " << progName
              << " [ -a, --all | -c, --closed | -p, --proper | --connected ]\n"
                 "    "
              << progName << " [ -v, --version | -h, --help ]\n\n";
    std::cerr
        << "    -a, --all      : Find all embedded submanifolds, regardless of "
           "boundary\n                     conditions (default)\n";
    std::cerr
        << "    -c, --closed   : Find only closed embedded submanifolds\n";
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
    std::cerr
        << "    -v, --version  : Show which version of Regina is being used\n";
    std::cerr << "    -h, --help     : Display this help\n";
    exit(1);
}

// template <int dim>
// void surfacesDetail(std::set<KnottedSurface<dim>> &surfaces,
//                     BoundaryCondition cond) {
//     std::cout << "--- "
//               << (cond == BoundaryCondition::all
//                       ? ""
//                       : (cond == BoundaryCondition::boundary
//                              ? "Proper "
//                              : (cond == BoundaryCondition::closed ? "Closed "
//                                                                  : "")))
//               << "Surfaces ---\n";
//     int surfaceCount = 0;
//     int closedCount = 0;
//
//     std::map<std::string, int> countMap;
//
//     for (const auto &surface : surfaces) {
//         ++surfaceCount;
//         if (surface.surface().isClosed()) {
//             ++closedCount;
//         }
//
//         // Sanity check: isProper() should agree with a direct scan of the
//         // surface's boundary image. This should never fire; if it does,
//         // it means isProper()'s bookkeeping has drifted from reality.
//         bool matchesIsProper = true;
//         for (const regina::BoundaryComponent<2> *comp :
//              surface.surface().boundaryComponents()) {
//             for (const regina::Edge<2> *edge : comp->edges()) {
//                 if (edge->isBoundary() && !surface.image(edge)->isBoundary())
//                 {
//                     matchesIsProper = false;
//                 }
//             }
//         }
//         if (!matchesIsProper) {
//             std::cout << "NOTE: " << surface.detail() << " is not
//             proper!!\n";
//         }
//
//         auto search = countMap.find(surface.detail());
//         if (search == countMap.end()) {
//             countMap.insert({surface.detail(), 1});
//         } else {
//             ++search->second;
//         }
//     }
//
//     std::cout << "\n";
//
//     std::vector<std::string> descriptions;
//     descriptions.reserve(countMap.size());
//     for (const auto &[description, count] : countMap) {
//         descriptions.push_back(description);
//     }
//
//     std::sort(descriptions.begin(), descriptions.end(),
//               [](const std::string &first, const std::string &second) {
//                   return first.size() < second.size() ||
//                          (first.size() == second.size() && first < second);
//               });
//
//     for (std::string &description : descriptions) {
//         int count = countMap.find(description)->second;
//         std::cout << "- " << count << (count == 1 ? " copy of " : " copies of
//         ")
//                   << description << "\n";
//     }
//
//     std::cout << "\n";
//     if (cond != BoundaryCondition::closed) {
//         std::cout << "[+] Total admissible surfaces = " << surfaceCount <<
//         "\n";
//     }
//     std::cout << "[+] Total closed surfaces = " << closedCount << "\n\n";
// }
} // namespace

int main(int argc, char *argv[]) {
    BoundaryCondition cond = BoundaryCondition::all;

    // Check for standard arguments:
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-h" || arg == "--help")
            usage(argv[0]);
        if (arg == "-v" || arg == "--version") {
            if (argc != 2)
                usage(argv[0], "Option --version cannot be used with "
                               "any other arguments.");
        }
        if (arg == "-a" || arg == "--all")
            cond = BoundaryCondition::all;
        else if (arg == "-c" || arg == "--closed")
            cond = BoundaryCondition::closed;
        else if (arg == "-p" || arg == "--proper")
            cond = BoundaryCondition::proper;
        else if (arg == "--connected")
            cond = BoundaryCondition::connected;
    }

    std::string pdcode_str = "[[4,2,5,1],[8,6,1,5],[6,3,7,4],[2,7,3,8]]";
    knotbuilder::PDCode pdcode = knotbuilder::parsePDCode(pdcode_str);
    auto [t2, edges0] = knotbuilder::buildLink(pdcode);
    Link l(t2, edges0);

    std::cout << "Link complement = " << l.buildComplement().isoSig() << "\n";
    // auto [t, edges] = knotbuilder::reduceVertices(t2, edges0);
    // std::cout << "Reduced vertices isosig = " << t.isoSig() << "\n";

    // Indices are preserved across CobordismBuilder's internal copy/reorder
    // of t2 (see CobordismBuilder::baseTriangulation()'s doc comment), so
    // edges0's indices (taken from t2) still identify the same edges in
    // cob.baseTriangulation().
    std::vector<int> edgeIndices;
    edgeIndices.reserve(edges0.size());
    for (const regina::Edge<3> *e : edges0)
        edgeIndices.push_back(static_cast<int>(e->index()));

    constexpr int numLayers = 1;
    constexpr int seedLayers = 1;
    CobordismBuilder<3> cob(t2);
    CollarBuilder collarBuilder(edgeIndices);
    for (int i = 0; i < numLayers; ++i) {
        cob.thicken();
        if (i < seedLayers)
            collarBuilder.addLayer(cob);
    }
    regina::Triangulation<4> tri = cob.getCobordism();
    //= cob.cone();
    // std::cout << tri.isoSig() << "\n";
    std::cerr << "Num triangles = " << tri.countTriangles() << "\n";

    std::vector<int> seedFaces;
    for (regina::Triangle<4> *t : collarBuilder.resolve())
        seedFaces.push_back(static_cast<int>(t->index()));
    std::cerr << "Collar seed faces = " << seedFaces.size() << "\n";

    Skeleton<4, 2> skel(tri);
    ////std::cout << skel << "\n";

    unsigned numThreads = 132;
    // unsigned numThreads = tri.countTriangles();
    if (numThreads == 0)
        numThreads = 1;
    std::cerr << "Running with " << numThreads
              << " threads, condition = " << boundaryConditionName(cond)
              << "\n\n";
    SurfaceSearch e(tri, seedFaces);
    e.search(numThreads, cond);

    // regina::Triangulation<4> tri;
    // BoundaryCondition cond = BoundaryCondition::boundary;
    //// std::vector<regina::Triangle<4> *> startingTriangles;

    // if (argc == 2) {
    //     regina::Triangulation<3> threeMfld;
    //     std::string pdcode_str = argv[1];
    //     knotbuilder::PDCode pdcode = knotbuilder::parsePDCode(pdcode_str);
    //     std::vector<const regina::Edge<3> *> linkEdges;
    //     knotbuilder::buildLink(threeMfld, pdcode, linkEdges);
    //     Link bdry(threeMfld, linkEdges);

    //    CobordismBuilder<3> cob(threeMfld);

    //    tri = cob.thicken();
    //} else if (argc == 3) {
    //    std::string arg = argv[1];
    //    std::string isoSig = argv[2];

    //    if (arg == "-a" || arg == "--all") {
    //        cond = BoundaryCondition::all;
    //    } else if (arg == "-b" || arg == "--boundary" || arg == "-l" ||
    //               arg == "--links") {
    //        cond = BoundaryCondition::boundary;
    //    } else if (arg == "-c" || arg == "--closed") {
    //        cond = BoundaryCondition::closed;
    //    } else {
    //        usage(argv[0], "Please specify a valid surface condition.");
    //    }

    //    tri = regina::Triangulation<4>(isoSig);
    //} else {

    //    usage(argv[0], "Please specify a surface condition and provide an "
    //                   "isomorphism signature.");
    //}

    //// regina::Triangulation<3> tri("caba");
    //// for (const regina::Tetrahedron<3> *tet : tri.tetrahedra()) {
    ////     std::cout << "Tetrahedron " << tet->index() << "\n";
    ////     for (int i = 0; i < 4; ++i) {
    ////         std::cout << "Vertex " << i << " = " << tet->vertex(i)->index()
    ////                   << "\n";
    ////     }
    ////     std::cout << "\n";
    //// }
    //// std::vector<int> edgeIndices = {2, 0, 7, 8, 5};
    //// Curve<3> c;
    //// for (int i : edgeIndices) {
    ////     c.push_back(tri.edge(i));
    //// }
    //// Knot k(tri, c);
    //// std::cout << k << "\n";
    //// k.simplify();
    //// std::cout << k;

    //// std::cout << "Original = " << tri.isoSig() << "\n";
    //// Knot k = {tri, {}};
    //// std::vector<regina::Tetrahedron<3> *> tets;
    //// for (regina::Tetrahedron<3> *tet : k.tri_.tetrahedra()) {
    ////     tets.push_back(tet);
    //// }

    //// k.subdivideSharedVertexSequence_(tets);
    //// std::cout << "Subdivided = " << k.tri_.isoSig() << "\n";
    ////  tri.newSimplex();

    //// regina::Triangulation<2> tri = regina::Example<2>::orientable(6, 0);
    //// tri.newSimplex();
    //// tri.subdivide();
    //// tri.subdivide();
    ////    tri.subdivide();

    // std::cout << "[*] Building gluing graph...\n";
    // SurfaceFinder finder(tri, cond);
    // std::cout << "[+] Total gluing graph nodes = " << finder.countNodes()
    //           << "\n"
    //           << "[+] Total gluing graph edges = " << finder.countEdges()
    //           << "\n";
    ////<< "[+] Graph = \n" << finder << "\n";

    // auto &surfaces = finder.findSurfaces();

    // surfacesDetail(surfaces, cond);

    //// for (const auto *comp : tri.boundaryComponents()) {
    ////     for (const auto *v : comp->vertices()) {
    ////         std::cout << v->index() << " ";
    ////     }
    ////     std::cout << "\n";
    //// }
    //// std::cout << "\n";

    // std::vector<std::tuple<size_t, Link, const KnottedSurface<4> *>>
    // boundaries; for (const auto &surface : surfaces) {
    //     if (surface.isClosed())
    //         continue;

    //    for (auto &[component, link] : surface.boundary()) {
    //        boundaries.emplace_back(component, std::move(link), &surface);
    //    }
    //}

    // std::cout << "[+] Found " << boundaries.size() << " total boundary
    // links\n"; std::cout << "[*] Sorting links by complexity...\n";

    // std::sort(boundaries.begin(), boundaries.end(),
    //           [](const auto &b1, const auto &b2) {
    //               return std::get<1>(b2) < std::get<1>(b1);
    //           });

    // int count = 0;
    // if (cond == BoundaryCondition::boundary) {
    //     for (const auto &[component, l, s] : boundaries) {
    //         if (count > 10)
    //             break;
    //         std::cout << s->detail() << " has boundary " << l
    //                   << " (boundary component " << component << ")\n";
    //         l.recognizeComplement();
    //         ++count;
    //     }
    // }

    return 0;
}

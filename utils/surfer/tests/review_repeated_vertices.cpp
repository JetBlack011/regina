//
//  review_repeated_vertices.cpp
//
//  Review utility (not a CTest test): for every row of one or more PD-code
//  tables, builds the ambient triangulation exactly as verifyslicegenus does
//  (knotbuilder::buildLink, then CobordismBuilder::thicken() with the collar,
//  no cone) and counts
//    - loop edges of T (an edge whose two ends are the same vertex);
//    - loop edges of the thickening;
//    - triangles of the thickening with a repeated vertex, and among them
//      those whose repeated vertex is interior.
//  A triangle with two corners at one interior vertex is where a petal
//  checked before all its corners are registered sees a partial trace
//  (tests/predicate_order_test.cpp); on 2026-09-26 no table row's thickening
//  had one (17,153 rows).
//
//  Usage:
//      review_repeated_vertices [--layers N] table.csv [table.csv ...]
//
//  Each table is "Name,PD Notation,..." with a header row, as in
//  cobordism-atlas/data/.

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <triangulation/dim3.h>
#include <triangulation/dim4.h>

#include "cobordismbuilder.h"
#include "collar.h"
#include "knotbuilder/knotbuilder.h"

namespace {
struct Counts {
    long rows = 0, failed = 0;
    long rowsWithLoopEdgeT = 0, rowsWithLoopEdge4 = 0;
    long rowsWithRepeatedTriangle = 0, rowsWithInteriorRepeat = 0;
    long repeatedTriangles = 0, interiorRepeats = 0;
};
} // namespace

int main(int argc, char *argv[]) {
    int layers = 2;
    int argi = 1;
    if (argi + 1 < argc && std::string(argv[argi]) == "--layers") {
        layers = std::atoi(argv[argi + 1]);
        argi += 2;
    }
    if (argi >= argc) {
        std::cerr << "usage: " << argv[0]
                  << " [--layers N] table.csv [table.csv ...]\n";
        return 2;
    }

    Counts c;
    for (; argi < argc; ++argi) {
        std::ifstream in(argv[argi]);
        if (!in) {
            std::cerr << "cannot open " << argv[argi] << "\n";
            return 2;
        }
        std::string line;
        std::getline(in, line); // header
        while (std::getline(in, line)) {
            if (line.empty())
                continue;
            const size_t c1 = line.find(',');
            const size_t c2 = line.find(',', c1 + 1);
            if (c1 == std::string::npos || c2 == std::string::npos)
                continue;
            const std::string name = line.substr(0, c1);
            const std::string pd = line.substr(c1 + 1, c2 - c1 - 1);
            ++c.rows;
            try {
                auto link = knotbuilder::buildLink(knotbuilder::parsePDCode(pd));
                auto &[t, edges, reversed] = link;

                long loopT = 0;
                for (const regina::Edge<3> *e : t.edges())
                    if (e->vertex(0) == e->vertex(1))
                        ++loopT;

                std::vector<int> edgeIndices;
                for (const regina::Edge<3> *e : edges)
                    edgeIndices.push_back(static_cast<int>(e->index()));
                CobordismBuilder<3> cob(t);
                CollarBuilder collar(edgeIndices);
                for (int i = 0; i < layers; ++i) {
                    cob.thicken();
                    collar.addLayer(cob);
                }
                const regina::Triangulation<4> &tri = cob.getCobordism();

                long loop4 = 0;
                for (const regina::Edge<4> *e : tri.edges())
                    if (e->vertex(0) == e->vertex(1))
                        ++loop4;

                long repeated = 0, interior = 0;
                for (const regina::Triangle<4> *f : tri.triangles()) {
                    const regina::Vertex<4> *v[3] = {f->vertex(0), f->vertex(1),
                                                     f->vertex(2)};
                    for (int a = 0; a < 3; ++a)
                        for (int b = a + 1; b < 3; ++b)
                            if (v[a] == v[b]) {
                                ++repeated;
                                if (!v[a]->isBoundary())
                                    ++interior;
                            }
                }

                if (loopT)
                    ++c.rowsWithLoopEdgeT;
                if (loop4)
                    ++c.rowsWithLoopEdge4;
                if (repeated)
                    ++c.rowsWithRepeatedTriangle;
                if (interior)
                    ++c.rowsWithInteriorRepeat;
                c.repeatedTriangles += repeated;
                c.interiorRepeats += interior;
                if (loopT || loop4 || repeated)
                    std::cout << name << ": loopEdgesT=" << loopT
                              << " loopEdges4=" << loop4
                              << " repeatedTriangles=" << repeated
                              << " interior=" << interior << "\n";
            } catch (const std::exception &e) {
                ++c.failed;
                std::cout << name << ": BUILD FAILED (" << e.what() << ")\n";
            }
        }
    }

    std::cout << "rows " << c.rows << ", build failures " << c.failed
              << ", rows with a loop edge in T " << c.rowsWithLoopEdgeT
              << ", in the thickening " << c.rowsWithLoopEdge4
              << ", rows with a repeated-vertex triangle "
              << c.rowsWithRepeatedTriangle << " (interior "
              << c.rowsWithInteriorRepeat << "), triangles "
              << c.repeatedTriangles << " (interior " << c.interiorRepeats
              << ")\n";
    return c.failed ? 1 : 0;
}

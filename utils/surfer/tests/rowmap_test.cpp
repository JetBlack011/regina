//
//  rowmap_test.cpp
//
//  The row's own link, as verifyslicegenus sets it up on the search side
//  of a real collar-seeded search:
//
//    1. buildRowOrientation(), pinned to the seed's own edges, lands on
//       exactly L x {0} -- every edge, one closed directed curve per
//       component -- even though the diagram's triangulation has
//       automorphisms (D3).
//    2. No searchable face other than the seed touches the search side, so
//       the search side can never change; an unprotected search, for
//       contrast, has plenty that do.
//    3. The bare collar's own boundary classifies as a MATCH. For a link its
//       collar is one annulus per component, each oriented independently;
//       comparing their signs globally (D2) rejected it for some variants.
//
//  Built from real PD codes through knotbuilder, CobordismBuilder and
//  CollarBuilder, exactly as verifyslicegenus does. An optional argument
//  -- a table CSV of Name,PD,... rows -- sweeps every row of it instead
//  (slow; not part of ctest):
//
//    ./rowmap_test ../../../../cobordism-atlas/data/links_..._pd_codes.csv 8
//

#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <unistd.h>
#include <unordered_map>
#include <vector>

#include <triangulation/dim3.h>
#include <triangulation/dim4.h>

#include "../cobordismbuilder.h"
#include "../cobordismgraph.h"
#include "../collar.h"
#include "../embeddedsubmanifold.h"
#include "../knotbuilder.h"
#include "../linkcomplement.h"
#include "../skeleton.h"
#include "../surfacesearch.h"

using namespace cobordismgraph;

static int passed = 0;
static int failed_count = 0;
static int divergedRows = 0; // isIsomorphicTo() would have mapped L differently

#define EXPECT_EQ(actual, expected, desc)                                     \
    do {                                                                      \
        auto _a = (actual);                                                   \
        auto _e = (expected);                                                 \
        if (_a == _e) {                                                       \
            ++passed;                                                         \
        } else {                                                              \
            std::cout << "  FAIL: " << (desc) << "\n        expected " << _e \
                      << ", got " << _a << "\n";                              \
            ++failed_count;                                                   \
        }                                                                     \
    } while (0)

namespace {

// The seed's edges in boundary component `bcIndex`, as sorted indices of
// that component's built triangulation -- the same computation as
// verifyslicegenus's seedEdgesOn().
std::vector<size_t> seedEdgesOn(const regina::Triangulation<4> &tri,
                                const std::vector<int> &seedFaces,
                                size_t bcIndex) {
    const regina::BoundaryComponent<4> *bc = tri.boundaryComponent(bcIndex);
    std::unordered_map<const regina::Edge<4> *, size_t> local;
    for (size_t k = 0; k < bc->countEdges(); ++k)
        local.emplace(bc->edge(k), k);
    std::set<size_t> edges;
    for (int f : seedFaces) {
        const regina::Triangle<4> *t = tri.triangle(f);
        for (int i = 0; i < 3; ++i)
            if (auto it = local.find(t->edge(i)); it != local.end())
                edges.insert(it->second);
    }
    return {edges.begin(), edges.end()};
}

void checkRow(const std::string &name, const std::string &pd) {
    auto link = knotbuilder::buildLink(knotbuilder::parsePDCode(pd));
    auto &[t2, edges2, reversed2] = link;
    const int components = Link(t2, edges2).countComponents();

    std::vector<int> edgeIndices;
    for (const regina::Edge<3> *e : edges2)
        edgeIndices.push_back(static_cast<int>(e->index()));

    CobordismBuilder<3> cob(t2);
    CollarBuilder collar(edgeIndices);
    for (int i = 0; i < 2; ++i) {
        cob.thicken();
        collar.addLayer(cob);
    }
    const size_t bc = cob.baseBoundaryComponent()->index();
    regina::Triangulation<4> tri = cob.getCobordism();
    std::vector<int> seedFaces;
    for (regina::Triangle<4> *t : collar.resolve())
        seedFaces.push_back(static_cast<int>(t->index()));

    const std::vector<size_t> rowEdges = seedEdgesOn(tri, seedFaces, bc);
    EXPECT_EQ(rowEdges.size(), edges2.size(),
              name + ": the seed holds every edge of L on the search side");

    RowOrientation row;
    try {
        row = buildRowOrientation(edges2, reversed2,
                                  tri.boundaryComponent(bc)->build(),
                                  &rowEdges);
    } catch (const regina::InvalidArgument &e) {
        std::cout << "  FAIL: " << name << ": buildRowOrientation threw: "
                  << e.what() << "\n";
        ++failed_count;
        return;
    }
    if (row.divergedFromDefaultIsomorphism) {
        ++divergedRows;
        std::cout << "  (" << name << ": the default isomorphism would have "
                  << "mapped L differently)\n";
    }
    EXPECT_EQ(row.edges == rowEdges, true,
              name + ": the row map lands exactly on the seed's edges");
    EXPECT_EQ(row.components, static_cast<size_t>(components),
              name + ": one closed directed curve per component");

    SurfaceSearch protectedSearch(tri, seedFaces, bc);
    EXPECT_EQ(protectedSearch.countSearchableFacesTouching(bc),
              static_cast<size_t>(0),
              name + ": no searchable non-seed face touches the search side");

    Skeleton<4, 2> skeleton(tri);
    KnottedSurface collarSurface(skeleton, seedFaces);
    std::vector<OrientedCurve> curves;
    for (auto &[c, cs] : collarSurface.orientedBoundaryLinks())
        if (c == bc)
            curves = std::move(cs);
    EXPECT_EQ(curves.size(), static_cast<size_t>(components),
              name + ": the collar's search-side boundary is L");
    EXPECT_EQ(classifyRowOrientation(row, curves,
                                     collarSurface
                                         .boundaryEdgeSurfaceComponent()) ==
                  OrientationVerdict::match,
              true,
              name + ": the bare collar (one annulus per component, each "
                     "oriented on its own) matches the row's orientation");
}

void test_row_map_battery() {
    const std::vector<std::pair<std::string, std::string>> rows = {
        {"3_1", "[[1;5;2;4];[3;1;4;6];[5;3;6;2]]"},
        {"4_1", "[[4;2;5;1];[8;6;1;5];[6;3;7;4];[2;7;3;8]]"},
        {"5_1", "[[2;8;3;7];[4;10;5;9];[6;2;7;1];[8;4;9;3];[10;6;1;5]]"},
        {"8_20", "[[1;7;2;6];[4;13;5;14];[5;9;6;8];[7;3;8;2];[10;15;11;16];"
                 "[12;9;13;10];[14;3;15;4];[16;11;1;12]]"},
        {"L2a1{0}", "PD[X[4; 1; 3; 2]; X[2; 3; 1; 4]]"},
        {"L2a1{1}", "PD[X[4; 2; 3; 1]; X[2; 4; 1; 3]]"},
        {"L4a1{0}", "PD[X[6; 1; 7; 2]; X[8; 3; 5; 4]; X[2; 5; 3; 6]; "
                    "X[4; 7; 1; 8]]"},
        {"L4a1{1}", "PD[X[6; 2; 7; 1]; X[8; 4; 5; 3]; X[2; 8; 3; 7]; "
                    "X[4; 6; 1; 5]]"},
        {"L6a3{0}", "PD[X[8; 1; 9; 2]; X[2; 9; 3; 10]; X[10; 3; 11; 4]; "
                    "X[12; 5; 7; 6]; X[6; 7; 1; 8]; X[4; 11; 5; 12]]"},
        {"L6a3{1}", "PD[X[10; 2; 11; 1]; X[2; 10; 3; 9]; X[8; 4; 9; 3]; "
                    "X[12; 6; 7; 5]; X[6; 12; 1; 11]; X[4; 8; 5; 7]]"},
        {"L6a4{0;0}", "PD[X[6; 1; 7; 2]; X[12; 8; 9; 7]; X[4; 12; 1; 11]; "
                      "X[10; 5; 11; 6]; X[8; 4; 5; 3]; X[2; 9; 3; 10]]"},
        {"L6a4{1;1}", "PD[X[6; 2; 7; 1]; X[12; 6; 9; 5]; X[4; 9; 1; 10]; "
                      "X[10; 7; 11; 8]; X[8; 3; 5; 4]; X[2; 12; 3; 11]]"},
        {"L6n1{1;0}", "PD[X[6; 2; 7; 1]; X[12; 5; 9; 6]; X[4; 12; 1; 11]; "
                      "X[7; 10; 8; 11]; X[3; 5; 4; 8]; X[9; 3; 10; 2]]"},
    };
    for (const auto &[name, pd] : rows)
        checkRow(name, pd);
}

void test_row_map_refuses_foreign_edges() {
    auto link = knotbuilder::buildLink(
        knotbuilder::parsePDCode("[[1;5;2;4];[3;1;4;6];[5;3;6;2]]"));
    auto &[t2, edges2, reversed2] = link;
    CobordismBuilder<3> cob(t2);
    cob.thicken();
    cob.thicken();
    regina::Triangulation<4> tri = cob.getCobordism();
    const size_t bc = cob.baseBoundaryComponent()->index();
    const std::vector<size_t> nonsense = {0};
    bool threw = false;
    try {
        buildRowOrientation(edges2, reversed2,
                            tri.boundaryComponent(bc)->build(), &nonsense);
    } catch (const regina::InvalidArgument &) {
        threw = true;
    }
    EXPECT_EQ(threw, true,
              "no isomorphism takes L onto an arbitrary edge set: refused, "
              "never approximated");

    EmbeddingSearch<4, 2> unprotected(tri);
    EXPECT_EQ(unprotected.countSearchableFacesTouching(bc) > 0, true,
              "without protection, plenty of faces touch the search side "
              "(so the zero asserted for a seeded search means something)");
}

void sweepTable(const std::string &path, int maxCrossings) {
    std::ifstream in(path);
    std::string line;
    std::getline(in, line); // header
    int rows = 0;
    while (std::getline(in, line)) {
        std::string name = line.substr(0, line.find(','));
        std::string rest = line.substr(line.find(',') + 1);
        // The PD field is everything up to the final ",<genus>" column.
        std::string pd = rest.substr(0, rest.rfind(','));
        int crossings = 0;
        for (char c : name) {
            if (c >= '0' && c <= '9')
                crossings = crossings * 10 + (c - '0');
            else if (crossings > 0)
                break;
        }
        if (maxCrossings > 0 && crossings > maxCrossings)
            continue;
        checkRow(name, pd);
        ++rows;
    }
    std::cout << "swept " << rows << " rows of " << path << "\n";
}

} // namespace

int main(int argc, char **argv) {
    if (argc >= 2) {
        sweepTable(argv[1], argc >= 3 ? std::stoi(argv[2]) : 0);
    } else {
        test_row_map_battery();
        test_row_map_refuses_foreign_edges();
    }
    std::cout << passed << " passed, " << failed_count << " failed; "
              << divergedRows << " rows where isIsomorphicTo() would have "
              << "mapped L differently\n";
    return failed_count == 0 ? 0 : 1;
}

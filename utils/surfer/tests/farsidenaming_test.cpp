//
//  farsidenaming_test.cpp
//
//  farside::DiagramNamer on real thickenings, the way verifyslicegenus uses
//  it:
//
//    1. The bare collar L x [0,2] has the row's own link as its far side,
//       so the namer must name it as the row: a table knot by its name, a
//       table link by its base name -- straight from the diagram, with no
//       complement drilled (the fallback counter stays at zero).
//    2. A curve around one triangle of the outgoing boundary is an unknot.
//    3. The signature tables refuse to come back empty.
//

#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <triangulation/dim3.h>
#include <triangulation/dim4.h>

#include "cobordismbuilder.h"
#include "collar.h"
#include "embeddedsubmanifold.h"
#include "farsidenaming.h"
#include "knotbuilder/knotbuilder.h"
#include "linkcomplement.h"
#include "skeleton.h"

static int passed = 0;
static int failed_count = 0;

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

const char *KNOTS = "farsidenaming_test_knots.csv";
const char *LINKS = "farsidenaming_test_links.csv";

void writeTables() {
    std::ofstream k(KNOTS);
    k << "Name,PD Notation,Genus-4D\n"
      << "3_1,[[1;5;2;4];[3;1;4;6];[5;3;6;2]],1\n"
      << "8_20,[[1;7;2;6];[4;13;5;14];[5;9;6;8];[7;3;8;2];[10;15;11;16];"
         "[12;9;13;10];[14;3;15;4];[16;11;1;12]],0\n";
    std::ofstream l(LINKS);
    l << "Name,PD Notation (KnotTheory),Genus-4D\n"
      << "L6a3{0},PD[X[8; 1; 9; 2]; X[2; 9; 3; 10]; X[10; 3; 11; 4]; "
         "X[12; 5; 7; 6]; X[6; 7; 1; 8]; X[4; 11; 5; 12]],1\n"
      << "L6a3{1},PD[X[10; 2; 11; 1]; X[2; 10; 3; 9]; X[8; 4; 9; 3]; "
         "X[12; 6; 7; 5]; X[6; 12; 1; 11]; X[4; 8; 5; 7]],0\n";
}

// The row thickened as verifyslicegenus does it, its collar, and the namer.
struct Row {
    knotbuilder::TriangulationWithLink link;
    CobordismBuilder<3> cob;
    std::vector<int> seed;
    regina::Triangulation<4> tri;

    explicit Row(const std::string &pd)
        : link(knotbuilder::buildLink(knotbuilder::parsePDCode(pd))),
          cob(link.tri) {
        std::vector<int> edgeIndices;
        for (const regina::Edge<3> *e : link.edges)
            edgeIndices.push_back(static_cast<int>(e->index()));
        CollarBuilder collar(edgeIndices);
        for (int i = 0; i < 2; ++i) {
            cob.thicken();
            collar.addLayer(cob);
        }
        tri = cob.getCobordism();
        for (regina::Triangle<4> *t : collar.resolve())
            seed.push_back(static_cast<int>(t->index()));
    }
};

void test_collar_far_side_is_the_row(const std::string &name, const std::string &pd,
                                     const std::string &want,
                                     const farside::SignatureTable &table) {
    Row row(pd);
    farside::DiagramNamer namer(row.link.tri, knotbuilder::parsePDCode(pd).size(),
                                row.cob, table);
    Skeleton<4, 2> skeleton(row.tri);
    KnottedSurface collar(skeleton, row.seed);
    std::string got = "<no far side>";
    for (const auto &[bc, link] : collar.boundaryLinks())
        if (namer.handles(bc)) got = namer.name(link);
    EXPECT_EQ(got, want, name + ": the collar's far side is the row itself");
    EXPECT_EQ(namer.stats().fallbacks.load(), 0LL,
              name + ": named from the diagram, no complement drilled");
}

void test_small_curve_is_unknot(const farside::SignatureTable &table) {
    Row row("[[1;5;2;4];[3;1;4;6];[5;3;6;2]]");
    farside::DiagramNamer namer(row.link.tri, 3, row.cob, table);
    size_t bc = row.tri.boundaryComponent(0)->index() == row.cob.baseBoundaryComponent()->index()
                    ? 1
                    : 0;
    regina::Triangulation<3> boundary = row.tri.boundaryComponent(bc)->build();
    const regina::Triangle<3> *t = boundary.triangle(0);
    std::vector<const regina::Edge<3> *> edges{t->edge(0), t->edge(1), t->edge(2)};
    Link curve(boundary, edges);
    EXPECT_EQ(namer.handles(bc), true, "the namer handles the outgoing boundary");
    EXPECT_EQ(namer.name(curve), std::string("Unknot"),
              "a curve around one triangle is an unknot, from its diagram");
    EXPECT_EQ(namer.stats().unknots.load(), 1LL, "counted as an unknot");
}

void test_tables_refuse_to_be_empty() {
    const char *empty = "farsidenaming_test_empty.csv";
    { std::ofstream e(empty); e << "Name,PD Notation,Genus-4D\n"; }
    bool threw = false;
    try {
        farside::SignatureTable::fromTables(empty, "");
    } catch (const regina::InvalidArgument &) {
        threw = true;
    }
    EXPECT_EQ(threw, true, "a knot table yielding no signatures is an error");
    threw = false;
    try {
        farside::SignatureTable::fromTables("/nonexistent/table.csv", "");
    } catch (const regina::InvalidArgument &) {
        threw = true;
    }
    EXPECT_EQ(threw, true, "a missing table is an error");
    std::remove(empty);
}

} // namespace

int main() {
    writeTables();
    farside::SignatureTable table = farside::SignatureTable::fromTables(KNOTS, LINKS);
    EXPECT_EQ(table.knots(), static_cast<size_t>(2), "both table knots signed");
    EXPECT_EQ(table.links(), static_cast<size_t>(1),
              "both orientations of L6a3 share one unoriented signature");

    test_collar_far_side_is_the_row("8_20", "[[1;7;2;6];[4;13;5;14];[5;9;6;8];[7;3;8;2];"
                                    "[10;15;11;16];[12;9;13;10];[14;3;15;4];[16;11;1;12]]",
                                    "8_20", table);
    test_collar_far_side_is_the_row("3_1", "[[1;5;2;4];[3;1;4;6];[5;3;6;2]]", "3_1", table);
    test_collar_far_side_is_the_row(
        "L6a3{1}", "PD[X[10; 2; 11; 1]; X[2; 10; 3; 9]; X[8; 4; 9; 3]; X[12; 6; 7; 5]; "
                   "X[6; 12; 1; 11]; X[4; 8; 5; 7]]",
        "L6a3", table);
    test_small_curve_is_unknot(table);
    test_tables_refuse_to_be_empty();

    std::remove(KNOTS);
    std::remove(LINKS);
    std::cout << passed << " passed, " << failed_count << " failed\n";
    return failed_count == 0 ? 0 : 1;
}

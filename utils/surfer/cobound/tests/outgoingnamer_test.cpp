//
//  outgoingnamer_test.cpp
//
//  outgoing::OutgoingNamer on real thickenings, the way verifyslicegenus uses
//  it:
//
//    1. The bare collar L x [0,2] has the incoming link itself as its
//       outgoing link, so the namer must name it as the incoming link: a
//       table knot by its name, a table link by its base name -- straight
//       from the diagram, with no complement drilled (the fallback counter
//       stays at zero).
//    2. A curve around one triangle of the outgoing boundary is an unknot.
//    3. The signature tables refuse to come back empty.
//    4. The complement namers share one dispatch, and name an unlink's
//       curves alike with or without the census.
//

#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include <triangulation/dim3.h>
#include <triangulation/dim4.h>

#include "diagramtriangulation/thickening/thickening.h"
#include "surfer/submanifold/submanifold.h"
#include "cobound/outgoing/outgoingnamer.h"
#include "diagramtriangulation/fromdiagram.h"
#include "linknaming/complement/linkcomplement.h"
#include "surfer/submanifold/skeleton.h"

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

const char *KNOTS = "outgoingnamer_test_knots.csv";
const char *LINKS = "outgoingnamer_test_links.csv";

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

// The incoming diagram thickened as a search does it (buildAmbient(): two
// layers, the collar through both).
struct Thickened : ThickenedLink {
    explicit Thickened(const std::string &pd) { buildAmbient(pd, 2, 2, *this); }
};

void test_collar_outgoing_is_the_incoming(const std::string &name, const std::string &pd,
                                     const std::string &want,
                                     const linknaming::SignatureTable &table) {
    Thickened thickened(pd);
    outgoing::OutgoingNamer namer(thickened.link.tri, diagramtriangulation::parsePDCode(pd).size(),
                                *thickened.cob, table);
    Skeleton<4, 2> skeleton(thickened.tri);
    KnottedSurface collar(skeleton, thickened.seedFaces);
    std::string got = "<no far side>";
    for (const auto &[bc, link] : collar.boundaryLinks())
        if (namer.handles(bc)) got = namer.name(link);
    EXPECT_EQ(got, want, name + ": the collar's far side is the row itself");
    EXPECT_EQ(namer.stats().fallbacks.load(), 0LL,
              name + ": named from the diagram, no complement drilled");
}

void test_small_curve_is_unknot(const linknaming::SignatureTable &table) {
    Thickened thickened("[[1;5;2;4];[3;1;4;6];[5;3;6;2]]");
    outgoing::OutgoingNamer namer(thickened.link.tri, 3, *thickened.cob, table);
    size_t bc = thickened.tri.boundaryComponent(0)->index() == thickened.cob->baseBoundaryComponent()->index()
                    ? 1
                    : 0;
    regina::Triangulation<3> boundary = thickened.tri.boundaryComponent(bc)->build();
    const regina::Triangle<3> *t = boundary.triangle(0);
    std::vector<const regina::Edge<3> *> edges{t->edge(0), t->edge(1), t->edge(2)};
    Link curve(boundary, edges);
    EXPECT_EQ(namer.handles(bc), true, "the namer handles the outgoing boundary");
    EXPECT_EQ(namer.name(curve), std::string("Unknot"),
              "a curve around one triangle is an unknot, from its diagram");
    EXPECT_EQ(namer.stats().unknots.load(), 1LL, "counted as an unknot");
}

void test_tables_refuse_to_be_empty() {
    const char *empty = "outgoingnamer_test_empty.csv";
    { std::ofstream e(empty); e << "Name,PD Notation,Genus-4D\n"; }
    bool threw = false;
    try {
        linknaming::SignatureTable::fromTables(empty, "");
    } catch (const regina::InvalidArgument &) {
        threw = true;
    }
    EXPECT_EQ(threw, true, "a knot table yielding no signatures is an error");
    threw = false;
    try {
        linknaming::SignatureTable::fromTables("/nonexistent/table.csv", "");
    } catch (const regina::InvalidArgument &) {
        threw = true;
    }
    EXPECT_EQ(threw, true, "a missing table is an error");
    std::remove(empty);
}

// The one complement dispatch (ComplementBoundaryNamer, phase 3): a lone
// curve by the knot route, several together by the link route, a curve on
// its own by the knot route. Then the two route pairs built on it, census-free
// (UnlinkBoundaryNamer) and census (ComplementNamer), on an unlink's curves:
// the genus and group checks name these before any census is consulted.
std::string knotRoute(const EdgeComplement &) { return "knot"; }
std::string linkRoute(const Link &) { return "link"; }

void test_complement_namers() {
    // The 2-component unlink, as unlinknaming_test's kUnlink2PD.
    const diagramtriangulation::TriangulationWithLink built =
        diagramtriangulation::buildLink({{0, 3, 1, 2}, {1, 3, 0, 2}});
    const Link both(built.tri, built.edges);
    EXPECT_EQ(both.countComponents(), 2, "fixture: two components");
    const Link one(built.tri, both.comps_[0].edges());
    EXPECT_EQ(one.countComponents(), 1, "fixture: one of them alone");

    const ComplementBoundaryNamer routes(knotRoute, linkRoute);
    EXPECT_EQ(routes.nameLink(0, one), std::string("knot"), "a lone curve goes by the knot route");
    EXPECT_EQ(routes.nameLink(0, both), std::string("link"),
              "two curves go together by the link route");
    EXPECT_EQ(routes.nameCurve(0, both.comps_[1]), std::string("knot"),
              "a curve named on its own goes by the knot route");

    const UnlinkBoundaryNamer censusFree{};
    const outgoing::ComplementNamer census{};
    const std::vector<std::pair<std::string, const BoundaryNamer *>> namers = {
        {"census-free", &censusFree}, {"census", &census}};
    for (const auto &[label, namer] : namers) {
        EXPECT_EQ(namer->nameLink(0, both), std::string("2-component unlink"),
                  label + ": the unlink's two curves together");
        EXPECT_EQ(namer->nameLink(0, one), std::string("Unknot"), label + ": one curve alone");
        EXPECT_EQ(namer->nameCurve(0, both.comps_[1]), std::string("Unknot"),
                  label + ": one curve on its own");
    }
}

} // namespace

int main() {
    writeTables();
    linknaming::SignatureTable table = linknaming::SignatureTable::fromTables(KNOTS, LINKS);
    EXPECT_EQ(table.knots(), static_cast<size_t>(2), "both table knots signed");
    EXPECT_EQ(table.links(), static_cast<size_t>(1),
              "both orientations of L6a3 share one unoriented signature");

    test_collar_outgoing_is_the_incoming("8_20", "[[1;7;2;6];[4;13;5;14];[5;9;6;8];[7;3;8;2];"
                                    "[10;15;11;16];[12;9;13;10];[14;3;15;4];[16;11;1;12]]",
                                    "8_20", table);
    test_collar_outgoing_is_the_incoming("3_1", "[[1;5;2;4];[3;1;4;6];[5;3;6;2]]", "3_1", table);
    test_collar_outgoing_is_the_incoming(
        "L6a3{1}", "PD[X[10; 2; 11; 1]; X[2; 10; 3; 9]; X[8; 4; 9; 3]; X[12; 6; 7; 5]; "
                   "X[6; 12; 1; 11]; X[4; 8; 5; 7]]",
        "L6a3", table);
    test_small_curve_is_unknot(table);
    test_tables_refuse_to_be_empty();
    test_complement_namers();

    std::remove(KNOTS);
    std::remove(LINKS);
    std::cout << passed << " passed, " << failed_count << " failed\n";
    return failed_count == 0 ? 0 : 1;
}

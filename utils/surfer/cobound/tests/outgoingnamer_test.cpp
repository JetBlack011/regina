//
//  outgoingnamer_test.cpp
//
//  outgoing::OutgoingNamer on real thickenings, the way a search uses it:
//
//    1. The bare collar L x [0,2] has the incoming link itself as its
//       outgoing link, so the namer must name it as the incoming link: a
//       table knot by its name; a table link, per edge set, by one of its
//       variants (the curves as the search carries them), and per surface,
//       oriented by the collar, as the incoming link's own variant -- all
//       straight from the diagram, with no complement drilled.
//    2. A curve around one triangle of the outgoing boundary is an unknot.
//    3. The complement namers share one dispatch, and name an unlink's
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

#include "cobound/outgoing/outgoingnamer.h"
#include "cobound/search/incoming.h"
#include "cobound/search/preconditions.h"
#include "diagramtriangulation/fromdiagram.h"
#include "diagramtriangulation/thickening/thickening.h"
#include "linknaming/complement/linkcomplement.h"
#include "linknaming/linknamer.h"
#include "linknaming/names.h"
#include "linknaming/tables.h"
#include "surfer/submanifold/skeleton.h"
#include "surfer/submanifold/submanifold.h"

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

void test_collar_outgoing_is_the_incoming(const std::string &name, const std::string &pd,
                                          const linknaming::Tables &tables) {
    // The incoming diagram thickened as a search does it (two layers, the
    // collar through both), with its incoming map.
    search::IncomingThickening thickened;
    search::buildIncoming(pd, 2, 2, thickened);
    outgoing::OutgoingNamer namer(thickened.link.tri, diagramtriangulation::parsePDCode(pd).size(),
                                  *thickened.cob, tables);
    Skeleton<4, 2> skeleton(thickened.tri);
    KnottedSurface collar(skeleton, thickened.seedFaces);
    const linknaming::LinkNamer canonical(tables);
    const std::string want = canonical.canonicalName(*tables.entry(name));
    std::string perEdgeSet = "<no outgoing link>";
    for (const auto &[bc, link] : collar.boundaryLinks())
        if (namer.handles(bc)) perEdgeSet = namer.name(link);
    if (thickened.componentCount == 1) {
        EXPECT_EQ(perEdgeSet, want, name + ": the collar's outgoing knot is the incoming knot");
    } else {
        EXPECT_EQ(perEdgeSet.rfind(linknaming::baseName(name) + "{", 0) == 0, true,
                  name + ": per edge set, a variant of the incoming link (" + perEdgeSet + ")");
        // Per surface, oriented by the collar against the incoming link: the
        // incoming link's own variant.
        const auto oriented = collar.orientedBoundaryLinks();
        const auto surfaceOf = collar.boundaryEdgeSurfaceComponent();
        std::vector<OrientedCurve> incomingCurves;
        for (const auto &[bc, curves] : oriented)
            if (bc == thickened.incomingBC) incomingCurves = curves;
        const search::IncomingOrientationJudgement judged =
            search::judgeIncomingOrientation(*thickened.orientation, incomingCurves, surfaceOf);
        std::string perSurface = "<no outgoing link>";
        for (const auto &[bc, curves] : oriented)
            if (namer.handles(bc)) perSurface = namer.orientedName(curves, surfaceOf, judged.flips);
        EXPECT_EQ(perSurface, want, name + ": per surface, oriented, the incoming link's variant");
    }
    EXPECT_EQ(namer.stats().fallbacks.load(), 0LL,
              name + ": named from the diagram, no complement drilled");
}

void test_small_curve_is_unknot(const linknaming::Tables &tables) {
    search::IncomingThickening thickened;
    search::buildIncoming("[[1;5;2;4];[3;1;4;6];[5;3;6;2]]", 2, 2, thickened);
    outgoing::OutgoingNamer namer(thickened.link.tri, 3, *thickened.cob, tables);
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
    const linknaming::Tables tables = linknaming::Tables::load(KNOTS, LINKS, "");
    test_collar_outgoing_is_the_incoming("8_20", "[[1;7;2;6];[4;13;5;14];[5;9;6;8];[7;3;8;2];"
                                                 "[10;15;11;16];[12;9;13;10];[14;3;15;4];[16;11;1;12]]",
                                         tables);
    test_collar_outgoing_is_the_incoming("3_1", "[[1;5;2;4];[3;1;4;6];[5;3;6;2]]", tables);
    test_collar_outgoing_is_the_incoming(
        "L6a3{1}", "PD[X[10; 2; 11; 1]; X[2; 10; 3; 9]; X[8; 4; 9; 3]; X[12; 6; 7; 5]; "
                   "X[6; 12; 1; 11]; X[4; 8; 5; 7]]",
        tables);
    test_collar_outgoing_is_the_incoming(
        "L6a3{0}", "PD[X[8; 1; 9; 2]; X[2; 9; 3; 10]; X[10; 3; 11; 4]; X[12; 5; 7; 6]; "
                   "X[6; 7; 1; 8]; X[4; 11; 5; 12]]",
        tables);
    test_small_curve_is_unknot(tables);
    test_complement_namers();

    std::remove(KNOTS);
    std::remove(LINKS);
    std::cout << passed << " passed, " << failed_count << " failed\n";
    return failed_count == 0 ? 0 : 1;
}

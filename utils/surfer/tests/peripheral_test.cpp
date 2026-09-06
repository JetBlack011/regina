// peripheral_test.cpp
//
// Tests for peripheral.h: drilling a link while retaining its meridians, and
// reading those meridians as slopes in SnapPea's peripheral basis.
//
// The ground truth is the one case where the whole loop closes without
// leaving C++: pinching a loop edge of Example3::lens(1,0) gives a 1-cusped
// ideal triangulation of the unknot exterior whose meridian is (1,-1) in
// SnapPea's basis, up to sign. That case has *no finite vertices*, which is
// why regina::SnapPeaTriangulation may be used on it here and nowhere else --
// see peripheral.h's warning, and the explicit assertion below that keeps the
// two situations from being confused.
//
// Realistic inputs (a knot built by knotbuilder, drilled out of a triangulated
// S^3 that keeps its other vertices) cannot get a SnapPea basis in-process for
// exactly that reason, so for those we check the meridian is a structurally
// valid closed curve on the right cusp, and leave the slope round trip to
// tools/identify_far_sides.py.

#include <iostream>
#include <map>
#include <string>
#include <unistd.h>
#include <vector>

#include <snappea/snappeatriangulation.h>
#include <triangulation/dim3.h>
#include <triangulation/example3.h>

#include "../knotbuilder.h"
#include "../linkcomplement.h"
#include "../peripheral.h"

static int passed = 0, failed_count = 0;

namespace {
bool colorEnabled() {
    static bool enabled = isatty(fileno(stdout));
    return enabled;
}
std::ostream &green(std::ostream &os) {
    return colorEnabled() ? os << "\033[32m" : os;
}
std::ostream &red(std::ostream &os) {
    return colorEnabled() ? os << "\033[31m" : os;
}
std::ostream &resetColor(std::ostream &os) {
    return colorEnabled() ? os << "\033[0m" : os;
}
} // namespace

#define EXPECT_EQ(actual, expected, desc)                                      \
    do {                                                                       \
        auto _a = (actual);                                                    \
        auto _e = (expected);                                                  \
        if (_a == _e) {                                                        \
            std::cout << green << "  PASS: " << resetColor << (desc) << "\n";  \
            ++passed;                                                          \
        } else {                                                               \
            std::cout << red << "  FAIL: " << (desc) << "\n"                   \
                      << "        expected " << _e << ", got " << _a           \
                      << resetColor << "\n";                                   \
            ++failed_count;                                                    \
        }                                                                      \
    } while (0)

#define EXPECT_TRUE(actual, desc) EXPECT_EQ(bool(actual), true, desc)

namespace {

/**
 * Whether `curve` is a closed curve on a cusp cross-section of `tri`: every
 * corner it meets has its crossings summing to zero (what goes in comes out),
 * and every crossing is matched by the opposite crossing on the other side of
 * the glued face (leaving one triangle is entering the next).
 *
 * This is the strongest statement available without a hyperbolic structure,
 * and it is what catches a mistabulated pinch annulus.
 */
bool isClosedCurve(const regina::Triangulation<3> &tri,
                   const peripheral::Curve &curve) {
    std::map<std::tuple<size_t, int, int>, int> value;
    for (const peripheral::Crossing &c : curve)
        value[{c.tet, c.vertex, c.face}] += c.sign;

    std::map<std::pair<size_t, int>, int> cornerSum;
    for (const auto &[key, v] : value) {
        const auto &[t, vert, face] = key;
        cornerSum[{t, vert}] += v;
    }
    for (const auto &[corner, sum] : cornerSum)
        if (sum != 0)
            return false;

    for (const auto &[key, v] : value) {
        const auto &[t, vert, face] = key;
        if (v == 0)
            continue;
        regina::Tetrahedron<3> *tet = tri.tetrahedron(t);
        regina::Tetrahedron<3> *adj = tet->adjacentTetrahedron(face);
        if (!adj)
            return false;
        const regina::Perm<4> g = tet->adjacentGluing(face);
        auto it = value.find({adj->index(), g[vert], g[face]});
        if (it == value.end() || it->second != -v)
            return false;
    }
    return true;
}

/** Every regina Vertex the curve's crossings sit at. */
std::set<const regina::Vertex<3> *>
verticesMet(const regina::Triangulation<3> &tri,
            const peripheral::Curve &curve) {
    std::set<const regina::Vertex<3> *> out;
    for (const peripheral::Crossing &c : curve)
        out.insert(tri.tetrahedron(c.tet)->vertex(c.vertex));
    return out;
}

/** The drilled unknot exterior used as ground truth, plus its meridian. */
peripheral::DrilledWithMeridians groundTruth() {
    regina::Triangulation<3> lens = regina::Example<3>::lens(1, 0);
    std::vector<const regina::Edge<3> *> loop;
    for (const regina::Edge<3> *e : lens.edges())
        if (e->vertex(0) == e->vertex(1)) {
            loop.push_back(e);
            break;
        }
    return peripheral::drillWithMeridians(lens, {loop});
}

void testCompleteBasis() {
    std::cout << "completeBasis()\n";
    for (auto [a, b] : std::vector<std::pair<long, long>>{
             {1, -1}, {-1, 1}, {1, 0}, {0, 1}, {3, 2}, {-5, 7}, {2, -3}}) {
        auto [c, d] = peripheral::completeBasis(a, b);
        EXPECT_EQ(a * d - b * c, 1L,
                  "det [(" + std::to_string(a) + "," + std::to_string(b) +
                      "),(" + std::to_string(c) + "," + std::to_string(d) +
                      ")] == 1");
    }
}

void testSnapPeaOriented() {
    std::cout << "snapPeaOriented()\n";
    auto drilled = groundTruth();
    const std::string text = peripheral::snapPeaOriented(drilled.tri);
    EXPECT_TRUE(text.find("oriented_manifold") != std::string::npos,
                "declares oriented_manifold");
    EXPECT_TRUE(text.find("unknown_orientability") == std::string::npos,
                "no longer declares unknown_orientability");
    // Everything else must be Regina's own output, byte for byte: the header
    // patch is the only intended difference.
    std::string regina = drilled.tri.snapPea();
    const size_t at = regina.find("unknown_orientability");
    regina.replace(at, std::string("unknown_orientability").size(),
                   "oriented_manifold");
    EXPECT_EQ(text, regina, "otherwise identical to Triangulation<3>::snapPea()");
}

void testGroundTruth() {
    std::cout << "ground truth: unknot exterior from a pinched loop edge\n";
    auto drilled = groundTruth();
    const regina::Triangulation<3> &tri = drilled.tri;

    EXPECT_EQ(drilled.meridians.size(), size_t(1), "one meridian");
    EXPECT_TRUE(tri.isIdeal(), "drilled complement is ideal");
    EXPECT_TRUE(tri.isOriented(), "drilled complement is oriented");
    EXPECT_EQ(tri.size(), regina::Example<3>::lens(1, 0).size() + 2,
              "two tetrahedra inserted by the pinch");
    EXPECT_TRUE(isClosedCurve(tri, drilled.meridians[0]),
                "meridian is a closed curve on the cusp");

    // This case, and only this case, has no finite vertices, so
    // regina::SnapPeaTriangulation will not retriangulate it. Assert that
    // rather than assume it: if it ever stops holding, the slope below would
    // silently describe a different manifold.
    EXPECT_TRUE(tri.countVertices() <= tri.countBoundaryComponents(),
                "no finite vertices, so SnapPea will not retriangulate");

    regina::SnapPeaTriangulation snappea(tri);
    EXPECT_TRUE(!snappea.isNull(), "SnapPea accepts it");
    EXPECT_EQ(snappea.size(), tri.size(), "SnapPea preserved the triangulation");

    peripheral::SnapPeaFile file = peripheral::parseSnapPea(snappea.snapPea());
    EXPECT_EQ(file.tets.size(), tri.size(), "parsed every tetrahedron");
    EXPECT_EQ(file.numOrCusps + file.numNonOrCusps, 1, "parsed one cusp");

    // SnapPea's own basis must behave like a basis under our intersection
    // pairing. This is the check that catches a wrong remaining_face table or
    // a wrong FLOW, independently of anything to do with meridians.
    const long mm = peripheral::intersectionNumber(file, 0, file.meridian,
                                                   file.meridian);
    const long ll = peripheral::intersectionNumber(file, 0, file.longitude,
                                                   file.longitude);
    const long ml = peripheral::intersectionNumber(file, 0, file.meridian,
                                                   file.longitude);
    EXPECT_EQ(mm, 0L, "<m, m> == 0");
    EXPECT_EQ(ll, 0L, "<l, l> == 0");
    EXPECT_TRUE(ml == 1 || ml == -1, "<m, l> == +-1");
    // std::pair is not streamable, and EXPECT_EQ prints what it compares, so
    // compare the coefficients rather than the pairs.
    auto [mA, mB] = peripheral::slope(file, 0, file.meridian);
    EXPECT_EQ(mA, 1L, "slope of SnapPea's m has a == 1");
    EXPECT_EQ(mB, 0L, "slope of SnapPea's m has b == 0");
    auto [lA, lB] = peripheral::slope(file, 0, file.longitude);
    EXPECT_EQ(lA, 0L, "slope of SnapPea's l has a == 0");
    EXPECT_EQ(lB, 1L, "slope of SnapPea's l has b == 1");

    auto [ourA, ourB] = peripheral::slope(
        file, 0, peripheral::toField(drilled.meridians[0], file.tets.size()));
    const long a = ourA, b = ourB;
    EXPECT_TRUE((a == 1 && b == -1) || (a == -1 && b == 1),
                "our meridian is (1, -1) up to sign [got (" +
                    std::to_string(a) + ", " + std::to_string(b) + ")]");
}

/** Drills the link `pd` builds and runs the structural checks on it. */
void testDiagram(const std::string &label, const std::string &pd,
                 size_t expectedComponents) {
    std::cout << label << "\n";
    knotbuilder::TriangulationWithLink built =
        knotbuilder::buildLink(knotbuilder::parsePDCode(pd));
    Link link(built.tri, built.edges);
    EXPECT_EQ(size_t(link.countComponents()), expectedComponents,
              "component count");

    const size_t before = built.tri.size();
    peripheral::DrilledWithMeridians drilled =
        link.buildComplementWithPeripheral();

    EXPECT_EQ(drilled.tri.size(), before + 2 * built.edges.size(),
              "two tetrahedra inserted per drilled edge");
    EXPECT_TRUE(drilled.tri.isOriented(), "drilled complement is oriented");
    EXPECT_TRUE(drilled.tri.isIdeal(), "drilled complement is ideal");
    EXPECT_EQ(drilled.meridians.size(), expectedComponents,
              "one meridian per component");

    // The realistic case: the ambient S^3 triangulation's other vertices
    // survive as finite vertices, which is exactly why this input cannot be
    // handed to regina::SnapPeaTriangulation.
    EXPECT_TRUE(drilled.tri.countVertices() >
                    drilled.tri.countBoundaryComponents(),
                "finite vertices present (so SnapPea would retriangulate)");

    std::set<const regina::Vertex<3> *> cusps;
    for (size_t i = 0; i < drilled.meridians.size(); ++i) {
        const peripheral::Curve &mu = drilled.meridians[i];
        EXPECT_TRUE(isClosedCurve(drilled.tri, mu),
                    "meridian " + std::to_string(i) + " is a closed curve");
        auto met = verticesMet(drilled.tri, mu);
        EXPECT_EQ(met.size(), size_t(1),
                  "meridian " + std::to_string(i) + " lies on one vertex");
        if (met.size() == 1) {
            const regina::Vertex<3> *v = *met.begin();
            EXPECT_TRUE(v->isIdeal(), "meridian " + std::to_string(i) +
                                          " lies on an ideal vertex");
            cusps.insert(v);
        }
    }
    EXPECT_EQ(cusps.size(), expectedComponents,
              "distinct components got distinct cusps");
}

/**
 * Reversing a component's traversal negates its meridian, and touches nothing
 * else.
 *
 * This is the whole content of the signed overload: a meridian is oriented by
 * `lk(mu, K) = +1`, so reversing K reverses mu. Without it every meridian
 * carries an arbitrary sign fixed by Regina's edge numbering, which is enough
 * to recognise an unoriented link and never enough to recognise an oriented
 * one.
 */
void testDirectionNegatesMeridian() {
    std::cout << "reversing a component negates its meridian\n";
    knotbuilder::TriangulationWithLink built = knotbuilder::buildLink(
        knotbuilder::parsePDCode("PD[X[4; 1; 3; 2]; X[2; 3; 1; 4]]"));
    Link link(built.tri, built.edges);

    std::vector<std::vector<peripheral::DirectedEdge>> forward;
    for (int c = 0; c < link.countComponents(); ++c) {
        std::vector<peripheral::DirectedEdge> edges;
        for (const regina::Edge<3> *e : link.comps_[c].edges())
            edges.push_back({e, false});
        forward.push_back(std::move(edges));
    }
    auto flipFirst = forward;
    for (auto &e : flipFirst[0])
        e.reversed = true;

    peripheral::DrilledWithMeridians a =
        peripheral::drillWithMeridians(built.tri, forward);
    peripheral::DrilledWithMeridians b =
        peripheral::drillWithMeridians(built.tri, flipFirst);

    EXPECT_EQ(a.meridians.size(), b.meridians.size(), "same component count");
    // Same drilling order either way, so the triangulations agree and the
    // curves are directly comparable crossing by crossing.
    EXPECT_EQ(a.tri.isoSig(), b.tri.isoSig(), "same drilled triangulation");

    bool negated = a.meridians[0].size() == b.meridians[0].size();
    for (size_t i = 0; negated && i < a.meridians[0].size(); ++i) {
        const peripheral::Crossing &x = a.meridians[0][i];
        const peripheral::Crossing &y = b.meridians[0][i];
        negated = x.tet == y.tet && x.vertex == y.vertex &&
                  x.face == y.face && x.sign == -y.sign;
    }
    EXPECT_TRUE(negated, "reversed component 0 has exactly the negated meridian");

    bool untouched = a.meridians[1].size() == b.meridians[1].size();
    for (size_t i = 0; untouched && i < a.meridians[1].size(); ++i) {
        const peripheral::Crossing &x = a.meridians[1][i];
        const peripheral::Crossing &y = b.meridians[1][i];
        untouched = x.tet == y.tet && x.vertex == y.vertex &&
                    x.face == y.face && x.sign == y.sign;
    }
    EXPECT_TRUE(untouched, "component 1 is untouched");
}

/**
 * The undirected overload is the directed one with every edge taken forwards.
 *
 * Kept because the undirected overload is still the right call wherever the
 * direction genuinely is not known, and it must not quietly drift away from
 * the directed one.
 */
void testUndirectedMatchesForward() {
    std::cout << "undirected drilling == every edge forwards\n";
    knotbuilder::TriangulationWithLink built = knotbuilder::buildLink(
        knotbuilder::parsePDCode("[[1;5;2;4];[3;1;4;6];[5;3;6;2]]"));
    Link link(built.tri, built.edges);

    std::vector<std::vector<const regina::Edge<3> *>> plain;
    std::vector<std::vector<peripheral::DirectedEdge>> forward;
    for (int c = 0; c < link.countComponents(); ++c) {
        std::vector<const regina::Edge<3> *> bare;
        std::vector<peripheral::DirectedEdge> directed;
        for (const regina::Edge<3> *e : link.comps_[c].edges()) {
            bare.push_back(e);
            directed.push_back({e, false});
        }
        plain.push_back(std::move(bare));
        forward.push_back(std::move(directed));
    }

    peripheral::DrilledWithMeridians a =
        peripheral::drillWithMeridians(built.tri, plain);
    peripheral::DrilledWithMeridians b =
        peripheral::drillWithMeridians(built.tri, forward);
    EXPECT_EQ(a.tri.isoSig(), b.tri.isoSig(), "same drilled triangulation");
    bool same = a.meridians.size() == b.meridians.size();
    for (size_t c = 0; same && c < a.meridians.size(); ++c) {
        same = a.meridians[c].size() == b.meridians[c].size();
        for (size_t i = 0; same && i < a.meridians[c].size(); ++i)
            same = a.meridians[c][i].sign == b.meridians[c][i].sign &&
                   a.meridians[c][i].tet == b.meridians[c][i].tet;
    }
    EXPECT_TRUE(same, "identical meridians");
}

/**
 * orient()'s vertex 2/3 swap reverses local edge {2,3}, and only that edge.
 *
 * This is the correction the directed overload folds in, and it is not
 * cosmetic: on real witness data it changes the recorded sign for 2,710 of the
 * meridian lines in a 40-witness dump. The old undirected code left the sign
 * at +1 there, which was harmless only because its consumer
 * (Isometry::extends_to_link()) tests meridian against +-meridian per cusp and
 * so never looked at the sign.
 */
void testSwap23ReversesOnlyEdge23() {
    std::cout << "orient()'s 2/3 swap reverses exactly local edge {2,3}\n";
    auto swap23 = [](int v) { return v == 2 ? 3 : (v == 3 ? 2 : v); };
    for (int e = 0; e < 6; ++e) {
        const int a = regina::Edge<3>::edgeVertex[e][0];
        const int b = regina::Edge<3>::edgeVertex[e][1];
        const bool reversed = swap23(a) > swap23(b);
        EXPECT_EQ(reversed, e == 5,
                  "local edge " + std::to_string(e) + " {" +
                      std::to_string(a) + "," + std::to_string(b) + "}" +
                      (e == 5 ? " reverses" : " keeps its direction"));
    }
}

} // namespace

int main() {
    testCompleteBasis();
    testSnapPeaOriented();
    testGroundTruth();
    testDiagram("trefoil 3_1 built by knotbuilder",
                "[[1;5;2;4];[3;1;4;6];[5;3;6;2]]", 1);
    testDiagram("Hopf link L2a1{0} built by knotbuilder",
                "PD[X[4; 1; 3; 2]; X[2; 3; 1; 4]]", 2);
    testDirectionNegatesMeridian();
    testUndirectedMatchesForward();
    testSwap23ReversesOnlyEdge23();

    std::cout << "\n" << passed << " passed, " << failed_count << " failed\n";
    return failed_count == 0 ? 0 : 1;
}
